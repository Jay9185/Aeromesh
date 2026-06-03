import argparse
import xml.etree.ElementTree as ET
import requests
import numpy as np
import trimesh
import osmnx as ox
import geopandas as gpd
from shapely.geometry import Polygon, Point
from scipy.interpolate import RegularGridInterpolator
import pandas as pd
import warnings
import gc
import os

warnings.filterwarnings('ignore')

def automate_flight_to_3d(kml_path, api_key, target_width_mm=150, padding_deg=0.02, z_scale_exaggeration=2.5, max_points=800):
    if not api_key:
        raise ValueError("CRITICAL: You must provide an OpenTopography API key.")

    print(f"Step 1: Parsing KML ({kml_path}) and auto-detecting bounding box...")
    tree = ET.parse(kml_path)
    root = tree.getroot()
    raw_coords = []
    
    for elem in root.iter():
        if 'coord' in elem.tag.lower() and elem.text:
            parts = elem.text.strip().replace(',', ' ').split()
            for i in range(0, len(parts) - 2, 3):
                try:
                    lon, lat, alt_m = float(parts[i]), float(parts[i+1]), float(parts[i+2])
                    raw_coords.append([lon, lat, alt_m])
                except ValueError:
                    continue

    if not raw_coords:
        raise ValueError("CRITICAL: No valid coordinates found.")
        
    if len(raw_coords) > max_points:
        step = len(raw_coords) // max_points
        raw_coords = raw_coords[::step]
        print(f"  -> Downsampled to {len(raw_coords)} points for optimization.")
            
    flight_path = np.array(raw_coords)
    lons, lats, alts_m = flight_path[:, 0], flight_path[:, 1], flight_path[:, 2]
    
    lon_min, lon_max = lons.min() - padding_deg, lons.max() + padding_deg
    lat_min, lat_max = lats.min() - padding_deg, lats.max() + padding_deg

    print("\nStep 2: Fetching Real Earth DEM from OpenTopography (SRTMGL3)...")
    url = "https://portal.opentopography.org/API/globaldem"
    params = {
        "demtype": "SRTMGL3", 
        "south": lat_min, "north": lat_max,
        "west": lon_min, "east": lon_max,
        "outputFormat": "AAIGrid",
        "API_Key": api_key
    }
    
    response = requests.get(url, params=params)
    if response.status_code != 200 or "Error" in response.text:
        raise ValueError(f"API Failed. Response: {response.text[:100]}")
        
    lines = response.text.splitlines()
    ncols, nrows = int(lines[0].split()[1]), int(lines[1].split()[1])
    xllcorner, yllcorner = float(lines[2].split()[1]), float(lines[3].split()[1])
    cellsize = float(lines[4].split()[1])
    
    elev_data = [float(v) for line in lines[6:] for v in line.split()]
    elevation_matrix_asc = np.flipud(np.array(elev_data).reshape((nrows, ncols)))
    grid_lats_asc = yllcorner + cellsize * np.arange(nrows)
    grid_lons = xllcorner + cellsize * np.arange(ncols)
    
    del response, lines, elev_data
    gc.collect()

    print("\nStep 3: Calculating geographic scale ratios...")
    lat_dist = (lat_max - lat_min) * 111000  
    lon_dist = (lon_max - lon_min) * 111000 * np.cos(np.radians((lat_min + lat_max) / 2))
    scale_x = target_width_mm / (lon_max - lon_min)
    scale_y = (target_width_mm * (lat_dist / lon_dist)) / (lat_max - lat_min)
    scale_z = (target_width_mm / lon_dist) * z_scale_exaggeration

    print("\nStep 4: Generating 3D Mountain Terrain...")
    vertices, faces = [], []
    for i in range(nrows):
        for j in range(ncols):
            z = max(elevation_matrix_asc[i, j], 0) * scale_z
            vertices.append([(grid_lons[j] - lon_min) * scale_x, (grid_lats_asc[i] - lat_min) * scale_y, z])
            
    for i in range(nrows - 1):
        for j in range(ncols - 1):
            v0 = i * ncols + j
            v1, v2, v3 = v0 + 1, (i + 1) * ncols + j, (i + 1) * ncols + j + 1
            faces.extend([[v0, v1, v2], [v1, v3, v2]])
            
    terrain_mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
    del vertices, faces
    gc.collect()

    print("\nStep 5: Interpolating ground altitudes...")
    interp = RegularGridInterpolator((grid_lats_asc, grid_lons), elevation_matrix_asc, bounds_error=False, fill_value=0)
    ground_z = interp(np.column_stack((lats, lons))) * scale_z
    path_points = np.column_stack(((lons - lon_min) * scale_x, (lats - lat_min) * scale_y, alts_m * scale_z))

    print("\nStep 6: Fetching OpenStreetMap Building Data...")
    building_meshes = []
    try:
        tags = {"building": True, "aeroway": ["hangar", "terminal"], "man_made": True}
        try:
            gdf = ox.features_from_bbox(bbox=(lon_min, lat_min, lon_max, lat_max), tags=tags)
        except TypeError:
            gdf = ox.features_from_bbox(north=lat_max, south=lat_min, east=lon_max, west=lon_min, tags=tags)
            
        if not gdf.empty:
            for _, row in gdf.iterrows():
                geom = row.geometry
                polygons = [geom] if isinstance(geom, Polygon) else getattr(geom, 'geoms', [])
                
                for poly in polygons:
                    if not isinstance(poly, Polygon): continue
                    x_coords, y_coords = poly.exterior.coords.xy
                    pts = np.column_stack(((np.array(x_coords) - lon_min) * scale_x, (np.array(y_coords) - lat_min) * scale_y))
                    b_ground_z = interp([[poly.centroid.y, poly.centroid.x]])[0] * scale_z
                    
                    height_m = 8 
                    if 'height' in row and pd.notnull(row['height']):
                        try: height_m = float(''.join(c for c in str(row['height']) if c.isdigit() or c == '.'))
                        except: pass
                    elif 'building:levels' in row and pd.notnull(row['building:levels']):
                        try: height_m = float(row['building:levels']) * 3.5
                        except: pass
                    
                    try:
                        b_mesh = trimesh.creation.extrude_polygon(Polygon(pts), height=height_m * scale_z)
                        b_mesh.apply_translation([0, 0, b_ground_z])
                        building_meshes.append(b_mesh)
                    except: continue
            
            if building_meshes:
                chunk_size = 250
                city_chunks = []
                for i in range(0, len(building_meshes), chunk_size):
                    city_chunks.append(trimesh.util.concatenate(building_meshes[i:i+chunk_size]))
                city_mesh = trimesh.util.concatenate(city_chunks)
                del building_meshes, city_chunks, gdf
                gc.collect()
            else:
                city_mesh = trimesh.Trimesh()
        else:
            city_mesh = trimesh.Trimesh()
    except Exception as e:
        print(f"  -> OSM fetch failed: {e}")
        city_mesh = trimesh.Trimesh()

    print("\nStep 7: Engineering Pylons and Track Tube...")
    pylon_meshes = []
    step = max(1, len(path_points) // 40) 
    
    for i in range(0, len(path_points), step):
        p_top = path_points[i]
        height = p_top[2] - ground_z[i]
        if height > 0:
            pylon = trimesh.creation.cylinder(radius=0.8, height=height)
            pylon.apply_translation([p_top[0], p_top[1], ground_z[i] + (height / 2)])
            pylon_meshes.append(pylon)
            
    all_pylons = trimesh.util.concatenate(pylon_meshes) if pylon_meshes else trimesh.Trimesh()

    tube_radius = 0.6  
    circle_profile = Point(0, 0).buffer(tube_radius)
    
    try:
        path_tube = trimesh.creation.sweep_polygon(polygon=circle_profile, path=path_points)
    except Exception as e:
        path_tube = trimesh.util.concatenate([trimesh.creation.cylinder(radius=tube_radius, segment=[path_points[i], path_points[i+1]]) for i in range(len(path_points)-1)])

    print("\nStep 8: Assembling 3MF File with Baked Colors...")
    terrain_mesh.metadata['name'] = '1_Terrain'
    terrain_mesh.visual.face_colors = [130, 130, 130, 255]
    
    if not city_mesh.is_empty:
        city_mesh.metadata['name'] = '2_Buildings'
        city_mesh.visual.face_colors = [240, 240, 240, 255]
    else:
        city_mesh = trimesh.creation.box([0.1, 0.1, 0.1]); city_mesh.metadata['name'] = '2_Buildings_EMPTY'
        
    if not all_pylons.is_empty:
        all_pylons.metadata['name'] = '3_Pylons'
        all_pylons.visual.face_colors = [30, 30, 30, 255]
    else:
        all_pylons = trimesh.creation.box([0.1, 0.1, 0.1]); all_pylons.metadata['name'] = '3_Pylons_EMPTY'
        
    path_tube.metadata['name'] = '4_FlightPath'
    path_tube.visual.face_colors = [255, 100, 0, 255]

    master_scene = trimesh.Scene([terrain_mesh, city_mesh, all_pylons, path_tube])
    master_filename = kml_path.replace('.kml', '_Assembly.3mf').split('/')[-1]
    master_scene.export(master_filename)
    print(f"\n[SUCCESS] Saved to {master_filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert ForeFlight KML to 3D Multi-Color Assembly")
    parser.add_argument("kml_file", help="Path to the .kml track log")
    parser.add_argument("--api-key", help="OpenTopography API Key (or set OPENTOPO_API_KEY env var)")
    
    args = parser.parse_args()
    api_key = args.api_key or os.getenv("OPENTOPO_API_KEY")
    
    automate_flight_to_3d(args.kml_file, api_key)
