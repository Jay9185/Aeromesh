"""
AeroMesh: Flight Telemetry to 3D Assembly
------------------------------------------
Converts a ForeFlight KML track log into a multi-color, Bambu-Lab-ready
.3mf assembly: real terrain (DEM), OSM building extrusions, engineered
support pylons, the flight path itself, and a printable base plate.
"""

import argparse
import os
import sys
import time
import warnings
from io import StringIO

import numpy as np
import pandas as pd
import requests
import trimesh
from scipy.interpolate import RegularGridInterpolator
from shapely.geometry import Point, Polygon

# Prefer defusedxml (hardened against XML entity-expansion attacks) if
# available, since KML files may come from an untrusted source.
try:
    import defusedxml.ElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
    print("WARNING: 'defusedxml' not installed - parsing KML with the "
          "standard library's XML parser, which is not hardened against "
          "malicious XML. Run 'pip install defusedxml' to fix this.")

# Only silence the specific warning categories that osmnx/geopandas are
# known to spam, rather than hiding every warning in the process.
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning, module='osmnx')

import osmnx as ox  # noqa: E402 (import after warning filters are set)
import geopandas as gpd  # noqa: E402,F401

# Bambu Lab printable volumes (mm), for the build-plate fit check.
BAMBU_BED_SIZES = {
    "a1_mini": (180, 180, 180),
    "a1": (256, 256, 256),
    "p1p": (256, 256, 256),
    "p1s": (256, 256, 256),
    "x1c": (256, 256, 256),
}
NODATA_THRESHOLD = -1000  # SRTM NODATA sentinel is typically -9999


# --------------------------------------------------------------------------
# Step 1: KML parsing
# --------------------------------------------------------------------------
def parse_kml_coordinates(kml_path, max_points):
    """Extract [lon, lat, alt_m] points from a ForeFlight KML, downsampled
    to at most max_points while always keeping the first and last point."""
    tree = ET.parse(kml_path)
    root = tree.getroot()
    raw_coords = []

    for elem in root.iter():
        if 'coord' in elem.tag.lower() and elem.text:
            parts = elem.text.strip().replace(',', ' ').split()
            for i in range(0, len(parts) - 2, 3):
                try:
                    lon, lat, alt_m = float(parts[i]), float(parts[i + 1]), float(parts[i + 2])
                    raw_coords.append([lon, lat, alt_m])
                except ValueError:
                    continue

    if not raw_coords:
        raise ValueError("CRITICAL: No valid coordinates found in KML.")

    if len(raw_coords) > max_points:
        # np.linspace guarantees the cap is respected exactly and keeps
        # the departure/arrival points, unlike a fixed stride slice.
        idx = np.linspace(0, len(raw_coords) - 1, max_points, dtype=int)
        raw_coords = [raw_coords[i] for i in idx]
        print(f"  -> Downsampled to {len(raw_coords)} points for optimization.")

    return np.array(raw_coords)


# --------------------------------------------------------------------------
# Step 2: DEM fetch (OpenTopography)
# --------------------------------------------------------------------------
def fetch_dem(lon_min, lat_min, lon_max, lat_max, api_key, timeout=60, retries=2):
    """Download an SRTMGL3 elevation grid and return it as a dict of arrays.
    NODATA cells are masked to 0 so they never poison interpolation later."""
    url = "https://portal.opentopography.org/API/globaldem"
    params = {
        "demtype": "SRTMGL3",
        "south": lat_min, "north": lat_max,
        "west": lon_min, "east": lon_max,
        "outputFormat": "AAIGrid",
        "API_Key": api_key,
    }

    last_error = None
    response = None
    for attempt in range(1, retries + 2):
        try:
            response = requests.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            if not response.text.strip().lower().startswith("ncols"):
                raise ValueError(f"Unexpected DEM response: {response.text[:150]}")
            break
        except (requests.RequestException, ValueError) as e:
            last_error = e
            response = None
            print(f"  -> DEM fetch attempt {attempt} failed ({e}); retrying...")
            time.sleep(2 * attempt)

    if response is None:
        raise ValueError(f"CRITICAL: DEM download failed after retries: {last_error}")

    lines = response.text.splitlines()
    header = {}
    for line in lines[:6]:
        key, val = line.split()[:2]
        header[key.lower()] = val

    ncols = int(header['ncols'])
    nrows = int(header['nrows'])
    xllcorner = float(header['xllcorner'])
    yllcorner = float(header['yllcorner'])
    cellsize = float(header['cellsize'])

    data_block = "\n".join(lines[6:])
    elev_flat = np.loadtxt(StringIO(data_block)).ravel()
    elevation_matrix = np.flipud(elev_flat.reshape((nrows, ncols)))

    # Mask NODATA (e.g. -9999) without clobbering legitimate low/negative
    # elevations (Death Valley, below-sea-level terrain, etc).
    elevation_matrix = np.where(elevation_matrix < NODATA_THRESHOLD, 0.0, elevation_matrix)

    grid_lats = yllcorner + cellsize * np.arange(nrows)
    grid_lons = xllcorner + cellsize * np.arange(ncols)

    return {
        "ncols": ncols, "nrows": nrows,
        "elevation_matrix": elevation_matrix,
        "grid_lats": grid_lats, "grid_lons": grid_lons,
    }


# --------------------------------------------------------------------------
# Step 3: geographic -> print-space scale factors
# --------------------------------------------------------------------------
def compute_scale_factors(lon_min, lon_max, lat_min, lat_max, target_width_mm, z_scale_exaggeration):
    lat_dist = (lat_max - lat_min) * 111000
    lon_dist = (lon_max - lon_min) * 111000 * np.cos(np.radians((lat_min + lat_max) / 2))
    scale_x = target_width_mm / (lon_max - lon_min)
    scale_y = (target_width_mm * (lat_dist / lon_dist)) / (lat_max - lat_min)
    scale_z = (target_width_mm / lon_dist) * z_scale_exaggeration
    width_x_mm = (lon_max - lon_min) * scale_x
    width_y_mm = (lat_max - lat_min) * scale_y
    return scale_x, scale_y, scale_z, width_x_mm, width_y_mm


# --------------------------------------------------------------------------
# Step 4: solid (watertight) terrain block
# --------------------------------------------------------------------------
def build_solid_terrain_mesh(dem, lon_min, lat_min, scale_x, scale_y, scale_z, z_floor):
    """Builds the terrain as a closed solid: a heightmap top surface, a
    flat bottom cap at z_floor, and vertical side walls stitching the two
    together. A bare heightmap (top surface only) is NOT watertight and
    most slicers, including Bambu Studio, will refuse or mis-slice it."""
    grid_lats, grid_lons = dem["grid_lats"], dem["grid_lons"]
    elevation = dem["elevation_matrix"]
    nrows, ncols = dem["nrows"], dem["ncols"]

    jj, ii = np.meshgrid(np.arange(ncols), np.arange(nrows))
    xs = (grid_lons[jj] - lon_min) * scale_x
    ys = (grid_lats[ii] - lat_min) * scale_y
    top_z = np.clip(elevation, 0, None) * scale_z + z_floor

    top_vertices = np.column_stack([xs.ravel(), ys.ravel(), top_z.ravel()])
    bottom_vertices = np.column_stack([xs.ravel(), ys.ravel(), np.full(xs.size, z_floor)])
    vertices = np.vstack([top_vertices, bottom_vertices])
    N = nrows * ncols

    idx = np.arange(N).reshape(nrows, ncols)
    v0, v1 = idx[:-1, :-1].ravel(), idx[:-1, 1:].ravel()
    v2, v3 = idx[1:, :-1].ravel(), idx[1:, 1:].ravel()

    top_faces = np.empty((2 * len(v0), 3), dtype=np.int64)
    top_faces[0::2] = np.column_stack([v0, v1, v2])
    top_faces[1::2] = np.column_stack([v1, v3, v2])

    b0, b1, b2, b3 = v0 + N, v1 + N, v2 + N, v3 + N
    bottom_faces = np.empty((2 * len(v0), 3), dtype=np.int64)
    bottom_faces[0::2] = np.column_stack([b0, b2, b1])
    bottom_faces[1::2] = np.column_stack([b2, b3, b1])

    top_row, right_col = idx[0, :], idx[:, -1]
    bottom_row, left_col = idx[-1, :][::-1], idx[:, 0][::-1]
    perimeter = np.concatenate([top_row, right_col[1:], bottom_row[1:], left_col[1:-1]])

    side_faces = []
    n_perim = len(perimeter)
    for k in range(n_perim):
        a, b = perimeter[k], perimeter[(k + 1) % n_perim]
        side_faces.append([a, b + N, b])
        side_faces.append([a, a + N, b + N])
    side_faces = np.array(side_faces, dtype=np.int64)

    faces = np.vstack([top_faces, bottom_faces, side_faces])
    return trimesh.Trimesh(vertices=vertices, faces=faces)


# --------------------------------------------------------------------------
# Step 5: building extrusion from OpenStreetMap
# --------------------------------------------------------------------------
def build_city_mesh(lon_min, lat_min, lon_max, lat_max, interp, scale_x, scale_y, scale_z, z_floor):
    try:
        tags = {"building": True, "aeroway": ["hangar", "terminal"], "man_made": True}
        try:
            gdf = ox.features_from_bbox(bbox=(lon_min, lat_min, lon_max, lat_max), tags=tags)
        except TypeError:
            gdf = ox.features_from_bbox(north=lat_max, south=lat_min, east=lon_max, west=lon_min, tags=tags)
    except Exception as e:
        print(f"  -> OSM fetch failed: {e}")
        return trimesh.Trimesh()

    if gdf.empty:
        return trimesh.Trimesh()

    building_meshes = []
    skipped = 0
    for _, row in gdf.iterrows():
        geom = row.geometry
        polygons = [geom] if isinstance(geom, Polygon) else getattr(geom, 'geoms', [])

        for poly in polygons:
            if not isinstance(poly, Polygon) or not poly.is_valid or poly.exterior is None:
                skipped += 1
                continue
            try:
                x_coords, y_coords = poly.exterior.coords.xy
                pts = np.column_stack(((np.array(x_coords) - lon_min) * scale_x,
                                        (np.array(y_coords) - lat_min) * scale_y))
                b_ground_z = interp([[poly.centroid.y, poly.centroid.x]])[0] * scale_z

                height_m = 8.0
                if 'height' in row and pd.notnull(row['height']):
                    cleaned = ''.join(c for c in str(row['height']) if c.isdigit() or c == '.')
                    if cleaned.count('.') <= 1 and cleaned not in ('', '.'):
                        height_m = float(cleaned)
                elif 'building:levels' in row and pd.notnull(row['building:levels']):
                    height_m = float(row['building:levels']) * 3.5

                b_mesh = trimesh.creation.extrude_polygon(Polygon(pts), height=height_m * scale_z)
                b_mesh.apply_translation([0, 0, b_ground_z + z_floor])
                building_meshes.append(b_mesh)
            except Exception:
                skipped += 1
                continue

    if skipped:
        print(f"  -> Skipped {skipped} malformed/unbuildable building footprints.")

    if not building_meshes:
        return trimesh.Trimesh()

    chunk_size = 250
    city_chunks = [trimesh.util.concatenate(building_meshes[i:i + chunk_size])
                   for i in range(0, len(building_meshes), chunk_size)]
    return trimesh.util.concatenate(city_chunks)


# --------------------------------------------------------------------------
# Step 6: pylons + flight path tube
# --------------------------------------------------------------------------
def build_pylons_and_tube(path_points, ground_z, z_floor, pylon_radius=0.8, tube_radius=0.6, max_pylons=40):
    pylon_meshes = []
    step = max(1, len(path_points) // max_pylons)

    for i in range(0, len(path_points), step):
        p_top = path_points[i]
        height = p_top[2] - ground_z[i]
        if height > 0:
            pylon = trimesh.creation.cylinder(radius=pylon_radius, height=height)
            pylon.apply_translation([p_top[0], p_top[1], ground_z[i] + z_floor + (height / 2)])
            pylon_meshes.append(pylon)

    all_pylons = trimesh.util.concatenate(pylon_meshes) if pylon_meshes else trimesh.Trimesh()

    shifted_points = path_points.copy()
    shifted_points[:, 2] += z_floor
    circle_profile = Point(0, 0).buffer(tube_radius)

    if len(shifted_points) < 2:
        raise ValueError("CRITICAL: Need at least 2 flight-path points to build a track tube.")

    try:
        path_tube = trimesh.creation.sweep_polygon(polygon=circle_profile, path=shifted_points)
    except Exception as e:
        print(f"  -> sweep_polygon failed ({e}); falling back to per-segment cylinders.")
        segments = [trimesh.creation.cylinder(radius=tube_radius, segment=[shifted_points[i], shifted_points[i + 1]])
                    for i in range(len(shifted_points) - 1)]
        path_tube = trimesh.util.concatenate(segments)

    return all_pylons, path_tube


# --------------------------------------------------------------------------
# Step 7: base plate
# --------------------------------------------------------------------------
def build_base_plate(width_x_mm, width_y_mm, base_thickness_mm, base_margin_mm):
    """A flat rectangular plinth sitting on the print bed (z=0 to
    base_thickness_mm) that the terrain/buildings/pylons sit on top of.
    Gives the model a stable footprint and first-layer adhesion area."""
    plate_w = width_x_mm + 2 * base_margin_mm
    plate_d = width_y_mm + 2 * base_margin_mm
    plate = trimesh.creation.box([plate_w, plate_d, base_thickness_mm])
    # trimesh boxes are centered at the origin; shift so it spans
    # x:[-margin, width_x+margin], y:[-margin, width_y+margin], z:[0, thickness]
    plate.apply_translation([width_x_mm / 2, width_y_mm / 2, base_thickness_mm / 2])
    return plate


def check_bambu_bed_fit(width_x_mm, width_y_mm, height_z_mm, bed_size_mm):
    if width_x_mm > bed_size_mm or width_y_mm > bed_size_mm or height_z_mm > bed_size_mm:
        print(f"  -> WARNING: model footprint is {width_x_mm:.0f} x {width_y_mm:.0f} x "
              f"{height_z_mm:.0f} mm, which exceeds the {bed_size_mm:.0f} mm Bambu Lab "
              f"build volume you specified. Reduce --target-width or increase --bed-size.")
    else:
        print(f"  -> Model footprint {width_x_mm:.0f} x {width_y_mm:.0f} x {height_z_mm:.0f} mm "
              f"fits within the {bed_size_mm:.0f} mm Bambu Lab build volume.")


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------
def automate_flight_to_3d(kml_path, api_key, target_width_mm=150, padding_deg=0.02,
                           z_scale_exaggeration=2.5, max_points=800,
                           base_thickness_mm=4.0, base_margin_mm=6.0,
                           bed_size_mm=256.0, add_base=True):
    if not api_key:
        raise ValueError("CRITICAL: You must provide an OpenTopography API key.")
    if not os.path.isfile(kml_path):
        raise ValueError(f"CRITICAL: KML file not found: {kml_path}")

    print(f"Step 1: Parsing KML ({kml_path}) and auto-detecting bounding box...")
    flight_path = parse_kml_coordinates(kml_path, max_points)
    lons, lats, alts_m = flight_path[:, 0], flight_path[:, 1], flight_path[:, 2]

    lon_min, lon_max = lons.min() - padding_deg, lons.max() + padding_deg
    lat_min, lat_max = lats.min() - padding_deg, lats.max() + padding_deg

    print("\nStep 2: Fetching Real Earth DEM from OpenTopography (SRTMGL3)...")
    dem = fetch_dem(lon_min, lat_min, lon_max, lat_max, api_key)

    print("\nStep 3: Calculating geographic scale ratios...")
    scale_x, scale_y, scale_z, width_x_mm, width_y_mm = compute_scale_factors(
        lon_min, lon_max, lat_min, lat_max, target_width_mm, z_scale_exaggeration)

    z_floor = base_thickness_mm if add_base else 0.0

    print("\nStep 4: Generating solid 3D Mountain Terrain block...")
    terrain_mesh = build_solid_terrain_mesh(dem, lon_min, lat_min, scale_x, scale_y, scale_z, z_floor)

    print("\nStep 5: Interpolating ground altitudes...")
    interp = RegularGridInterpolator((dem["grid_lats"], dem["grid_lons"]), dem["elevation_matrix"],
                                      bounds_error=False, fill_value=0)
    ground_z = interp(np.column_stack((lats, lons))) * scale_z
    path_points = np.column_stack(((lons - lon_min) * scale_x, (lats - lat_min) * scale_y, alts_m * scale_z))

    print("\nStep 6: Fetching OpenStreetMap Building Data...")
    city_mesh = build_city_mesh(lon_min, lat_min, lon_max, lat_max, interp, scale_x, scale_y, scale_z, z_floor)

    print("\nStep 7: Engineering Pylons and Track Tube...")
    all_pylons, path_tube = build_pylons_and_tube(path_points, ground_z, z_floor)

    print("\nStep 8: Assembling 3MF File with Baked Colors...")
    terrain_mesh.metadata['name'] = '1_Terrain'
    terrain_mesh.visual.face_colors = [130, 130, 130, 255]

    scene_parts = [terrain_mesh]

    if not city_mesh.is_empty:
        city_mesh.metadata['name'] = '2_Buildings'
        city_mesh.visual.face_colors = [240, 240, 240, 255]
        scene_parts.append(city_mesh)

    if not all_pylons.is_empty:
        all_pylons.metadata['name'] = '3_Pylons'
        all_pylons.visual.face_colors = [30, 30, 30, 255]
        scene_parts.append(all_pylons)

    path_tube.metadata['name'] = '4_FlightPath'
    path_tube.visual.face_colors = [255, 100, 0, 255]
    scene_parts.append(path_tube)

    height_z_mm = float(terrain_mesh.vertices[:, 2].max())

    if add_base:
        print("\nStep 9: Adding print bed base plate...")
        base_plate = build_base_plate(width_x_mm, width_y_mm, base_thickness_mm, base_margin_mm)
        base_plate.metadata['name'] = '0_BasePlate'
        base_plate.visual.face_colors = [60, 60, 60, 255]
        scene_parts.insert(0, base_plate)

    check_bambu_bed_fit(width_x_mm + 2 * base_margin_mm, width_y_mm + 2 * base_margin_mm,
                         height_z_mm, bed_size_mm)

    master_scene = trimesh.Scene(scene_parts)
    base_name = os.path.splitext(os.path.basename(kml_path))[0]
    master_filename = f"{base_name}_Assembly.3mf"
    master_scene.export(master_filename)
    print(f"\n[SUCCESS] Saved to {master_filename}")
    return master_filename


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert ForeFlight KML to a Bambu-Lab-ready 3D Multi-Color Assembly")
    parser.add_argument("kml_file", help="Path to the .kml track log")
    parser.add_argument("--api-key", help="OpenTopography API Key (or set OPENTOPO_API_KEY env var)")
    parser.add_argument("--target-width", type=float, default=150, help="Model width in mm (default: 150)")
    parser.add_argument("--base-thickness", type=float, default=4.0, help="Base plate thickness in mm (default: 4)")
    parser.add_argument("--base-margin", type=float, default=6.0, help="Base plate margin around terrain in mm (default: 6)")
    parser.add_argument("--bed-size", type=float, default=256.0,
                         help="Printer bed size in mm for the fit check, e.g. 180 for A1 mini, 256 for X1C/P1S/A1 (default: 256)")
    parser.add_argument("--no-base", action="store_true", help="Skip adding the base plate")

    args = parser.parse_args()
    api_key = args.api_key or os.getenv("OPENTOPO_API_KEY")

    try:
        automate_flight_to_3d(
            args.kml_file, api_key,
            target_width_mm=args.target_width,
            base_thickness_mm=args.base_thickness,
            base_margin_mm=args.base_margin,
            bed_size_mm=args.bed_size,
            add_base=not args.no_base,
        )
    except ValueError as e:
        print(str(e))
        sys.exit(1)
