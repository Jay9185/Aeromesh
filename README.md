# AeroMesh: Flight Telemetry to 3D Assembly 🛩️ 🏔️

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

A Python geospatial pipeline that transforms 2D aviation track logs (ForeFlight KMLs) into physical, multi-color 3D printable models. 

<img width="10240" height="5400" alt="new render1" src="https://github.com/user-attachments/assets/bfe4687f-f310-4ae0-bade-fc688c952c76" />



## The Concept
Looking at a 2D line on an iPad screen doesn't capture altitude changes or the physical scale of the terrain you are flying over. Generating 3D prints from GPS tracks usually involves heavy manual CAD work, slicing solid blocks of plastic, and tedious hand-painting. 

**AeroMesh automates this.** It reads a flight log, downloads the real-world Earth below it, extrudes the city skyline, engineers its own intermittent support pylons based on true AGL (Altitude Above Ground Level), and bakes AMS colors directly into a `.3mf` file.

You just drag it into Slicer, assign your spools, and hit print.

## How It Works Under the Hood
1. **Coordinate Normalization:** Parses WGS 84 (Lat/Lon/Alt) telemetry and dynamically frames a bounding box, mapping the spherical Earth to a flat Cartesian plane in millimeters.
2. **Topography (DEM):** Hits the OpenTopography API (SRTMGL3) to download a 90m-resolution Digital Elevation Model of your exact flight area.
3. **The Skyline:** Uses `osmnx` to query OpenStreetMap building footprints and extrudes them to their real-world structural heights.
4. **Intermittent Pylons:** Instead of relying on slicer-generated supports, it calculates exact AGL using Scipy's `RegularGridInterpolator` and drops support cylinders directly to the mountain ridges.
5. **Memory Management:** Implements aggressive downsampling and chunked mesh concatenation to ensure it doesn't OOM crash on dense urban flights.
6. **3MF Export:** Bakes RGBA colors (Grey Terrain, White Buildings, Black Pylons, Orange Track) directly into a master `.3mf` container.

## ⚙️ How It Works Under the Hood
1. **Coordinate Normalization:** Parses WGS 84 (Lat/Lon/Alt) telemetry and dynamically frames a bounding box, mapping the spherical Earth to a flat Cartesian plane in millimeters.
2. **Topography (DEM):** Hits the OpenTopography API (SRTMGL3) to download a 90m-resolution Digital Elevation Model of your exact flight area.
3. **The Skyline:** Uses `osmnx` to query OpenStreetMap building footprints and extrudes them to their real-world structural heights.
4. **Intermittent Pylons:** Instead of relying on slicer-generated supports, it calculates exact AGL using Scipy's `RegularGridInterpolator` and drops support cylinders directly to the mountain ridges.
5. **3MF Export:** Bakes RGBA colors (Grey Terrain, White Buildings, Black Pylons, Orange Track) directly into a master `.3mf` container.

---

## 🚀 Quickstart Guide

### 1. Prerequisites
You need a free OpenTopography API key to fetch the elevation data.
* Go to [OpenTopography](https://portal.opentopography.org/login) and create a free account.
* Navigate to "My Account" -> "MyOpenTopo Authorizations and Quotas" and request an API Key.

### 2. Installation
Clone the repository and install the required geospatial libraries. *(Note: Using a virtual environment is highly recommended).*

git clone git@github.com:Jay9185/Aeromesh.git
cd Aeromesh
pip install -r requirements.txt
### 3. Usage
Run the script via the command line, pointing it to your `.kml` flight log and passing your API key.

```bash
python Aeromesh.py path/to/your/flight.kml --api-key "YOUR_API_KEY_HERE"
```

**Optional Arguments for Customization:**
You can customize the size and thickness of the final model to fit your specific 3D printer bed.
* `--target-width`: Model width in mm (default: 150)
* `--base-thickness`: Base plate thickness in mm (default: 4.0)
* `--bed-size`: Printer bed size in mm to check if it fits (default: 256.0 for Bambu X1/P1)
* `--no-base`: Skips generating the flat base plate

**Example with custom size:**
```bash
python Aeromesh.py flight.kml --api-key "YOUR_KEY" --target-width 200 --bed-size 256
```

---

## 🖨️ Slicer Instructions (For Bambu Studio / PrusaSlicer)

Because AeroMesh generates a `.3mf` assembly with embedded objects, setting up the multi-color print is frictionless.

1. **Import:** Drag and drop the generated `_Assembly.3mf` into your slicer. 
2. **Assign Colors:** Look at the "Objects" or "Process" tab in your slicer. You will see distinct parts (`1_Terrain`, `2_Buildings`, `3_Pylons`, `4_FlightPath`). Assign your AMS/MMU filaments to each part. 
3. **Supports:** Depending on your printer's bridging capabilities, you *may* need to enable standard slicer supports for the sections of the flight path that bridge between the engineered pylons. 

## 🤝 Contributing
I coded the mesh generation and AGL math entirely "blind" without a 3D printer currently in my studio. If you have an AMS/IDEX setup and run your flight logs through this, please open an issue or submit a PR with photos of the physical results or any physics/overhang bugs you encounter!
