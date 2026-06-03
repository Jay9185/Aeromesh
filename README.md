# AeroMesh: Flight Telemetry to 3D Assembly 🛩️ 🏔️

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

A Python geospatial pipeline that transforms 2D aviation track logs (ForeFlight KMLs) into physical, multi-color 3D printable models. 

<img width="1920" height="1080" alt="model" src="https://github.com/user-attachments/assets/1b895ae3-fd6a-4344-a69e-e26e4cf27ed1" />

## The Concept
Looking at a 2D line on an iPad screen doesn't capture altitude changes or the physical scale of the terrain you are flying over. Generating 3D prints from GPS tracks usually involves heavy manual CAD work, slicing solid blocks of plastic, and tedious hand-painting. 

**AeroMesh automates this.** It reads a flight log, downloads the real-world Earth below it, extrudes the city skyline, engineers its own intermittent support pylons based on true AGL (Altitude Above Ground Level), and bakes AMS colors directly into a `.3mf` file.

You just drag it into Bambu Studio, assign your spools, and hit print.

## How It Works Under the Hood
1. **Coordinate Normalization:** Parses WGS 84 (Lat/Lon/Alt) telemetry and dynamically frames a bounding box, mapping the spherical Earth to a flat Cartesian plane in millimeters.
2. **Topography (DEM):** Hits the OpenTopography API (SRTMGL3) to download a 90m-resolution Digital Elevation Model of your exact flight area.
3. **The Skyline:** Uses `osmnx` to query OpenStreetMap building footprints and extrudes them to their real-world structural heights.
4. **Intermittent Pylons:** Instead of relying on slicer-generated supports, it calculates exact AGL using Scipy's `RegularGridInterpolator` and drops support cylinders directly to the mountain ridges.
5. **Memory Management:** Implements aggressive downsampling and chunked mesh concatenation to ensure it doesn't OOM crash on dense urban flights.
6. **3MF Export:** Bakes RGBA colors (Grey Terrain, White Buildings, Black Pylons, Orange Track) directly into a master `.3mf` container.

## Installation

Clone the repository and install the geospatial requirements:

```bash
git clone [https://github.com/YOUR_USERNAME/AeroMesh.git](https://github.com/YOUR_USERNAME/AeroMesh.git)
cd AeroMesh
pip install -r requirements.txt
