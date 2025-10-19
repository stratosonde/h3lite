#!/usr/bin/env python3
"""
Visual Grid Test for H3Lite Region Coverage

This script creates a world map visualization showing where h3lite
successfully finds regions vs where it doesn't (using the actual C library).
Uses the same plotting method as regions_h3_conversion.py.
"""

import subprocess
import sys
import os

# Install required packages if missing
def install_if_missing(package):
    try:
        __import__(package)
    except ImportError:
        print(f"Installing {package}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])

required_packages = ['h3', 'geopandas', 'contextily', 'matplotlib', 'shapely']
for package in required_packages:
    install_if_missing(package)

import h3
import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as cx
from shapely.geometry import Polygon

# Region colors (matching regions_h3_conversion.py style)
REGION_COLORS = {
    'Not Found': '#CCCCCC',  # Light Gray - no region found
    'Unknown': '#808080',     # Dark Gray - Unknown LoRaWAN region
    'US915': '#FF5733',       # Red
    'EU868': '#33FF57',       # Green
    'AU915': '#3357FF',       # Blue
    'AS923-1': '#FF33E9',     # Pink
    'AS923-2': '#E933FF',     # Purple
    'AS923-3': '#33E9FF',     # Cyan
    'AS923-4': '#E9FF33',     # Yellow
    'KR920': '#33FFF6',       # Light cyan
    'IN865': '#F6FF33',       # Light yellow
    'RU864': '#7F33FF',       # Purple
    'CN470': '#FF8333',       # Orange
    'EU433': '#83FF33',       # Lime
    'CD900-1A': '#FF6B9D',    # Pink
}

def check_crosses_idl(polygon):
    """Check if a polygon crosses the International Date Line (180°/-180° longitude).
    From regions_h3_conversion.py"""
    if polygon is None:
        return False
    
    if hasattr(polygon, 'exterior'):
        coords = list(polygon.exterior.coords)
    else:
        return False
    
    for i in range(len(coords) - 1):
        lon1, _ = coords[i]
        lon2, _ = coords[i + 1]
        
        if ((lon1 > 150 and lon2 < -150) or 
            (lon1 < -150 and lon2 > 150)):
            return True
    
    return False

def split_polygon_at_dateline(polygon):
    """Split a polygon that crosses the international date line into two separate polygons.
    From regions_h3_conversion.py"""
    from shapely.geometry import MultiPolygon
    
    if not check_crosses_idl(polygon):
        return polygon
    
    coords = list(polygon.exterior.coords)
    
    west_coords = []
    east_coords = []
    
    for x, y in coords:
        if x < 0:  # Western hemisphere
            west_coords.append((x, y))
        else:      # Eastern hemisphere
            east_coords.append((x, y))
    
    polygons = []
    
    if len(west_coords) >= 3:
        if west_coords[0] != west_coords[-1]:
            west_coords.append(west_coords[0])
        try:
            west_poly = Polygon(west_coords)
            if west_poly.is_valid:
                polygons.append(west_poly)
        except:
            pass
    
    if len(east_coords) >= 3:
        if east_coords[0] != east_coords[-1]:
            east_coords.append(east_coords[0])
        try:
            east_poly = Polygon(east_coords)
            if east_poly.is_valid:
                polygons.append(east_poly)
        except:
            pass
    
    if not polygons:
        return polygon
    
    if len(polygons) > 1:
        return MultiPolygon(polygons)
    else:
        return polygons[0]

def plot_df(df, column=None, ax=None, title=None):
    """Plot based on the `geometry` column of a GeoPandas dataframe
    Same as regions_h3_conversion.py"""
    df = df.copy()
    df = df.to_crs(epsg=3857)  # web mercator

    if ax is None:
        _, ax = plt.subplots(figsize=(12, 10))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    df.plot(
        ax=ax,
        alpha=0.5, edgecolor='k',
        column=column, categorical=True,
        legend=True, legend_kwds={'loc': 'upper left'}
    )
    cx.add_basemap(ax, crs=df.crs, source=cx.providers.CartoDB.Positron)
    
    if title:
        ax.set_title(title)

def get_region_from_h3lite(lat, lng):
    """Call h3lite to get the region for a lat/lng point"""
    code = f"""
#include "include/h3lite.h"
#include <stdio.h>
int main() {{
    h3liteInit();
    NearestRegionsInfo info = findNearestRegions({lat}, {lng}, 3);
    if (info.numRegions > 0) {{
        printf("%s\\n", info.regions[0].regionName);
    }} else {{
        printf("Unknown\\n");
    }}
    return 0;
}}
"""
    
    with open('/tmp/test_point.c', 'w') as f:
        f.write(code)
    
    try:
        compile_result = subprocess.run([
            'gcc', '-O3', '-I.', '/tmp/test_point.c',
            './obj/h3lite.o', './obj/h3lite_faceijk.o', './obj/h3lite_basecells.o',
            './obj/h3lite_neighbor.o', './obj/h3lite_regions_table.o',
            '-o', '/tmp/test_point', '-lm'
        ], capture_output=True, text=True, cwd='.')
        
        if compile_result.returncode != 0:
            return 'Unknown'
        
        result = subprocess.run(['/tmp/test_point'], capture_output=True, text=True)
        if result.returncode != 0:
            return 'Unknown'
            
        region = result.stdout.strip().split('\n')[-1] if result.stdout.strip() else 'Unknown'
        return region if region else 'Unknown'
    except Exception:
        return 'Unknown'

def create_visual_grid(resolution=2):
    """Create a visual grid showing region coverage using H3 cells
    
    Gets all H3 cells at the specified resolution and tests each one exactly once.
    This ensures no duplicate cells in the visualization.
    
    Args:
        resolution: H3 resolution to use (default 2)
    """
    print(f"Creating visual grid at H3 resolution {resolution}...")
    print(f"Getting all H3 cells at resolution {resolution}...")
    
    # Get all H3 cells at this resolution
    all_cells = h3.get_res0_cells()
    cells_to_test = set()
    
    # Get all cells at target resolution
    for base_cell in all_cells:
        try:
            children = h3.cell_to_children(base_cell, resolution)
            cells_to_test.update(children)
        except:
            continue
    
    print(f"Found {len(cells_to_test)} total H3 cells at resolution {resolution}")
    
    # Collect H3 cells by region
    region_cells = {}
    
    tested = 0
    total_cells = len(cells_to_test)
    
    for h3_cell in cells_to_test:
        # Get center point of this H3 cell
        try:
            lat, lng = h3.cell_to_latlng(h3_cell)
        except:
            continue
        
        # Get region from h3lite for this cell's center
        region = get_region_from_h3lite(lat, lng)
        if region == 'Unknown':
            region = 'Not Found'
        
        # Add to region's cell list
        if region not in region_cells:
            region_cells[region] = set()
        region_cells[region].add(h3_cell)
        
        tested += 1
        if tested % 1000 == 0:
            print(f"Tested {tested}/{total_cells} cells ({100*tested/total_cells:.1f}%)")
    
    # Create geometries for each region
    all_geometries = []
    all_regions = []
    
    for region, cells in region_cells.items():
        cells_list = list(cells)
        print(f"Region {region}: {len(cells_list)} H3 cells")
        
        # Skip "Not Found" region (likely ocean) but keep "Unknown" (found but unknown LoRaWAN region)
        if region == 'Not Found':
            print(f"  Skipping {len(cells_list)} 'Not Found' cells from visualization (likely ocean)")
            continue
        
        # Convert H3 cells to polygons
        for cell in cells_list:
            try:
                # Get cell boundary and create polygon
                # cell_to_boundary returns [(lat, lng), ...] - need to swap to (lng, lat) for Polygon
                boundary = h3.cell_to_boundary(cell)
                # Swap lat/lng to lng/lat for proper polygon creation
                boundary_swapped = [(lng, lat) for lat, lng in boundary]
                poly = Polygon(boundary_swapped)
                
                # Check if polygon crosses dateline and split if needed
                if check_crosses_idl(poly):
                    split_geom = split_polygon_at_dateline(poly)
                    # If split into MultiPolygon, add each part separately
                    if hasattr(split_geom, 'geoms'):
                        for part in split_geom.geoms:
                            all_geometries.append(part)
                            all_regions.append(region)
                    else:
                        all_geometries.append(split_geom)
                        all_regions.append(region)
                else:
                    all_geometries.append(poly)
                    all_regions.append(region)
            except Exception as e:
                print(f"Error converting cell {cell}: {e}")
    
    # Create GeoDataFrame
    if not all_geometries:
        print("No geometries to plot!")
        return
    
    gdf = gpd.GeoDataFrame({
        'geometry': all_geometries,
        'region': all_regions
    }, crs='EPSG:4326')
    
    # Create the plot using the same method as regions_h3_conversion.py
    plt.figure(figsize=(16, 12))
    plot_df(gdf, column='region', title=f'H3Lite Region Coverage Map (H3 Resolution: {resolution})')
    plt.tight_layout()
    plt.savefig('h3lite_coverage_visual.png', dpi=300, bbox_inches='tight')
    print(f"\nSaved visualization to h3lite_coverage_visual.png")
    
    # Print statistics
    print("\n=== Region Coverage Statistics ===")
    total = len(all_regions)
    
    from collections import Counter
    region_counts = Counter(all_regions)
    for region_name, count in sorted(region_counts.items()):
        print(f"{region_name}: {count} H3 cells ({100*count/total:.1f}%)")

if __name__ == "__main__":
    print("H3Lite Visual Grid Coverage Test")
    print("=================================\n")
    
    # Use resolution 3 for a good balance of detail and performance
    # Resolution 2: ~5,882 cells
    # Resolution 3: ~41,162 cells
    # Resolution 4: ~288,122 cells
    create_visual_grid(resolution=3)
    
    print("\nDone!")
