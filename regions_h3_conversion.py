#!/usr/bin/env python3
# regions_h3_conversion.py - Convert regions.geojson to H3 Hexagons
#
# This script demonstrates how to convert the regions.geojson file to Uber's H3 
# hexagonal grid system and visualize each region at different resolutions (2-6).

import sys
import subprocess
import os
import gc
import json

# Function to install dependencies if they're missing
def install_if_missing(package):
    try:
        __import__(package)
        print(f"{package} is already installed")
    except ImportError:
        print(f"Installing {package}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        print(f"{package} installed successfully")

# Install required packages
required_packages = ['h3', 'geopandas', 'contextily', 'matplotlib', 'shapely', 'tqdm', 'fiona', 'pandas']
for package in required_packages:
    install_if_missing(package)

print("All dependencies checked/installed.")

# Import necessary libraries
import h3
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import contextily as cx
from shapely.geometry import shape, Polygon, MultiPolygon
from tqdm import tqdm
import numpy as np

# Helper functions for plotting
def plot_df(df, column=None, ax=None, title=None):
    """Plot based on the `geometry` column of a GeoPandas dataframe"""
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

def plot_shape(shape, ax=None, title=None):
    """Plot a single shape"""
    df = gpd.GeoDataFrame({'geometry': [shape]}, crs='EPSG:4326')
    plot_df(df, ax=ax, title=title)

def plot_cells(cells, ax=None, title=None):
    """Plot H3 cells with special handling for cells across the international date line"""
    if not cells:  # Check if cells list is empty
        print("Warning: No cells to plot")
        if ax is None:
            _, ax = plt.subplots(figsize=(12, 10))
        ax.text(0.5, 0.5, "No H3 cells generated", 
                horizontalalignment='center',
                verticalalignment='center',
                transform=ax.transAxes,
                fontsize=14)
        if title:
            ax.set_title(title)
        return
    
    try:
        # Group cells by whether they're likely on the east or west side of the date line
        # This is a heuristic approach based on examining the cell center coordinates
        east_cells = []
        west_cells = []
        
        for cell in cells:
            # Get the center coordinates of the cell
            lat, lng = h3.cell_to_latlng(cell)
            # Classify based on longitude - use a wider band to catch all cells near the date line
            if lng > 150 or lng < -150:  # Near the international date line
                if lng > 0:
                    east_cells.append(cell)
                else:
                    west_cells.append(cell)
            else:
                # For cells away from the date line, just add to both groups
                # They'll be handled correctly when converted to shapes
                east_cells.append(cell)
                west_cells.append(cell)
        
        # Create a plot for each side if needed
        if ax is None:
            _, ax = plt.subplots(figsize=(12, 10))
        
        if east_cells and west_cells and (set(east_cells) != set(west_cells)):
            print(f"Splitting {len(cells)} cells into East ({len(east_cells)}) and West ({len(west_cells)}) groups")
            
            # Create shapes for each side
            try:
                if east_cells:
                    east_shape = h3.cells_to_h3shape(east_cells)
                    east_gdf = gpd.GeoDataFrame({'geometry': [east_shape]}, crs='EPSG:4326')
                    east_gdf = east_gdf.to_crs(epsg=3857)  # web mercator
                    east_gdf.plot(ax=ax, alpha=0.5, edgecolor='k', color='blue')
                
                if west_cells:
                    west_shape = h3.cells_to_h3shape(west_cells)
                    west_gdf = gpd.GeoDataFrame({'geometry': [west_shape]}, crs='EPSG:4326')
                    west_gdf = west_gdf.to_crs(epsg=3857)  # web mercator
                    west_gdf.plot(ax=ax, alpha=0.5, edgecolor='k', color='blue')
                
                # Add basemap
                cx.add_basemap(ax, crs='EPSG:3857', source=cx.providers.CartoDB.Positron)
            except Exception as e:
                print(f"Error plotting split cells: {e}")
                # Fall back to simple plotting
                for cell in cells:
                    bounds = h3.cell_to_boundary(cell, geo_json=True)
                    polygon = Polygon(bounds)
                    gdf = gpd.GeoDataFrame({'geometry': [polygon]}, crs='EPSG:4326')
                    gdf = gdf.to_crs(epsg=3857)  # web mercator
                    gdf.plot(ax=ax, alpha=0.2, edgecolor='k', color='green')
        else:
            # If cells aren't split across the date line, use standard plotting
            try:
                shape = h3.cells_to_h3shape(cells)
                plot_shape(shape, ax=ax, title=None)  # Title will be set later
            except Exception as e:
                print(f"Error in standard cell plotting: {e}")
                # Fall back to individual cell plotting
                for cell in cells:
                    try:
                        bounds = h3.cell_to_boundary(cell, geo_json=True)
                        polygon = Polygon(bounds)
                        gdf = gpd.GeoDataFrame({'geometry': [polygon]}, crs='EPSG:4326')
                        gdf = gdf.to_crs(epsg=3857)  # web mercator
                        gdf.plot(ax=ax, alpha=0.2, edgecolor='k', color='green')
                    except Exception as cell_e:
                        print(f"Error plotting individual cell: {cell_e}")
        
        # Set title if provided
        if title:
            ax.set_title(title)
    except Exception as e:
        print(f"Error plotting cells: {e}")
        if ax is None:
            _, ax = plt.subplots(figsize=(12, 10))
        ax.text(0.5, 0.5, f"Error plotting cells: {e}", 
                horizontalalignment='center',
                verticalalignment='center',
                transform=ax.transAxes,
                fontsize=10,
                wrap=True)
        if title:
            ax.set_title(title)

def check_crosses_idl(polygon):
    """Check if a polygon crosses the International Date Line (180°/-180° longitude).
    
    This function checks for consecutive coordinates that jump from one side of
    the date line to the other, which is a clear indicator of date line crossing.
    """
    if polygon is None:
        return False
    
    # Extract coordinates from the polygon exterior
    if hasattr(polygon, 'exterior'):
        coords = list(polygon.exterior.coords)
    else:
        return False
    
    # Specifically check for segments that cross from +175 to -175 longitude (or vice versa)
    # which clearly indicates crossing the international date line
    for i in range(len(coords) - 1):
        lon1, _ = coords[i]
        lon2, _ = coords[i + 1]
        
        # Check if this segment crosses the date line
        if ((lon1 > 150 and lon2 < -150) or 
            (lon1 < -150 and lon2 > 150)):
            return True
    
    return False

def split_polygon_at_dateline(polygon):
    """Split a polygon that crosses the international date line into two separate polygons.
    
    Returns a MultiPolygon containing the split parts.
    This is a simple approach that handles the date line crossing.
    """
    if not check_crosses_idl(polygon):
        return polygon
    
    # Get coordinates
    coords = list(polygon.exterior.coords)
    
    # Simply create two separate polygons:
    # 1. One with all western hemisphere coordinates (-180 to 0)
    # 2. One with all eastern hemisphere coordinates (0 to 180)
    west_coords = []
    east_coords = []
    
    # First pass: classify all points directly
    for x, y in coords:
        if x < 0:  # Western hemisphere
            west_coords.append((x, y))
        else:      # Eastern hemisphere
            east_coords.append((x, y))
    
    # Ensure we have closed polygons
    if west_coords and west_coords[0] != west_coords[-1]:
        west_coords.append(west_coords[0])
        
    if east_coords and east_coords[0] != east_coords[-1]:
        east_coords.append(east_coords[0])
    
    # Create polygons for each side if they have enough points
    polygons = []
    
    # Close the polygons by adding the first point again at the end
    if len(west_coords) >= 3:
        if west_coords[0] != west_coords[-1]:
            west_coords.append(west_coords[0])
        try:
            west_poly = Polygon(west_coords)
            if west_poly.is_valid:
                polygons.append(west_poly)
        except Exception as e:
            print(f"Error creating west polygon: {e}")
    
    if len(east_coords) >= 3:
        if east_coords[0] != east_coords[-1]:
            east_coords.append(east_coords[0])
        try:
            east_poly = Polygon(east_coords)
            if east_poly.is_valid:
                polygons.append(east_poly)
        except Exception as e:
            print(f"Error creating east polygon: {e}")
    
    # Return original polygon if splitting failed
    if not polygons:
        return polygon
    
    # Return MultiPolygon if we have multiple parts, otherwise return the single Polygon
    if len(polygons) > 1:
        return MultiPolygon(polygons)
    else:
        return polygons[0]

def convert_geometry_to_h3(geom, resolution=7, verbose=False):
    """Enhanced convert_geometry_to_h3 function with improved error handling and dateline crossing fix
    
    Important: H3 requires coordinates in (lat, lng) order while GeoJSON uses (lng, lat)
    """
    cells = []
    
    if geom is None:
        print("Warning: Received None geometry")
        return cells
    
    if verbose:
        print(f"Converting {geom.geom_type} to H3 cells (resolution {resolution})")
    
    # Check if we need to handle date line crossing
    date_line_fixed = False
        
    try:
        if geom.geom_type == 'MultiPolygon':
            if verbose:
                print(f"Processing MultiPolygon with {len(geom.geoms)} parts")
            
            # Handle MultiPolygon by processing each polygon separately
            for poly_idx, poly in enumerate(geom.geoms):
                if verbose:
                    print(f"  Processing polygon #{poly_idx+1} with {len(poly.exterior.coords)} points")
                
                # Skip polygons with too few points
                if len(poly.exterior.coords) < 4:  # < 4 because first and last points are the same
                    print(f"  Warning: Skipping polygon with only {len(poly.exterior.coords)} points")
                    continue
                
                # Check if this polygon crosses the international date line
                if check_crosses_idl(poly):
                    if verbose:
                        print(f"  Polygon #{poly_idx+1} crosses the international date line - splitting")
                    
                    # Split the polygon at the date line
                    try:
                        split_geom = split_polygon_at_dateline(poly)
                        date_line_fixed = True
                        
                        # If the split resulted in a MultiPolygon, process each part
                        if isinstance(split_geom, MultiPolygon):
                            if verbose:
                                print(f"  Split into {len(split_geom.geoms)} parts")
                            
                            # Process each part of the split polygon
                            for split_idx, split_poly in enumerate(split_geom.geoms):
                                # Extract coordinates for H3 (swapping xy to yx)
                                try:
                                    split_exterior = [(y, x) for x, y in split_poly.exterior.coords]
                                    h3_poly = h3.LatLngPoly(split_exterior)
                                    split_cells = h3.h3shape_to_cells(h3_poly, res=resolution)
                                    cells.extend(split_cells)
                                    if verbose:
                                        print(f"  Successfully converted split part #{split_idx+1}: got {len(split_cells)} cells")
                                except Exception as e:
                                    print(f"  Error converting split part #{split_idx+1}: {e}")
                            
                            # Skip to the next polygon in the MultiPolygon
                            continue
                        # If it's still a Polygon, process it normally below
                        elif isinstance(split_geom, Polygon):
                            poly = split_geom  # Use the fixed polygon
                            if verbose:
                                print(f"  Using split polygon with {len(poly.exterior.coords)} points")
                    except Exception as e:
                        print(f"  Error splitting polygon at date line: {e}")
                        # Continue with original polygon if splitting failed
                
                # Extract exterior and interior coordinates (swapping xy to yx for H3)
                exterior = [(y, x) for x, y in poly.exterior.coords]
                
                if verbose:
                    print(f"  First few coordinates (lat, lng format for H3):\n    {exterior[:3]}")
                    if len(exterior) > 3:
                        print(f"  Last coordinate: {exterior[-1]}")
                
                # Process each polygon individually
                try:
                    # Standard approach - using the polygon directly
                    h3_poly = h3.LatLngPoly(exterior)
                    poly_cells = h3.h3shape_to_cells(h3_poly, res=resolution)
                    cells.extend(poly_cells)
                    if verbose:
                        print(f"  Successfully converted polygon #{poly_idx+1}: got {len(poly_cells)} cells")
                        
                except Exception as e:
                    print(f"  Error converting polygon #{poly_idx+1}: {e}")
                    
                    # Try approaches for problematic polygons
                    approaches = [
                        ("reversed coordinates", lambda: h3.h3shape_to_cells(h3.LatLngPoly(exterior[::-1]), res=resolution)),
                        ("simplified geometry", lambda: h3.h3shape_to_cells(h3.LatLngPoly([(y, x) for x, y in poly.simplify(0.001).exterior.coords]), res=resolution)),
                        ("buffered geometry", lambda: h3.h3shape_to_cells(h3.LatLngPoly([(y, x) for x, y in poly.buffer(0.0001).exterior.coords]), res=resolution)),
                        ("convex hull", lambda: h3.h3shape_to_cells(h3.LatLngPoly([(y, x) for x, y in poly.convex_hull.exterior.coords]), res=resolution))
                    ]
                    
                    for approach_name, approach_func in approaches:
                        try:
                            print(f"  Trying {approach_name}...")
                            poly_cells = approach_func()
                            cells.extend(poly_cells)
                            print(f"  Success with {approach_name}: got {len(poly_cells)} cells")
                            break
                        except Exception as alt_e:
                            print(f"  {approach_name} failed: {alt_e}")
        
        elif geom.geom_type == 'Polygon':
            if verbose:
                print(f"Processing Polygon with {len(geom.exterior.coords)} points")
                print(f"Area: {geom.area:.6f} square degrees")
            
            # Skip polygons with too few points
            if len(geom.exterior.coords) < 4:  # < 4 because first and last points are the same
                print(f"Warning: Skipping polygon with only {len(geom.exterior.coords)} points")
                return cells
            
            # Check if this polygon crosses the international date line
            if check_crosses_idl(geom):
                if verbose:
                    print(f"Polygon crosses the international date line - splitting")
                
                # Split the polygon at the date line
                try:
                    split_geom = split_polygon_at_dateline(geom)
                    date_line_fixed = True
                    
                    # If the split resulted in a MultiPolygon, process each part separately
                    if isinstance(split_geom, MultiPolygon):
                        if verbose:
                            print(f"Split into {len(split_geom.geoms)} parts")
                        
                        all_split_cells = []
                        # Process each part of the split polygon
                        for split_idx, split_poly in enumerate(split_geom.geoms):
                            if verbose:
                                print(f"Processing split part #{split_idx+1}")
                            
                            # Extract coordinates for H3 (swapping xy to yx)
                            try:
                                split_exterior = [(y, x) for x, y in split_poly.exterior.coords]
                                h3_poly = h3.LatLngPoly(split_exterior)
                                split_cells = h3.h3shape_to_cells(h3_poly, res=resolution)
                                all_split_cells.extend(split_cells)
                                if verbose:
                                    print(f"Successfully converted split part #{split_idx+1}: got {len(split_cells)} cells")
                            except Exception as e:
                                print(f"Error converting split part #{split_idx+1}: {e}")
                        
                        # Return the cells from all split parts
                        if all_split_cells:
                            cells = all_split_cells
                            if verbose:
                                print(f"Successfully processed all split parts: got {len(cells)} total cells")
                            return cells
                        
                    # If it's still a Polygon or splitting didn't yield cells, process it normally below
                    elif isinstance(split_geom, Polygon):
                        geom = split_geom  # Use the fixed polygon
                        if verbose:
                            print(f"Using split polygon with {len(geom.exterior.coords)} points")
                except Exception as e:
                    print(f"Error splitting polygon at date line: {e}")
                    # Continue with original polygon if splitting failed
            
            # Handle simple Polygon (swapping xy to yx for H3)
            exterior = [(y, x) for x, y in geom.exterior.coords]
            
            if verbose:
                print(f"First few coordinates (lat, lng format for H3):\n  {exterior[:3]}")
                if len(exterior) > 3:
                    print(f"Last coordinate: {exterior[-1]}")
            
            # Try standard conversion
            try:
                h3_poly = h3.LatLngPoly(exterior)
                cells = h3.h3shape_to_cells(h3_poly, res=resolution)
                if verbose:
                    print(f"Successfully converted polygon: got {len(cells)} cells")
            except Exception as e:
                print(f"Error converting polygon: {e}")
                
                # Try various approaches for problematic polygons
                approaches = [
                    ("reversed coordinates", lambda: h3.h3shape_to_cells(h3.LatLngPoly(exterior[::-1]), res=resolution)),
                    ("simplified geometry", lambda: h3.h3shape_to_cells(h3.LatLngPoly([(y, x) for x, y in geom.simplify(0.001).exterior.coords]), res=resolution)),
                    ("buffered geometry", lambda: h3.h3shape_to_cells(h3.LatLngPoly([(y, x) for x, y in geom.buffer(0.0001).exterior.coords]), res=resolution)),
                    ("convex hull", lambda: h3.h3shape_to_cells(h3.LatLngPoly([(y, x) for x, y in geom.convex_hull.exterior.coords]), res=resolution))
                ]
                
                for approach_name, approach_func in approaches:
                    try:
                        print(f"Trying {approach_name}...")
                        cells = approach_func()
                        print(f"Success with {approach_name}: got {len(cells)} cells")
                        break
                    except Exception as alt_e:
                        print(f"{approach_name} failed: {alt_e}")
                        
                if not cells:
                    print("All approaches failed to convert this polygon!")
        else:
            print(f"Warning: Unsupported geometry type: {geom.geom_type}")
            
    except Exception as e:
        print(f"Unexpected error in convert_geometry_to_h3: {e}")
    
    if verbose:
        print(f"Total cells generated: {len(cells)}")
    
    return cells

def main():
    print("\n## Load and Explore the regions.geojson File")
    # Memory-efficient loading for the GeoJSON file
    try:
        # Try standard loading first
        regions_gdf = gpd.read_file('hplans/regions.geojson')
        print(f"Successfully loaded file with {len(regions_gdf)} features")
    except Exception as e:
        print(f"Standard loading failed with error: {e}")
        print("Attempting memory-efficient loading...")
        
        # Alternative approach using fiona for chunked reading
        import fiona
        
        all_geometries = []
        all_properties = []
        
        # Open the file with fiona and read chunks
        with fiona.open('hplans/regions.geojson') as src:
            total_features = len(src)
            print(f"Total features: {total_features}")
            
            # Read in chunks of 100 features at a time
            chunk_size = 100
            for i in range(0, total_features, chunk_size):
                chunk = list(src.items(i, min(i + chunk_size, total_features)))
                
                # Extract geometries and properties
                for _, feature in chunk:
                    all_geometries.append(shape(feature['geometry']))
                    all_properties.append(feature['properties'])
                
                print(f"Processed features {i} to {min(i + chunk_size, total_features)}")
        
        # Create geodataframe from collected data
        regions_gdf = gpd.GeoDataFrame(
            all_properties, 
            geometry=all_geometries,
            crs=fiona.crs.from_epsg(4326)
        )
        print(f"Successfully created GeoDataFrame with {len(regions_gdf)} features")

    # Display basic information
    print(f"Number of regions: {len(regions_gdf)}")
    print(f"Geometry types: {regions_gdf.geometry.geom_type.unique()}")

    # Check the properties/columns to find region identifiers
    print("\nColumns in the dataset:")
    print(regions_gdf.columns.tolist())

    # Try to find a suitable column for identifying each region
    region_id_column = None
    for col in regions_gdf.columns:
        if col != 'geometry':
            # Check if this column could be a good identifier
            unique_values = regions_gdf[col].nunique()
            if unique_values == len(regions_gdf):
                region_id_column = col
                print(f"\nUsing '{col}' as the region identifier column")
                print(f"Unique regions: {regions_gdf[col].unique().tolist()}")
                break

    # If no unique identifier found create one
    if region_id_column is None:
        # Try to find the best column to use
        best_col = None
        best_count = 0
        for col in regions_gdf.columns:
            if col != 'geometry':
                unique_count = regions_gdf[col].nunique()
                if unique_count > best_count:
                    best_count = unique_count
                    best_col = col
        
        if best_col:
            region_id_column = best_col
            print(f"\nUsing '{best_col}' as region identifier column")
            print(f"Values: {regions_gdf[best_col].unique().tolist()}")
        else:
            # No suitable column create a new one
            regions_gdf['region_id'] = [f"Region_{i}" for i in range(len(regions_gdf))]
            region_id_column = 'region_id'
            print(f"\nCreated 'region_id' column as identifier")
            print(f"Values: {regions_gdf[region_id_column].unique().tolist()}")

    # Plot the regions to visualize them
    plt.figure(figsize=(16, 12))
    plot_df(regions_gdf, column=region_id_column, title="Regions from regions.geojson")
    plt.tight_layout()
    plt.savefig('regions_original.png', bbox_inches='tight', dpi=300)
    print("Saved original regions plot to regions_original.png")

    # Create a directory to store the region plots
    output_dir = 'region_plots'
    os.makedirs(output_dir, exist_ok=True)

    # Plot each region separately
    for i, region_id in enumerate(regions_gdf[region_id_column].unique()):
        region_data = regions_gdf[regions_gdf[region_id_column] == region_id]
        
        plt.figure(figsize=(12, 10))
        plot_df(region_data, title=f"Region: {region_id}")
        plt.tight_layout()
        plt.savefig(f'{output_dir}/region_{i}_{region_id}.png', bbox_inches='tight', dpi=300)
        plt.close()
        
    print(f"Individual region plots saved to {output_dir}/")

    print("\n## Convert Each Region to H3 Hexagons at Different Resolutions")
    # Create a directory for H3 hexagon plots
    h3_output_dir = 'region_h3_plots'
    os.makedirs(h3_output_dir, exist_ok=True)

    # Dictionary to store results
    region_cells = {}

    # Process each region at different resolutions
    resolutions = list(range(2, 5))  # Resolutions 2-4

    for region_idx, region_id in enumerate(regions_gdf[region_id_column].unique()):
        # Get the region data
        region_data = regions_gdf[regions_gdf[region_id_column] == region_id]
        region_cells[region_id] = {}
        
        print(f"\nProcessing region: {region_id}")
        
        # Create a plot with subplots for each resolution
        fig, axes = plt.subplots(len(resolutions) + 1, 1, figsize=(12, 5 * (len(resolutions) + 1)))
        
        # Plot original region
        plot_df(region_data, ax=axes[0], title=f"Original Region: {region_id}")
        
        # Process each resolution
        for i, res in enumerate(resolutions):
            print(f"  Processing resolution {res}...")
            all_cells = []
            
            # Process each feature in the region
            for _, row in tqdm(region_data.iterrows(), total=len(region_data), desc=f"Res {res}"):
                try:
                    cells = convert_geometry_to_h3(row.geometry, resolution=res)
                    all_cells.extend(cells)
                except Exception as e:
                    print(f"Error processing feature: {e}")
                    try:
                        simplified_geom = row.geometry.simplify(0.0001)
                        cells = convert_geometry_to_h3(simplified_geom, resolution=res)
                        all_cells.extend(cells)
                    except Exception as e2:
                        print(f"Simplified conversion also failed: {e2}")
            
            # Save cells for this region and resolution
            all_cells = list(set(all_cells))  # Remove duplicates
            region_cells[region_id][res] = all_cells
            print(f"  Resolution {res}: Generated {len(all_cells)} unique H3 cells")
            
            # Save to JSON file
            os.makedirs(f"{h3_output_dir}/region_{region_idx}_{region_id}", exist_ok=True)
            with open(f"{h3_output_dir}/region_{region_idx}_{region_id}/h3_res{res}.json", 'w') as f:
                json.dump(all_cells, f)
            
            # Plot the cells
            plot_cells(all_cells, ax=axes[i+1], title=f"H3 Hexagons (Resolution {res})")
        
        # Save the combined plot
        plt.tight_layout()
        plt.savefig(f"{h3_output_dir}/region_{region_idx}_{region_id}/all_resolutions.png", bbox_inches='tight', dpi=300)
        plt.close()
        
        # Generate single resolution plots for better detail
        for res in resolutions:
            plt.figure(figsize=(12, 10))
            plot_cells(region_cells[region_id][res], title=f"Region {region_id} - H3 Hexagons (Resolution {res})")
            plt.tight_layout()
            plt.savefig(f"{h3_output_dir}/region_{region_idx}_{region_id}/resolution_{res}.png", bbox_inches='tight', dpi=300)
            plt.close()
        
        # Force garbage collection
        gc.collect()

    print(f"\nAll region H3 hexagon plots saved to {h3_output_dir}/")

    print("\n## Compare Cell Counts for Different Regions")
    # Create a summary table of cell counts
    cell_counts = []
    for region_id in region_cells:
        for res in region_cells[region_id]:
            cell_counts.append({
                'Region': region_id,
                'Resolution': res,
                'Cell Count': len(region_cells[region_id][res])
            })

    cell_count_df = pd.DataFrame(cell_counts)
    print("Summary of H3 Cell Counts by Region and Resolution:")
    summary_table = cell_count_df.pivot_table(values='Cell Count', index='Region', columns='Resolution')
    print(summary_table)

    # Create combined plots for each resolution with all regions in different colors
    print("\n## Creating Combined Region Plots with All Regions in Different Colors")
    
    # For each resolution, create a plot with all regions
    for res in resolutions:
        print(f"Creating combined plot for resolution {res}...")
        plt.figure(figsize=(16, 12))
        
        # Create a GeoDataFrame to hold all H3 cells for this resolution
        all_geometries = []
        all_regions = []
        
        for region_id in region_cells:
            if not region_cells[region_id][res]:
                continue
                
            cells = region_cells[region_id][res]
            
            # For visualization, examine each cell's boundary and completely drop any cell
            # that has boundary points on both sides of the date line
            filtered_cells = []
            problematic_cells = []
            
            for cell in cells:
                # Get all boundary points of the cell
                boundary = h3.cell_to_boundary(cell)
                
                # Check if this cell has points on both sides of the longitude line
                has_positive_lng = False
                has_negative_lng = False
                
                for lat, lng in boundary:
                    if lng > 0:
                        has_positive_lng = True
                    if lng < 0:
                        has_negative_lng = True
                
                # If this cell has boundary points on both sides, drop it completely
                if has_positive_lng and has_negative_lng:
                    problematic_cells.append(cell)
                else:
                    filtered_cells.append(cell)
            
            # Report on filtering
            if problematic_cells:
                print(f"Region {region_id}: Filtered out {len(problematic_cells)} cells near the date line")
                print(f"  Kept {len(filtered_cells)} cells for visualization")
            
            # Only visualize the filtered cells
            if filtered_cells:
                try:
                    # Convert cells to a shape and add to the geometries list
                    shape = h3.cells_to_h3shape(filtered_cells)
                    all_geometries.append(shape)
                    all_regions.append(region_id)
                except Exception as e:
                    print(f"Error creating shape for region {region_id}: {e}")
                    # Fall back to individual cell plotting
                    print(f"Falling back to individual cell plotting for {region_id}")
                    for cell in filtered_cells:
                        try:
                            # Convert each cell to a polygon
                            bounds = h3.cell_to_boundary(cell, geo_json=True)
                            poly = Polygon(bounds)
                            all_geometries.append(poly)
                            all_regions.append(region_id)
                        except Exception as cell_e:
                            print(f"Error with individual cell in {region_id}: {cell_e}")
        
        # Create GeoDataFrame and plot
        if all_geometries:
            combined_gdf = gpd.GeoDataFrame({
                'geometry': all_geometries,
                'region': all_regions
            }, crs='EPSG:4326')
            
            # Plot with different colors for each region
            plot_df(combined_gdf, column='region', title=f"All Regions - H3 Hexagons (Resolution {res})")
            plt.tight_layout()
            plt.savefig(f'combined_regions_res{res}.png', bbox_inches='tight', dpi=300)
            plt.close()
            print(f"Saved combined plot to combined_regions_res{res}.png")
        else:
            print(f"No valid shapes to plot for resolution {res}")
    
    # Create pie charts showing the proportion of cells for each region at each resolution
    print("\n## Creating Pie Charts for Cell Distribution")
    
    # Get a consistent color map to match the region plots
    from matplotlib.colors import ListedColormap
    import matplotlib.cm as cm
    
    # Filter out regions to exclude (CN470, EU433, CD900-1A)
    excluded_regions = ['CN470', 'EU433', 'CD900-1A']
    print(f"Excluding regions: {', '.join(excluded_regions)}")
    
    # Generate colors using the same categorical colormap that GeoPandas uses by default
    region_ids = [r for r in region_cells.keys() if r not in excluded_regions]
    num_regions = len(region_ids)
    color_map = plt.get_cmap('tab10' if num_regions <= 10 else 'tab20')
    
    # Create a color mapping for regions with "Unknown" regions fixed as red
    colors = []
    for i, r in enumerate(region_ids):
        if r == 'Unknown':
            colors.append('red')  # Make Unknown region red
        else:
            colors.append(color_map(i % color_map.N))
    
    # Create a pie chart for each resolution
    for res in resolutions:
        plt.figure(figsize=(14, 10))
        
        # Calculate cell counts for each region at this resolution
        counts = []
        labels = []
        
        total_cells = 0
        for i, region_id in enumerate(region_ids):
            count = len(region_cells[region_id][res])
            if count > 0:  # Only include regions with cells
                counts.append(count)
                # Use shorter label format to avoid overlapping
                labels.append(f"{region_id}")
                total_cells += count
        
        # Create the pie chart with better label placement
        # Move labels outside the pie with connecting lines
        wedges, texts, autotexts = plt.pie(
            counts,
            labels=None,  # No direct labels - will use a legend instead
            colors=colors[:len(counts)],
            autopct='%1.1f%%',
            startangle=90,
            pctdistance=0.85,  # Move percentage labels inward
            wedgeprops={'edgecolor': 'w', 'linewidth': 1},
            explode=[0.05] * len(counts)  # Slight explosion for all slices
        )
        
        # Set properties for percentage labels
        plt.setp(autotexts, size=9, weight="bold", color="white")
        
        # Add a legend outside the pie chart for better readability
        # Calculate percentage for each region to show in legend
        legend_labels = [f"{labels[i]} ({counts[i]:,} cells, {counts[i]/total_cells*100:.1f}%)" 
                         for i in range(len(labels))]
        
        # Place legend to the right of the chart with appropriate font size
        plt.legend(wedges, legend_labels, 
                  title="Regions",
                  loc="center left", 
                  bbox_to_anchor=(1, 0, 0.5, 1),
                  fontsize=9)
        
        plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle
        
        # Add a title with the total cell count
        plt.suptitle(f"H3 Cell Distribution - Resolution {res}", fontsize=16, y=0.98)
        plt.title(f"Total: {total_cells:,} cells", fontsize=12)
        
        # Save the pie chart
        plt.tight_layout()
        plt.savefig(f'cell_distribution_res{res}_pie.png', bbox_inches='tight', dpi=300)
        plt.close()
        print(f"Saved pie chart for resolution {res}")
    
    # Also create a summary bar chart showing total cell counts by resolution
    plt.figure(figsize=(12, 8))
    
    # Calculate total cells per resolution
    total_by_res = []
    for res in resolutions:
        total = sum(len(region_cells[region_id][res]) for region_id in region_cells)
        total_by_res.append(total)
    
    # Create bars
    bars = plt.bar(
        [f"Resolution {res}" for res in resolutions],
        total_by_res,
        color='steelblue'
    )
    
    # Add the total count on top of each bar
    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width()/2.,
            height + 0.05 * max(total_by_res),
            f'{int(height):,}',
            ha='center', 
            va='bottom',
            fontweight='bold'
        )
    
    plt.title("Total H3 Cell Count by Resolution")
    plt.ylabel("Number of Cells")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig('total_cells_by_resolution.png', bbox_inches='tight', dpi=300)
    print("Saved total cell count chart to total_cells_by_resolution.png")

if __name__ == "__main__":
    main()
