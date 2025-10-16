#!/usr/bin/env python3
"""
H3Lite Lookup Table Generator

This script generates a compact lookup table for the H3Lite library 
based on LoRaWAN region definitions from GeoJSON files.

It:
1. Loads region boundary definitions from GeoJSON files
2. Converts each region to H3 cells at the specified resolution
3. Generates a compact lookup table for use in the H3Lite C implementation
4. Outputs C code for the lookup table

This is optimized for STM32 microcontrollers with limited memory
"""

import sys
import json
import os
import h3
import geopandas as gpd
from shapely.geometry import shape
import numpy as np
from tqdm import tqdm

# Configuration
DEFAULT_RESOLUTION = 3  # Resolution 3 gives a good balance of accuracy vs size
OUTPUT_C_FILE = "h3lite_regions_table.c"
OUTPUT_H_FILE = "h3lite_regions_table.h"

# Region ID mapping (customize these for your application)
REGION_IDS = {
    "US915": 1,
    "EU868": 2,
    "AU915": 3,
    "AS923-1": 4,
    "AS923-2": 5,
    "AS923-3": 6,
    "AS923-4": 7,
    "KR920": 8,
    "IN865": 9,
    "RU864": 10,
    "CN470": 11,
    "EU433": 12,
    # Add other regions as needed
}

def load_regions(directory=".."):  # Look in parent directory for GeoJSON files
    """Load all GeoJSON files that match region names from the current directory"""
    regions = {}
    files = [f for f in os.listdir(directory) if f.endswith('.geojson')]
    
    for filename in files:
        # Extract region name from filename (strip .geojson extension)
        region_name = os.path.splitext(filename)[0]
        
        # Check if this is a region we're interested in
        found_match = False
        for key in REGION_IDS.keys():
            if key in region_name:
                region_name = key
                found_match = True
                break
        
        if not found_match:
            continue
        
        # Load the GeoJSON
        try:
            gdf = gpd.read_file(os.path.join(directory, filename))
            regions[region_name] = gdf
            print(f"Loaded region {region_name} from {filename}")
        except Exception as e:
            print(f"Error loading {filename}: {e}")
    
    return regions

def convert_region_to_h3(gdf, resolution=DEFAULT_RESOLUTION):
    """Convert a GeoDataFrame to H3 cells at the specified resolution
    Based on the implementation in regions_h3_conversion.py
    """
    all_h3_cells = []
    
    for _, row in tqdm(gdf.iterrows(), desc="Converting to H3", total=len(gdf)):
        geometry = row.geometry
        
        try:
            if geometry.geom_type == 'MultiPolygon':
                # Handle MultiPolygon by processing each polygon separately
                for poly in geometry.geoms:
                    # Skip polygons with too few points
                    if len(poly.exterior.coords) < 4:
                        print(f"Warning: Skipping polygon with only {len(poly.exterior.coords)} points")
                        continue
                    
                    # Convert (lng, lat) to (lat, lng) for H3
                    exterior = [(y, x) for x, y in poly.exterior.coords]
                    
                    try:
                        # Standard approach
                        h3_poly = h3.LatLngPoly(exterior)
                        cells = h3.h3shape_to_cells(h3_poly, res=resolution)
                        all_h3_cells.extend(cells)
                    except Exception as e:
                        print(f"Error converting polygon: {e}")
                        
                        # Try fallback approaches for problematic polygons
                        approaches = [
                            ("reversed coordinates", lambda: h3.h3shape_to_cells(h3.LatLngPoly(exterior[::-1]), res=resolution)),
                            ("simplified geometry", lambda: h3.h3shape_to_cells(h3.LatLngPoly(
                                [(y, x) for x, y in poly.simplify(0.001).exterior.coords]), res=resolution)),
                            ("buffered geometry", lambda: h3.h3shape_to_cells(h3.LatLngPoly(
                                [(y, x) for x, y in poly.buffer(0.0001).exterior.coords]), res=resolution)),
                            ("convex hull", lambda: h3.h3shape_to_cells(h3.LatLngPoly(
                                [(y, x) for x, y in poly.convex_hull.exterior.coords]), res=resolution))
                        ]
                        
                        for approach_name, approach_func in approaches:
                            try:
                                print(f"Trying {approach_name}...")
                                cells = approach_func()
                                all_h3_cells.extend(cells)
                                print(f"Success with {approach_name}: got {len(cells)} cells")
                                break
                            except Exception as alt_e:
                                print(f"{approach_name} failed: {alt_e}")
            
            elif geometry.geom_type == 'Polygon':
                # Skip polygons with too few points
                if len(geometry.exterior.coords) < 4:
                    print(f"Warning: Skipping polygon with only {len(geometry.exterior.coords)} points")
                    continue
                
                # Convert (lng, lat) to (lat, lng) for H3
                exterior = [(y, x) for x, y in geometry.exterior.coords]
                
                try:
                    # Standard approach
                    h3_poly = h3.LatLngPoly(exterior)
                    cells = h3.h3shape_to_cells(h3_poly, res=resolution)
                    all_h3_cells.extend(cells)
                except Exception as e:
                    print(f"Error converting polygon: {e}")
                    
                    # Try fallback approaches for problematic polygons
                    approaches = [
                        ("reversed coordinates", lambda: h3.h3shape_to_cells(h3.LatLngPoly(exterior[::-1]), res=resolution)),
                        ("simplified geometry", lambda: h3.h3shape_to_cells(h3.LatLngPoly(
                            [(y, x) for x, y in geometry.simplify(0.001).exterior.coords]), res=resolution)),
                        ("buffered geometry", lambda: h3.h3shape_to_cells(h3.LatLngPoly(
                            [(y, x) for x, y in geometry.buffer(0.0001).exterior.coords]), res=resolution)),
                        ("convex hull", lambda: h3.h3shape_to_cells(h3.LatLngPoly(
                            [(y, x) for x, y in geometry.convex_hull.exterior.coords]), res=resolution))
                    ]
                    
                    for approach_name, approach_func in approaches:
                        try:
                            print(f"Trying {approach_name}...")
                            cells = approach_func()
                            all_h3_cells.extend(cells)
                            print(f"Success with {approach_name}: got {len(cells)} cells")
                            break
                        except Exception as alt_e:
                            print(f"{approach_name} failed: {alt_e}")
                    
                    if not cells:
                        print("All approaches failed to convert this polygon!")
            
            else:
                print(f"Warning: Unsupported geometry type: {geometry.geom_type}")
        
        except Exception as e:
            print(f"Unexpected error processing geometry: {e}")
    
    # Remove duplicates
    all_h3_cells = list(set(all_h3_cells))
    print(f"Generated {len(all_h3_cells)} unique H3 cells at resolution {resolution}")
    
    return all_h3_cells

def extract_h3_components(h3_index):
    """Extract base cell and partial index from an H3 index"""
    # Convert h3_index to integer
    h3_int = int(h3_index, 16)
    
    # Extract resolution (using bit operations instead of API calls)
    res = (h3_int >> 52) & 0xF
    
    # Extract base cell (bits 45-52)
    base_cell = (h3_int >> 45) & 0x7F
    
    # Extract the first 3 resolution digits directly from the bit pattern
    # Each digit in H3 uses 3 bits, starting from bit 45 and going down
    partial_index = 0
    for r in range(1, min(4, res+1)):
        # Calculate position in the bit pattern based on resolution
        # First digit is at bits 42-44, second at 39-41, third at 36-38
        shift = 45 - (3 * r)
        digit = (h3_int >> shift) & 0x7
        partial_index = (partial_index * 8) + digit
    
    # For debugging
    print(f"H3: {h3_index}, BC: {base_cell}, PI: {partial_index}")
    
    return base_cell, partial_index

def generate_lookup_table(regions, max_resolution=DEFAULT_RESOLUTION):
    """Generate lookup tables for all regions using a multi-resolution approach"""
    entries = []
    processed_cells = set()  # Keep track of areas we've already covered
    
    # Process resolutions from lowest to highest for better coverage with fewer cells
    for resolution in range(1, max_resolution + 1):
        print(f"\nProcessing resolution {resolution}...")
        
        for region_name, gdf in regions.items():
            if region_name not in REGION_IDS:
                print(f"Skipping region {region_name} - not in REGION_IDS")
                continue
                
            region_id = REGION_IDS[region_name]
            print(f"Processing region {region_name} (ID: {region_id}) at resolution {resolution}")
            
            # For resolution 1, process all regions completely
            if resolution == 1:
                h3_cells = convert_region_to_h3(gdf, resolution)
                
                for cell in h3_cells:
                    base_cell, partial_index = extract_h3_components(cell)
                    cell_key = (base_cell, partial_index, resolution)
                    
                    if cell_key not in processed_cells:
                        processed_cells.add(cell_key)
                        entries.append({
                            'h3': cell,
                            'baseCell': base_cell,
                            'partialIndex': partial_index,
                            'regionId': region_id,
                            'regionName': region_name
                        })
            
            # For higher resolutions, only process areas near region boundaries or not covered by lower resolutions
            else:
                # Convert region to H3 cells at current resolution
                h3_cells = convert_region_to_h3(gdf, resolution)
                
                for cell in h3_cells:
                    base_cell, partial_index = extract_h3_components(cell)
                    cell_key = (base_cell, partial_index, resolution)
                    
                    # Check if this cell or its parent is already processed
                    parent_already_processed = False
                    
                    # Calculate parent's partial index (divide by 8 for each level up)
                    parent_partial_index = partial_index // 8
                    if (base_cell, parent_partial_index, resolution-1) in processed_cells:
                        parent_already_processed = True
                    
                    # Only add cells that aren't already covered by a parent
                    if not parent_already_processed and cell_key not in processed_cells:
                        processed_cells.add(cell_key)
                        entries.append({
                            'h3': cell,
                            'baseCell': base_cell,
                            'partialIndex': partial_index,
                            'regionId': region_id,
                            'regionName': region_name
                        })
    
    print(f"Generated {len(entries)} lookup table entries")
    return entries

def generate_c_code(entries, output_c_file, output_h_file):
    """Generate C code for the lookup table"""
    # Sort by base cell and partial index for binary search
    entries.sort(key=lambda x: (x['baseCell'], x['partialIndex']))
    
    # Generate header file
    with open(output_h_file, 'w') as f:
        f.write("""/*
 * H3Lite Region Lookup Table
 * Generated by generate_lookup_table.py
 */

#ifndef H3LITE_REGIONS_TABLE_H
#define H3LITE_REGIONS_TABLE_H

#include <stdint.h>
#include "h3lite.h"

// Region IDs
""")
        
        # Write region ID definitions
        for name, region_id in REGION_IDS.items():
            f.write(f"#define REGION_{name.replace('-', '_')} {region_id}\n")
        
        f.write(f"""
// Number of entries in the lookup table
#define REGION_ENTRY_COUNT {len(entries)}

// Lookup entry structure
typedef struct {{
    uint8_t baseCell;      // Base cell (0-121)
    uint16_t partialIndex; // Partial index for first few resolutions
    RegionId regionId;     // Region ID (1-{max(REGION_IDS.values())})
}} RegionEntry;

// Region lookup table (sorted by baseCell and partialIndex)
extern const RegionEntry regionLookup[REGION_ENTRY_COUNT];

// Region names array
extern const char* regionNames[{max(REGION_IDS.values()) + 1}];

// Binary search the region lookup table
RegionId findRegion(uint8_t baseCell, uint16_t partialIndex);

#endif /* H3LITE_REGIONS_TABLE_H */
""")
    
    # Generate C file
    with open(output_c_file, 'w') as f:
        f.write(f"""/*
 * H3Lite Region Lookup Table Implementation
 * Generated by generate_lookup_table.py
 */

#include "../include/h3lite_regions_table.h"

// Region names array
const char* regionNames[{max(REGION_IDS.values()) + 1}] = {{
    "Unknown", // ID 0
""")
        
        # Write region names
        for name, region_id in sorted(REGION_IDS.items(), key=lambda x: x[1]):
            f.write(f'    "{name}", // ID {region_id}\n')
        
        f.write(f"""
}};

// Region lookup table
const RegionEntry regionLookup[REGION_ENTRY_COUNT] = {{
""")
        
        # Write lookup table entries
        for entry in entries:
            f.write(f'    {{ {entry["baseCell"]}, {entry["partialIndex"]}, {entry["regionId"]} }}, // {entry["regionName"]}\n')
        
        f.write(f"""
}};

// Binary search the region lookup table
RegionId findRegion(uint8_t baseCell, uint16_t partialIndex) {{
    int low = 0;
    int high = REGION_ENTRY_COUNT - 1;
    
    while (low <= high) {{
        int mid = (low + high) / 2;
        
        // Check base cell first (more significant)
        if (regionLookup[mid].baseCell < baseCell) {{
            low = mid + 1;
        }} else if (regionLookup[mid].baseCell > baseCell) {{
            high = mid - 1;
        }} else {{
            // Base cells match, check partial index
            if (regionLookup[mid].partialIndex < partialIndex) {{
                low = mid + 1;
            }} else if (regionLookup[mid].partialIndex > partialIndex) {{
                high = mid - 1;
            }} else {{
                // Found a match
                return regionLookup[mid].regionId;
            }}
        }}
    }}
    
    // No match found
    return 0; // Unknown/invalid region
}}
""")
    
    print(f"Generated C code in {output_c_file} and {output_h_file}")
    
    # Print some stats
    table_size = len(entries) * (1 + 2 + 1)  # baseCell + partialIndex + regionId
    print(f"Lookup table size: {table_size} bytes")
    print(f"Number of unique regions: {len(set(e['regionId'] for e in entries))}")
    
    # Count cells by region
    cells_by_region = {}
    for entry in entries:
        region = entry['regionName']
        cells_by_region[region] = cells_by_region.get(region, 0) + 1
    
    print("\nCells per region:")
    for region, count in sorted(cells_by_region.items(), key=lambda x: x[1], reverse=True):
        print(f"  {region}: {count} cells")

def main():
    # Load regions
    print("Loading region definitions from GeoJSON files...")
    regions = load_regions()
    
    if not regions:
        print("No region files found or loaded. Exiting.")
        return
    
    # Generate lookup tables
    print(f"\nGenerating lookup tables at resolution {DEFAULT_RESOLUTION}...")
    entries = generate_lookup_table(regions)
    
    # Generate C code
    print("\nGenerating C code...")
    generate_c_code(entries, 
                   os.path.join('src', OUTPUT_C_FILE),
                   os.path.join('include', OUTPUT_H_FILE))
    
    print("\nLookup table generation complete!")

if __name__ == "__main__":
    main()
