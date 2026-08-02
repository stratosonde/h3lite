#!/usr/bin/env python3
"""
H3Lite Lookup Table Generator

Generates the packed, resolution-aware region lookup table for the H3Lite
library from LoRaWAN region GeoJSON files.

H3-1/H3-4: entries are packed uint32_t carrying (baseCell, resolution,
partialIndex, regionId); the device probes res 3 -> res 2 -> res 1,
most-specific-first. The compaction below is the legacy region-blind one
(H3-2 replaces it with uniform-only compaction).

Dependencies: h3, shapely, tqdm   (no geopandas/numpy)
"""

import sys
import json
import os
import argparse
import h3
from shapely.geometry import shape
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
    "AS923-1B": 5,
    "AS923-1C": 6,
    "AS923-2": 7,
    "AS923-3": 8,
    "AS923-4": 9,
    "KR920": 10,
    "IN865": 11,
    "RU864": 12,
    "CN470": 13,
    "EU433": 14,
    "CD900-1A": 15,
    # Add other regions as needed
}


def load_regions(directory=".."):
    """Load all GeoJSON files that match region names.

    Returns {region_name: [shapely geometry, ...]}.
    """
    regions = {}
    files = [f for f in os.listdir(directory) if f.endswith('.geojson')]

    for filename in files:
        region_name = os.path.splitext(filename)[0]

        if region_name not in REGION_IDS:
            # Substring match, longest keys first so "AS923-1B" beats "AS923-1"
            found_match = False
            for key in sorted(REGION_IDS.keys(), key=len, reverse=True):
                if key in region_name:
                    region_name = key
                    found_match = True
                    break
            if not found_match:
                continue

        try:
            with open(os.path.join(directory, filename), encoding='utf-8') as f:
                gj = json.load(f)
            geoms = []
            if gj.get('type') == 'FeatureCollection':
                for feat in gj.get('features', []):
                    g = shape(feat['geometry'])
                    if not g.is_empty:
                        geoms.append(g)
            elif gj.get('type') == 'Feature':
                geoms.append(shape(gj['geometry']))
            else:
                geoms.append(shape(gj))
            regions[region_name] = geoms
            print(f"Loaded region {region_name} from {filename} "
                  f"({len(geoms)} geometries)")
        except Exception as e:
            print(f"Error loading {filename}: {e}")

    return regions


def _polygons(geom):
    """Yield individual Polygon parts of a (Multi)Polygon geometry."""
    if geom.geom_type == 'MultiPolygon':
        yield from geom.geoms
    elif geom.geom_type == 'Polygon':
        yield geom
    else:
        print(f"Warning: Unsupported geometry type: {geom.geom_type}")


def convert_region_to_h3(geoms, resolution=DEFAULT_RESOLUTION):
    """Convert a list of shapely geometries to H3 cells at the resolution."""
    all_h3_cells = []

    for geom in tqdm(geoms, desc="Converting to H3"):
        for poly in _polygons(geom):
            if len(poly.exterior.coords) < 4:
                print(f"Warning: Skipping polygon with only "
                      f"{len(poly.exterior.coords)} points")
                continue

            # Convert (lng, lat) to (lat, lng) for H3
            exterior = [(y, x) for x, y in poly.exterior.coords]

            cells = None
            try:
                cells = h3.h3shape_to_cells(h3.LatLngPoly(exterior),
                                            res=resolution)
            except Exception as e:
                print(f"Error converting polygon: {e}")
                approaches = [
                    ("reversed coordinates",
                     lambda: h3.h3shape_to_cells(
                         h3.LatLngPoly(exterior[::-1]), res=resolution)),
                    ("simplified geometry",
                     lambda: h3.h3shape_to_cells(h3.LatLngPoly(
                         [(y, x) for x, y in
                          poly.simplify(0.001).exterior.coords]),
                         res=resolution)),
                    ("buffered geometry",
                     lambda: h3.h3shape_to_cells(h3.LatLngPoly(
                         [(y, x) for x, y in
                          poly.buffer(0.0001).exterior.coords]),
                         res=resolution)),
                    ("convex hull",
                     lambda: h3.h3shape_to_cells(h3.LatLngPoly(
                         [(y, x) for x, y in
                          poly.convex_hull.exterior.coords]),
                         res=resolution)),
                ]
                for approach_name, approach_func in approaches:
                    try:
                        print(f"Trying {approach_name}...")
                        cells = approach_func()
                        print(f"Success with {approach_name}: "
                              f"got {len(cells)} cells")
                        break
                    except Exception as alt_e:
                        print(f"{approach_name} failed: {alt_e}")
                if not cells:
                    print("All approaches failed to convert this polygon!")

            if cells:
                all_h3_cells.extend(cells)

    all_h3_cells = list(set(all_h3_cells))
    print(f"Generated {len(all_h3_cells)} unique H3 cells "
          f"at resolution {resolution}")
    return all_h3_cells


def extract_h3_components(h3_index):
    """Extract (base_cell, partial_index, resolution) from an H3 index.

    partial_index uses min(3, res) digits: 1 digit at res 1 (0-7),
    2 digits at res 2 (0-63), 3 digits at res 3 (0-511).
    """
    h3_int = int(h3_index, 16)
    res = (h3_int >> 52) & 0xF
    base_cell = (h3_int >> 45) & 0x7F

    partial_index = 0
    num_digits = min(3, res)
    for r in range(1, num_digits + 1):
        shift = 45 - (3 * r)
        digit = (h3_int >> shift) & 0x7
        partial_index = (partial_index * 8) + digit

    return base_cell, partial_index, res


def pack_entry(base_cell, resolution, partial_index, region_id):
    """Pack an entry into the uint32 format (H3-1/H3-4)."""
    return ((base_cell & 0x7F) << 25 | (resolution & 0x03) << 23 |
            (partial_index & 0x1FF) << 14 | (region_id & 0x0F) << 10)


def generate_lookup_table(regions, max_resolution=DEFAULT_RESOLUTION):
    """Legacy multi-resolution compaction (H3-2 replaces this).

    Returns entries with keys: baseCell, resolution, partialIndex,
    regionId, regionName.
    """
    entries = []
    processed_cells = set()

    for resolution in range(1, max_resolution + 1):
        print(f"\nProcessing resolution {resolution}...")

        for region_name, geoms in regions.items():
            if region_name not in REGION_IDS:
                print(f"Skipping region {region_name} - not in REGION_IDS")
                continue

            region_id = REGION_IDS[region_name]
            print(f"Processing region {region_name} (ID: {region_id}) "
                  f"at resolution {resolution}")

            h3_cells = convert_region_to_h3(geoms, resolution)

            for cell in h3_cells:
                base_cell, partial_index, _ = extract_h3_components(cell)
                cell_key = (base_cell, partial_index, resolution)

                if resolution > 1:
                    # Suppress children whose parent was already emitted
                    parent_partial_index = partial_index // 8
                    if (base_cell, parent_partial_index,
                            resolution - 1) in processed_cells:
                        continue

                if cell_key not in processed_cells:
                    processed_cells.add(cell_key)
                    entries.append({
                        'h3': cell,
                        'baseCell': base_cell,
                        'resolution': resolution,
                        'partialIndex': partial_index,
                        'regionId': region_id,
                        'regionName': region_name,
                    })

    print(f"\nGenerated {len(entries)} lookup table entries")
    return entries


def generate_c_code(entries, output_c_file, output_h_file):
    """Generate C code for the packed lookup table."""
    # Numeric sort of the packed value == tuple sort by
    # (baseCell, resolution, partialIndex, regionId)
    entries.sort(key=lambda x: pack_entry(x['baseCell'], x['resolution'],
                                          x['partialIndex'], x['regionId']))

    max_id = max(REGION_IDS.values())

    # ---- header ----
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
        for name, region_id in REGION_IDS.items():
            f.write(f"#define REGION_{name.replace('-', '_')} {region_id}\n")

        f.write(f"""
// Number of entries in the lookup table
#define REGION_ENTRY_COUNT {len(entries)}

/* Packed entry (H3-1/H3-4):
 *   bits [31:25] baseCell (0-121)
 *   bits [24:23] resolution (1-3)
 *   bits [22:14] partialIndex (0-511)
 *   bits [13:10] regionId (0-{max_id})
 *   bits [9:0]   reserved (zero)
 * Sorted ascending by (baseCell, resolution, partialIndex) — which is
 * numeric order of the packed value. */
typedef uint32_t RegionEntry;

#define RE_BASECELL(e)  ((uint8_t)(((e) >> 25) & 0x7F))
#define RE_RES(e)       ((uint8_t)(((e) >> 23) & 0x03))
#define RE_PARTIAL(e)   ((uint16_t)(((e) >> 14) & 0x1FF))
#define RE_REGION(e)    ((RegionId)(((e) >> 10) & 0x0F))

// Region lookup table (sorted by baseCell, resolution, partialIndex)
extern const RegionEntry regionLookup[REGION_ENTRY_COUNT];

// Region names array
extern const char* regionNames[{max_id + 1}];

// Binary search the region lookup table for an exact
// (baseCell, resolution, partialIndex) key; returns 0 (Unknown) on miss.
RegionId findRegion(uint8_t baseCell, uint8_t res, uint16_t partialIndex);

#endif /* H3LITE_REGIONS_TABLE_H */
""")

    # ---- C file ----
    with open(output_c_file, 'w') as f:
        f.write(f"""/*
 * H3Lite Region Lookup Table Implementation
 * Generated by generate_lookup_table.py
 */

#include "../include/h3lite_regions_table.h"

// Region names array
const char* regionNames[{max_id + 1}] = {{
    "Unknown", // ID 0
""")
        for name, region_id in sorted(REGION_IDS.items(), key=lambda x: x[1]):
            f.write(f'    "{name}", // ID {region_id}\n')

        f.write(f"""
}};

// Region lookup table (packed uint32, sorted by packed value)
const RegionEntry regionLookup[REGION_ENTRY_COUNT] = {{
""")
        for e in entries:
            packed = pack_entry(e['baseCell'], e['resolution'],
                                e['partialIndex'], e['regionId'])
            f.write(f"    0x{packed:08X}u, // bc={e['baseCell']} "
                    f"res={e['resolution']} pi={e['partialIndex']} "
                    f"{e['regionName']}\n")

        f.write("""
};

/* Binary search on the (baseCell, resolution, partialIndex) key.
 * The packed layout puts those fields in the high bits, so masking off
 * the low 14 bits (regionId + reserved) leaves a directly comparable key. */
#define RE_KEY_MASK 0xFFFFC000u

RegionId findRegion(uint8_t baseCell, uint8_t res, uint16_t partialIndex) {
    uint32_t probe = ((uint32_t)baseCell << 25) |
                     ((uint32_t)res << 23) |
                     ((uint32_t)partialIndex << 14);
    int low = 0;
    int high = REGION_ENTRY_COUNT - 1;

    while (low <= high) {
        int mid = (low + high) / 2;
        uint32_t key = regionLookup[mid] & RE_KEY_MASK;

        if (key < probe) {
            low = mid + 1;
        } else if (key > probe) {
            high = mid - 1;
        } else {
            return RE_REGION(regionLookup[mid]);
        }
    }

    return 0; // Unknown/invalid region
}
""")

    print(f"Generated C code in {output_c_file} and {output_h_file}")
    print(f"Lookup table size: {len(entries) * 4} bytes "
          f"({len(entries)} packed entries)")
    print(f"Number of unique regions: "
          f"{len(set(e['regionId'] for e in entries))}")

    cells_by_region = {}
    for e in entries:
        cells_by_region[e['regionName']] = \
            cells_by_region.get(e['regionName'], 0) + 1
    print("\nCells per region:")
    for region, count in sorted(cells_by_region.items(),
                                key=lambda x: x[1], reverse=True):
        print(f"  {region}: {count} cells")


def main():
    parser = argparse.ArgumentParser(
        description='Generate H3Lite lookup tables from GeoJSON region files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 generate_lookup_table.py --geojson-dir ../hplans
  python3 generate_lookup_table.py --geojson-dir /path/to/regions
        """
    )
    parser.add_argument('--geojson-dir', type=str, default='../hplans',
                        help='Directory containing GeoJSON region files '
                             '(default: ../hplans)')
    parser.add_argument('--resolution', type=int, default=DEFAULT_RESOLUTION,
                        help=f'H3 resolution level (default: '
                             f'{DEFAULT_RESOLUTION})')
    args = parser.parse_args()

    print(f"Loading region definitions from GeoJSON files in "
          f"{args.geojson_dir}...")
    regions = load_regions(directory=args.geojson_dir)

    if not regions:
        print(f"No region files found or loaded in {args.geojson_dir}. "
              f"Exiting.")
        return

    print(f"\nGenerating lookup tables at resolution {args.resolution}...")
    entries = generate_lookup_table(regions, max_resolution=args.resolution)

    print("\nGenerating C code...")
    generate_c_code(entries,
                    os.path.join('src', OUTPUT_C_FILE),
                    os.path.join('include', OUTPUT_H_FILE))

    print("\nLookup table generation complete!")


if __name__ == "__main__":
    main()
