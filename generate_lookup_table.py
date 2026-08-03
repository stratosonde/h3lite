#!/usr/bin/env python3
"""
H3Lite Lookup Table Generator

Generates the packed, resolution-aware region lookup table for the H3Lite
library from LoRaWAN region GeoJSON files.

H3-1/H3-4: entries are packed uint32_t carrying (baseCell, resolution,
partialIndex, regionId); the device probes res 3 -> res 2 -> res 1,
most-specific-first.
H3-2: compaction is uniform-only — a parent is emitted only when all of
its children are present and map to the same region; contested res-3
cells are resolved deterministically (largest intersected area wins,
lowest region ID breaks ties).

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
    # ID 15 is RESTRICTED (no transmission allowed). It was repurposed from
    # the CD900-1A test plan, which is not used in production. The packed
    # entry gives regionId only 4 bits, so 15 is the largest encodable ID;
    # REGION_RESTRICTED can therefore never be 255 in the table.
    "RESTRICTED": 15,
    # Add other regions as needed
}


def load_regions(directory):
    """Load all GeoJSON files that match region names.

    Returns {region_name: [shapely geometry, ...]}.
    """
    regions = {}
    # Sorted for determinism: os.listdir() order is filesystem-dependent,
    # and two files mapping to the same region now MERGE (setdefault)
    # instead of the later one silently overwriting the earlier one.
    files = sorted(f for f in os.listdir(directory) if f.endswith('.geojson'))

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
            regions.setdefault(region_name, []).extend(geoms)
            print(f"Loaded region {region_name} from {filename} "
                  f"({len(geoms)} geometries)")
        except Exception as e:
            print(f"Error loading {filename}: {e}")

    return regions


def _crosses_antimeridian(poly):
    """True if a polygon's longitude span is only explicable by wrapping.

    A polygon whose exterior spans more than 180 degrees of longitude is
    almost certainly a dateline-crossing shape stored as -179..+179 rather
    than a genuinely hemispheric one. h3shape_to_cells has no notion of
    wrapping and will fill the *complement* of such a shape, which is how
    Fiji, Tuvalu and Kiribati ended up entirely absent from the table.
    """
    lons = [x for x, _ in poly.exterior.coords]
    return (max(lons) - min(lons)) > 180.0


def _split_at_antimeridian(poly):
    """Split a dateline-crossing polygon into eastern and western parts."""
    from shapely.geometry import box
    from shapely.ops import unary_union

    # Shift western-hemisphere vertices into 0..360, clip, then shift back.
    shifted = type(poly)(
        [(x + 360.0 if x < 0 else x, y) for x, y in poly.exterior.coords],
        [[(x + 360.0 if x < 0 else x, y) for x, y in ring.coords]
         for ring in poly.interiors])
    if not shifted.is_valid:
        shifted = shifted.buffer(0)

    east = shifted.intersection(box(-180.0, -90.0, 180.0, 90.0))
    west = shifted.intersection(box(180.0, -90.0, 540.0, 90.0))
    parts = []
    if not east.is_empty:
        parts.append(east)
    if not west.is_empty:
        # translate back into -180..0
        from shapely.affinity import translate
        parts.append(translate(west, xoff=-360.0))
    return unary_union(parts) if parts else poly


def _polygons(geom):
    """Yield individual Polygon parts of a (Multi)Polygon geometry.

    Dateline-crossing parts are split so that each yielded polygon lies
    wholly within -180..180 without wrapping.
    """
    if geom.geom_type == 'MultiPolygon':
        parts = list(geom.geoms)
    elif geom.geom_type == 'Polygon':
        parts = [geom]
    else:
        print(f"Warning: Unsupported geometry type: {geom.geom_type}")
        return

    for part in parts:
        if _crosses_antimeridian(part):
            split = _split_at_antimeridian(part)
            if split.geom_type == 'MultiPolygon':
                yield from split.geoms
            else:
                yield split
        else:
            yield part


def convert_region_to_h3(geoms, resolution=DEFAULT_RESOLUTION):
    """Convert a list of shapely geometries to H3 cells at the resolution."""
    all_h3_cells = []

    for geom in tqdm(geoms, desc="Converting to H3"):
        for poly in _polygons(geom):
            if len(poly.exterior.coords) < 4:
                print(f"Warning: Skipping polygon with only "
                      f"{len(poly.exterior.coords)} points")
                continue

            # Convert (lng, lat) to (lat, lng) for H3, honouring holes.
            # Dropping interior rings fills enclaves solid (Lesotho inside
            # South Africa, Vatican inside Italy, ...).
            exterior = [(y, x) for x, y in poly.exterior.coords]
            holes = [[(y, x) for x, y in ring.coords]
                     for ring in poly.interiors]

            cells = None
            try:
                cells = h3.h3shape_to_cells(h3.LatLngPoly(exterior, *holes),
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
                    # NO convex-hull fallback: a hull can cover territory
                    # that was never in the source polygon. Silently
                    # enlarging a regulatory boundary is worse than
                    # failing — so we fail (below) instead.
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
                    # Fail closed: a region whose geometry cannot be
                    # converted must abort the whole generation rather
                    # than silently ship a table with missing territory.
                    raise RuntimeError(
                        "All repair approaches failed to convert this "
                        "polygon; aborting generation. Fix the source "
                        "GeoJSON rather than shipping a partial table.")

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
    """Pack an entry into the uint32 format (H3-1/H3-4).

    Raises rather than masking: the field widths are tight (regionId is
    4 bits and IDs 1-15 are already allocated), so a silent `& 0x0F` would
    turn a newly added region 16 into 0 = Unknown, and `& 0x03` would turn
    a res-4 table into res 0.
    """
    if not 0 <= base_cell <= 121:
        raise ValueError(f"baseCell {base_cell} out of range 0..121")
    if not 1 <= resolution <= 3:
        raise ValueError(
            f"resolution {resolution} does not fit the 2-bit field (1..3); "
            f"the packed entry format must be widened first")
    if not 0 <= partial_index <= 511:
        raise ValueError(f"partialIndex {partial_index} out of range 0..511")
    if not 1 <= region_id <= 15:
        raise ValueError(
            f"regionId {region_id} does not fit the 4-bit field (1..15); "
            f"all 15 slots are allocated, so the packed entry format must "
            f"be widened before adding another region")
    return (base_cell << 25 | resolution << 23 |
            partial_index << 14 | region_id << 10)


def _cell_polygon(cell):
    """Shapely polygon for an H3 cell boundary (lng/lat order)."""
    boundary = h3.cell_to_boundary(cell)  # [(lat, lng), ...]
    from shapely.geometry import Polygon
    return Polygon([(lng, lat) for lat, lng in boundary])


def generate_lookup_table(regions, max_resolution=DEFAULT_RESOLUTION):
    """H3-2: uniform-only compaction with deterministic conflict resolution.

    1. Polyfill every region at res 3 only; build cell -> region.
    2. Conflicts (two regions claim the same res-3 cell): largest
       intersected polygon area wins; tie-break lowest region ID.
       Every resolution is logged.
    3. Bottom-up compaction: a group of res-3 siblings collapses to its
       res-2 parent only if ALL children are present and ALL map to the
       same region. Repeat res-2 -> res-1 under the same rule.
    4. Assert zero duplicate (baseCell, resolution, partialIndex) keys.

    Returns entries with keys: baseCell, resolution, partialIndex,
    regionId, regionName.
    """
    from shapely.ops import unary_union

    id_to_name = {v: k for k, v in REGION_IDS.items()}

    # ---- 1. polyfill every region at res 3 -------------------------------
    # h3 polyfill only keeps cells whose CENTER is inside the polygon, so
    # coastal cells erode (Sydney, Lagos). Fix: also polyfill a buffered
    # polygon, but keep only buffered cells that NO region claims truly —
    # i.e. the buffer extends coverage seaward only, never across a land
    # border into a neighbor.
    SEAWARD_BUFFER_DEG = 0.6  # ~ one res-3 cell radius (~60 km)

    true_cells = {}   # region_id -> set of cells from the true polygon
    unions = {}       # region_id -> unary union of true geometry
    for region_name, geoms in regions.items():
        if region_name not in REGION_IDS:
            print(f"Skipping region {region_name} - not in REGION_IDS")
            continue
        region_id = REGION_IDS[region_name]
        print(f"Polyfilling {region_name} (ID: {region_id}) at res 3")
        true_cells[region_id] = set(convert_region_to_h3(geoms,
                                                         max_resolution))
        unions[region_id] = unary_union(geoms)

    # Total planar extent per region, used as the specificity tie-break.
    region_extent = {rid: u.area for rid, u in unions.items()}

    all_true = set().union(*true_cells.values()) if true_cells else set()

    claims = {}  # cell -> list of region_ids claiming it
    for region_name, geoms in regions.items():
        if region_name not in REGION_IDS:
            continue
        region_id = REGION_IDS[region_name]
        buffered = [g.buffer(SEAWARD_BUFFER_DEG) for g in geoms]
        seaward = set(convert_region_to_h3(buffered, max_resolution)) \
            - all_true
        if seaward:
            print(f"  {region_name}: +{len(seaward)} seaward cells")
        for cell in true_cells[region_id] | seaward:
            claims.setdefault(cell, []).append(region_id)

    # ---- 2. deterministic conflict resolution ----------------------------
    cell_region = {}
    n_conflicts = 0
    for cell, region_ids in claims.items():
        if len(region_ids) == 1:
            cell_region[cell] = region_ids[0]
            continue
        n_conflicts += 1
        cpoly = _cell_polygon(cell)
        areas = []
        for rid in region_ids:
            area = unions[rid].intersection(cpoly).area
            # Tie-break on TOTAL region extent, smallest first: when two
            # plans both cover a cell completely (a national carve-out
            # nested inside a continental catch-all) the intersected areas
            # are equal and the more specific plan must win. The old
            # tie-break was "lowest region ID", which handed every interior
            # cell of the DRC to EU868 (ID 2) instead of CD900-1A (ID 15)
            # and left CD900-1A with 2 entries for a 2.34 Mkm2 country.
            areas.append((area, region_extent[rid], rid))
        areas.sort(key=lambda t: (-t[0], t[1], t[2]))
        winner = areas[0][2]
        cell_region[cell] = winner
        print(f"  conflict {cell}: " +
              ", ".join(f"{id_to_name[r]}={a:.6f}" for a, _, r in areas) +
              f" -> {id_to_name[winner]}")
    print(f"Resolved {n_conflicts} contested res-3 cells "
          f"(largest intersected area wins, smallest total region extent "
          f"breaks ties)")

    # ---- 3. bottom-up uniform compaction ---------------------------------
    for parent_res in (max_resolution - 1, max_resolution - 2):
        if parent_res < 1:
            break
        child_res = parent_res + 1
        # group cells currently at child_res by their parent
        by_parent = {}
        for cell, rid in cell_region.items():
            if h3.get_resolution(cell) != child_res:
                continue
            parent = h3.cell_to_parent(cell, parent_res)
            by_parent.setdefault(parent, {})[cell] = rid
        n_compacted = 0
        for parent, kids in by_parent.items():
            full = set(h3.cell_to_children(parent, child_res))
            if set(kids) != full:
                continue  # not all children present
            if len(set(kids.values())) != 1:
                continue  # not uniform
            rid = next(iter(kids.values()))
            for cell in kids:
                del cell_region[cell]
            cell_region[parent] = rid
            n_compacted += 1
        print(f"Compacted {n_compacted} uniform res-{child_res} groups "
              f"to res-{parent_res} parents")

    # ---- 4. emit + assert -------------------------------------------------
    entries = []
    seen = set()
    for cell, rid in cell_region.items():
        base_cell, partial_index, res = extract_h3_components(cell)
        key = (base_cell, res, partial_index)
        assert key not in seen, f"duplicate key {key} ({cell})"
        seen.add(key)
        entries.append({
            'h3': cell,
            'baseCell': base_cell,
            'resolution': res,
            'partialIndex': partial_index,
            'regionId': rid,
            'regionName': id_to_name[rid],
        })

    print(f"\nGenerated {len(entries)} lookup table entries "
          f"(0 duplicate keys asserted)")
    return entries


def _provenance(geojson_dir, resolution):
    """A stamp so a committed table can be traced back to its inputs."""
    import datetime
    import hashlib
    digests = []
    if geojson_dir and os.path.isdir(geojson_dir):
        for name in sorted(os.listdir(geojson_dir)):
            if not name.endswith('.geojson'):
                continue
            with open(os.path.join(geojson_dir, name), 'rb') as fh:
                digests.append(f" *   {name}  sha256:"
                               f"{hashlib.sha256(fh.read()).hexdigest()[:16]}")
    stamp = datetime.datetime.now(datetime.timezone.utc).isoformat(
        timespec='seconds')
    lines = [f" * Generated {stamp} by generate_lookup_table.py",
             f" * h3 python package: {getattr(h3, '__version__', 'unknown')}",
             f" * resolution: {resolution}",
             f" * source GeoJSON ({len(digests)} files):"]
    return "\n".join(lines + (digests or [" *   (none found)"]))


def generate_c_code(entries, output_c_file, output_h_file,
                    geojson_dir=None, resolution=DEFAULT_RESOLUTION):
    """Generate C code for the packed lookup table."""
    prov = _provenance(geojson_dir, resolution)
    # Numeric sort of the packed value == tuple sort by
    # (baseCell, resolution, partialIndex, regionId)
    entries.sort(key=lambda x: pack_entry(x['baseCell'], x['resolution'],
                                          x['partialIndex'], x['regionId']))

    max_id = max(REGION_IDS.values())

    # ---- header ----
    with open(output_h_file, 'w') as f:
        f.write("""/*
 * H3Lite Region Lookup Table
""" + prov + """
 * DO NOT EDIT BY HAND.
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
extern const char* const regionNames[{max_id + 1}];

// Binary search the region lookup table for an exact
// (baseCell, resolution, partialIndex) key; returns 0 (Unknown) on miss.
RegionId findRegion(uint8_t baseCell, uint8_t res, uint16_t partialIndex);

#endif /* H3LITE_REGIONS_TABLE_H */
""")

    # ---- C file ----
    with open(output_c_file, 'w') as f:
        f.write(f"""/*
 * H3Lite Region Lookup Table Implementation
{prov}
 * DO NOT EDIT BY HAND.
 */

#include "../include/h3lite_regions_table.h"

// Region names array
const char* const regionNames[{max_id + 1}] = {{
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

    # The packed entry gives resolution 2 bits, extract_h3_components()
    # keeps at most 3 digits, and the runtime probes res 3 -> 2 -> 1 only.
    # Any other --resolution would silently truncate (res 4 packs as 0).
    if args.resolution != 3:
        parser.error("The current packed H3Lite region table supports "
                     "resolution 3 only")

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
                    os.path.join('include', OUTPUT_H_FILE),
                    geojson_dir=args.geojson_dir,
                    resolution=args.resolution)

    print("\nLookup table generation complete!")


if __name__ == "__main__":
    main()
