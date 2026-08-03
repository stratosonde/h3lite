# H3Lite

A minimalist implementation of Uber's H3 Hierarchical Geospatial Indexing System, specifically optimized for resource-constrained embedded systems like STM32 microcontrollers.

## Overview

H3Lite provides essential H3 functionality for converting geographic coordinates (latitude/longitude) to H3 indexes and performing fast region lookups. The implementation is optimized for:

- **Minimal flash memory usage** (~5-10KB for core code)
- **Small RAM footprint** (~1KB)
- **Fixed-point math** options for platforms without FPU
- **Simplified algorithms** focused on core use cases
- **Pre-computed lookup tables** for efficient region detection

H3Lite is ideal for embedded applications where the full H3 library would be too large, but you still need reliable geospatial indexing capabilities.

## Use Case: LoRaWAN Region Detection

The primary use case is automatic LoRaWAN regional parameter detection for mobile devices (e.g., radiosondes, asset trackers). As a device moves globally, it converts its GPS coordinates to an H3 index and quickly looks up the appropriate radio region (US915, EU868, AS923-1, etc.) from a pre-computed lookup table, enabling automatic frequency plan selection without manual configuration.

## Features

- **Geographic to H3 Conversion**: Convert lat/lng coordinates to H3 indexes at resolutions 0-4 (bit-exact with reference H3, validated over a 0.5° global grid)
- **Region Lookup**: Fast binary search through pre-computed region tables
- **LoRaWAN Region Support**: 14 LoRaWAN region plans plus a RESTRICTED marker (ID 15, transmission prohibited by policy)
- **Embeddable**: Designed for STM32 and other ARM Cortex-M microcontrollers
- **Portable**: Standard C99 with minimal dependencies

## Project Structure

```
h3lite_github/
├── include/                     # Header files
│   ├── h3lite.h                 # Main public API
│   ├── h3lite_constants.h       # Constants and configuration
│   ├── h3lite_faceijk.h         # Coordinate conversion utilities
│   └── h3lite_regions_table.h   # Generated region lookup tables
├── src/                         # Implementation files
│   ├── h3lite.c                 # Core H3 implementation
│   ├── h3lite_faceijk.c         # Coordinate conversion
│   ├── h3lite_basecells.c       # Base-cell / icosa-face tables
│   ├── h3lite_neighbor.c        # k-ring traversal (h3GetRing)
│   └── h3lite_regions_table.c   # Generated region tables
├── test/                        # Test programs
│   ├── h3lite_grid_test.c       # Grid coverage tests
│   ├── h3lite_nearest_test.c    # Nearest neighbor tests
│   ├── visual_grid_test.py      # Visualization tools
│   ├── visualize_table.py       # Table inspection utilities
│   └── simple_map.py            # Simple mapping utilities
├── obj/                         # Build artifacts (generated)
├── archive/                     # Development history & debug tools
├── generate_lookup_table.py     # Table generation script
├── regions_h3_conversion.py     # Region processing utilities
├── analyze_regions.py           # Region analysis tools
├── Makefile                     # Build system
├── LICENSE                      # Apache 2.0 License
└── README.md                    # This file
```

## How It Works

### Core Components

1. **Geographic Coordinate Conversion**: Simplified algorithm to convert lat/lng to H3 indexes
2. **Pre-computed Region Tables**: Compact lookup tables mapping H3 cells to regions
3. **Binary Search**: Fast O(log n) region lookup
4. **Region API**: Simple functions for coordinate-to-region conversion

### Memory Optimization Strategy

Instead of implementing the full H3 specification, H3Lite uses:

- **Pre-computed lookup tables** for region boundaries at resolution 3
- **Simplified coordinate math** optimized for embedded systems
- **Binary search** for O(log n) lookups instead of spatial queries
- **Compact encoding**: 4 bytes per table entry (packed uint32: baseCell, resolution, partialIndex, regionId)
- **Regional focus**: Only includes cells needed for region detection

## Building and Using

### Standard Build

```bash
# Build the static library
make lib

# Build and run the test programs (non-zero exit on failure)
make test
./bin/h3lite_grid_test
./bin/h3lite_nearest_test
```

### STM32 Cross-Compilation

```bash
# Build for STM32
make stm32lib
```

### Generating Lookup Tables

To generate the compact region lookup tables:

```bash
# Ensure you have the required Python packages
pip install h3 shapely tqdm

# Generate tables from GeoJSON region definitions
python generate_lookup_table.py
```

The script will process GeoJSON region files in the working directory and generate `src/h3lite_regions_table.c` and `include/h3lite_regions_table.h` files.

## API Usage

```c
#include "h3lite.h"

// Initialize the library
h3liteInit();

// Convert latitude/longitude to H3 index
double lat = 37.775;
double lng = -122.418;
H3Index h3 = latLngToH3(lat, lng, 3);  // Resolution 3

// Look up which region contains this H3 index
RegionId region = h3ToRegion(h3);

// Direct lat/lng to region lookup
RegionId region = latLngToRegion(lat, lng);

// Get region name
const char* regionName = getRegionName(region);  // e.g., "US915"
```

## Memory Usage

### Flash Memory (STM32)

- **Core Code**: ~2-4KB
- **Geometry / base-cell tables**: ~4KB (int8/uint8-narrowed const tables)
- **Region Lookup Table**: 18,256 bytes (4,564 packed entries, mixed res 1-3)
  - Each entry: 4 bytes packed uint32 (baseCell[31:25] resolution[24:23]
    partialIndex[22:14] regionId[13:10])
- **Total**: ~24-26KB

### RAM Usage

- **Static Data**: ~100 bytes (all large tables are `const` and live in flash;
  `regionNames` is `const char* const` so the pointer array is flash-resident too)
- **Stack**: ~50-100 bytes per operation; `findNearestRegions` uses a
  42-entry (336 byte) stack buffer worst case
- **Total**: <1KB

### Supported Regions

14 LoRaWAN region plans:
- EU868, US915, CN470, AU915
- AS923-1, AS923-1B, AS923-1C, AS923-2, AS923-3, AS923-4
- KR920, IN865, RU864, EU433
- RESTRICTED (ID 15 — no transmission allowed; enforced by the application)
- Unknown (ID 0, default/fallback)

The regionId field in the packed entry is 4 bits (IDs 0-15). ID 15 was
repurposed from the CD900-1A test plan, which is not used in production.
`REGION_RESTRICTED` is therefore 15 — a value of 255 could never be
emitted by the table (it would silently truncate to 15).

## STM32 Integration

### Step 1: Add to Project

Copy the following to your STM32 project:
- `include/` directory → Your project's include path
- ALL five source files → Your source files:
  `h3lite.c`, `h3lite_faceijk.c`, `h3lite_basecells.c`,
  `h3lite_neighbor.c`, `h3lite_regions_table.c`
  (omitting `h3lite_basecells.c`/`h3lite_neighbor.c` builds fine but
  link-fails the moment `findNearestRegions()`/`h3GetRing()` is called)

### Step 2: Configure Build

Add to your Makefile or IDE:
```makefile
CFLAGS += -Ipath/to/h3lite/include
SOURCES += h3lite.c h3lite_faceijk.c h3lite_basecells.c \
           h3lite_neighbor.c h3lite_regions_table.c
```

### Step 3: Use in Code

```c
#include "h3lite.h"

// In your initialization
h3liteInit();

// In your GPS handler
void onGPSUpdate(double lat, double lng) {
    RegionId region = latLngToRegion(lat, lng);
    const char* name = getRegionName(region);
    
    // Configure LoRaWAN radio for this region
    lorawan_set_region(region);
}
```

## Development

### Building Tests

```bash
make test          # Build and run basic tests
make grid_test     # Build grid coverage test
make nearest_test  # Build nearest neighbor test
```

### Generating New Lookup Tables

If you need to regenerate the region lookup tables (e.g., for different regions or resolution):

```bash
# Install dependencies
pip install h3 shapely tqdm

# Generate tables from GeoJSON files
python generate_lookup_table.py

# This will regenerate:
#   - src/h3lite_regions_table.c
#   - include/h3lite_regions_table.h
```

### Analyzing Regions

```bash
python analyze_regions.py      # Analyze region coverage
python regions_h3_conversion.py # Process GeoJSON regions
```

## Limitations

- **Simplified H3 implementation**: Only core functionality (lat/lng → H3 index)
- **Limited API surface**: Doesn't implement the full H3 specification
- **Resolution constraint**: Optimized for resolutions 0-4 (higher resolutions possible but require more memory)
- **Precision trade-off**: Optimized for memory efficiency over maximum precision
- **Region-focused**: Designed specifically for region lookup, not general H3 operations

## Contributing

Contributions are welcome! Areas for improvement:
- Additional region definitions
- Memory optimizations
- Extended H3 API support
- Platform-specific optimizations

## License

Based on Uber's H3 library, licensed under the Apache License, Version 2.0.

## Acknowledgments

- **Uber H3**: Original H3 geospatial indexing system
- **LoRaWAN Alliance**: Regional parameter definitions
- Built with support from the embedded systems community

## Known limitations

- **Pentagon rings — corrected statement.** `h3GetRing()` fails cleanly
  (returns negative) ONLY for `k=1` rings whose origin is itself one of
  the 12 pentagon cells (H3's legitimate pentagon distortion case). Every
  other origin — including hexagons adjacent to pentagons — produces
  rings bit-exact with reference H3.
  The earlier claim that failures "cluster near pentagons" and that
  `findNearestRegions()` "searches the next ring outward" was WRONG on
  both counts: a predicate bug (`_h3IsPentagon` tested only the base
  cell, not the all-zero-digits condition) made ~6% of the globe report
  ring failure, and a ring-1 failure aborted the whole search instead of
  continuing. Both are fixed; regression coverage lives in
  `test/t_pentagon.py` (differential vs reference h3 around all 12
  pentagons). The 12 pentagon-centred regions include the Bohai
  Sea/China, Gulf of Alaska, Norwegian Sea, Arabian Sea, Gulf of Guinea,
  and waters off Argentina and NW Australia — so this mattered for
  worldwide operation.

- **Region border granularity.** The table polyfills cell CENTRES at
  res 3 (~69 km edge, ~120 km between centres) plus a ~0.6° seaward
  buffer. Border placement is therefore good to roughly one cell —
  tens of km, not metres. Points within ~one cell of a land border may
  resolve to the neighbour plan.

- **Coverage gaps in the source data (hplans), as of 2026-08.**
  Probed and verified against the hplans polygons; all of these fall
  through to `Unknown` (0) because NO region polygon contains them —
  this is a data gap, not a library defect:
  Mongolia, Fiji, Tuvalu, Kiribati, Marshall Islands, Micronesia (FSM),
  Nauru, Palau, Maldives. Mitigation: `findNearestRegions()` picks up
  the closest plan offshore of these countries; closing the gap
  requires extending the hplans source polygons (tracked separately —
  see the firmware repo's h3lite verification report).

- **First-productive-ring semantics.** `findNearestRegions()` returns
  regions from the FIRST ring containing any known region (up to 3),
  not the globally nearest three. `distanceKm` is approximate at
  ~120 km per ring at res 3 (previously mis-stated as ~65 km/ring).


