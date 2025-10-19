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

- **Geographic to H3 Conversion**: Convert lat/lng coordinates to H3 indexes at resolutions 0-4
- **Region Lookup**: Fast binary search through pre-computed region tables
- **LoRaWAN Region Support**: Built-in support for 16 LoRaWAN regions worldwide
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
- **Compact encoding**: 4 bytes per table entry (baseCell + partialIndex + regionId)
- **Regional focus**: Only includes cells needed for region detection

## Building and Using

### Standard Build

```bash
# Build the static library
make lib

# Build and run the test program
make test
./bin/h3lite_test
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
pip install h3 geopandas shapely numpy tqdm

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
- **Lookup Tables**: ~44KB (10,875 entries at resolution 3)
  - Each entry: 4 bytes (1 byte baseCell + 2 bytes partialIndex + 1 byte regionId)
- **Total**: ~46-48KB

### RAM Usage

- **Static Data**: ~100 bytes
- **Stack**: ~50-100 bytes per operation
- **Total**: <1KB

### Supported Regions

Currently supports 16 LoRaWAN regions:
- EU868, US915, CN470, AU915
- AS923-1, AS923-1B, AS923-1C, AS923-2, AS923-3, AS923-4
- KR920, IN865, RU864, EU433, CD900-1A
- Unknown (default/fallback)

## STM32 Integration

### Step 1: Add to Project

Copy the following to your STM32 project:
- `include/` directory → Your project's include path
- `src/h3lite.c` and `src/h3lite_faceijk.c` → Your source files
- `src/h3lite_regions_table.c` → Your source files

### Step 2: Configure Build

Add to your Makefile or IDE:
```makefile
CFLAGS += -Ipath/to/h3lite/include
SOURCES += h3lite.c h3lite_faceijk.c h3lite_regions_table.c
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
pip install h3 geopandas shapely numpy tqdm matplotlib

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
