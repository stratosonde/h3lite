# H3Lite

H3Lite is a minimalist implementation of Uber's H3 Hierarchical Geospatial Indexing System, specifically optimized for resource-constrained STM32 microcontrollers.

## Overview

This library provides essential H3 functionality for converting geographic coordinates (latitude/longitude) to H3 indexes and looking up geographic regions based on those indexes. The implementation is optimized for:

- Minimal flash memory usage (~5-10KB)
- Small RAM footprint (~1KB)
- Fixed-point math options for platforms without FPU
- Simplified algorithms for core use cases

H3Lite is particularly suited for embedded applications where the full H3 library would be too large, but you still need geospatial indexing capabilities.

## Use Case: LoRaWAN Region Detection

The primary use case for this implementation is a radiosonde device that needs to determine which LoRaWAN regional parameters to use based on its current geographic location. As the device travels, it converts its lat/lng coordinates to an H3 index and quickly looks up the appropriate radio region (US915, EU868, etc.) from a pre-computed table.

## Project Structure

```
h3lite/
├── include/               # Header files
│   ├── h3lite.h           # Main public API
│   ├── h3lite_constants.h # Constants and configuration
│   ├── h3lite_faceijk.h   # Coordinate conversion utilities
│   └── h3lite_regions_table.h # Generated region lookups
├── src/                   # Implementation files
│   ├── h3lite.c           # Core implementation
│   ├── h3lite_faceijk.c   # Coordinate conversion code
│   └── h3lite_regions_table.c # Generated lookup tables
├── test/                  # Test code (see archive/ for development tests)
├── generate_lookup_table.py # Python script to generate lookup tables
├── Makefile               # Build system
└── README.md              # This file
```

## How It Works

H3Lite includes three main components:

1. **Geographic Coordinate Conversion**: Converts lat/lng to an H3 index (simplified algorithm)
2. **Region Lookup Table**: Maps H3 indexes to geographic regions using a compact, pre-computed representation
3. **Region Identification API**: Simple functions to convert coordinates to regions with minimal processing

### Memory Optimization Approach

Rather than implementing the full H3 algorithm (which requires complex geometry calculations), H3Lite uses:

- Pre-computed lookup tables for region boundaries
- Simplified math for coordinate conversion
- Binary search for fast lookup
- Compact region representation

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

### Estimated Flash Memory Usage

- Core Code: ~2-4KB
- Lookup Tables: ~44KB (10,875 entries at resolution 3 for 12 LoRaWAN regions)
  - Each entry: 4 bytes (1 byte baseCell + 2 bytes partialIndex + 1 byte regionId)

### Estimated RAM Usage

- Static Data: ~100 bytes
- Runtime: ~50-100 bytes stack

## STM32 Integration

To integrate with your STM32 firmware:

1. Copy the include and src directories to your project
2. Include the appropriate headers in your code
3. Call `h3liteInit()` during your initialization
4. Use `latLngToRegion()` to determine the LoRaWAN region based on coordinates


## Limitations

- Only supports the core H3 functionality (lat/lng to H3 index)
- Doesn't implement the majority of the H3 API
- Optimized for low memory usage at the cost of some precision
- Resolution limited to 0-4 (configurable, with higher resolutions requiring more memory)

## License

Based on Uber's H3 library, licensed under the Apache License, Version 2.0.
