/*
 * H3Lite - Core implementation
 * Minimalist version of Uber's H3 Hexagonal Hierarchical Geospatial Indexing System
 * Optimized for STM32 microcontrollers with limited resources
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../include/h3lite.h"
#include "../include/h3lite_constants.h"
#include "../include/h3lite_faceijk.h"
#include "../include/h3lite_regions_table.h"

// Initialization flag
static bool initialized = false;

/**
 * Base cell lookup table 
 * This is a simplified version of the H3 base cell table with key representative values
 */
const int baseCellTable[12][12] = {
    // Face 0 (North pole)
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
    // Face 1 (North)
    {12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23},
    // Face 2 (North)
    {24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35},
    // Face 3 (North)
    {36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47},
    // Face 4 (North)
    {48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59},
    // Face 5 (North)
    {60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71},
    // Face 6 (South)
    {72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83},
    // Face 7 (South)
    {84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95},
    // Face 8 (South)
    {96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107},
    // Face 9 (South)
    {108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119},
    // Face 10 (South)
    {120, 121, 109, 108, 19, 78, 79, 80, 36, 35, 92, 37},
    // Face 11 (South pole)
    {104, 103, 102, 0, 1, 2, 3, 4, 5, 6, 7, 8}
};

// The region lookup table is now imported from h3lite_regions_table.h and 
// defined in h3lite_regions_table.c that was auto-generated

/**
 * Initialize the H3Lite library
 */
bool h3liteInit(void) {
    if (initialized) {
        return true;  // Already initialized
    }
    
    // Simple initialization for now
    // In a real implementation, this might load lookup data from flash
    initialized = true;
    return initialized;
}

/**
 * Convert lat/lng to H3 index
 * 
 * Simplified H3 implementation optimized for embedded systems.
 * Produces H3-compatible indexes for geolocation and region lookup.
 */
H3Index latLngToH3(double lat, double lng, int resolution) {
    // Constrain resolution (H3 supports 0-15, but we limit to lower resolutions)
    if (resolution < 0 || resolution > H3LITE_MAX_RESOLUTION) {
        return 0;  // Invalid resolution
    }
    
    // Convert lat/lng to radians
    double latRad = lat * M_PI_180;
    double lngRad = lng * M_PI_180;
    
    // Create FaceIJK structure
    FaceIJK fijk;
    
    // Convert geo coordinates to face and ijk coordinates
    LatLng g = {latRad, lngRad};
    _geoToFaceIjk(&g, resolution, &fijk);
    
    // Convert face/ijk coordinates to H3 index
    H3Index h3Index = faceIjkToH3(&fijk, resolution);
    
    return h3Index;
}

/**
 * Find region that contains the given H3 index
 */
RegionId h3ToRegion(H3Index h3) {
    if (!initialized || h3 == 0) {
        return INVALID_REGION;
    }
    
    // Extract the resolution
    int res = H3_GET_RESOLUTION(h3);
    
    // Extract the base cell
    int baseCell = H3_GET_BASE_CELL(h3);
    
    // Extract partial indexing information (using 3 digits)
    uint32_t partialIndex = 0;
    int numDigits = (res < 3) ? res : 3;  // Use up to 3 digits
    for (int r = 1; r <= numDigits; r++) {
        int digit = H3_GET_INDEX_DIGIT(h3, r);
        partialIndex = (partialIndex * 8) + digit;  // Base-8 encoding
    }
    
    // Debug output
    printf("DEBUG h3ToRegion: h3=0x%llx baseCell=%d partialIndex=%d\n", 
           (unsigned long long)h3, baseCell, partialIndex);
    
    // Search for matching region
    RegionId region = findRegion(baseCell, partialIndex);
    printf("DEBUG h3ToRegion: Search result = %d (%s)\n", 
           region, getRegionName(region));
    
    return region;
}

/**
 * Directly convert lat/lng to region ID
 */
RegionId latLngToRegion(double lat, double lng) {
    // Use the target resolution (hardcoded to 3 for now)
    H3Index h3 = latLngToH3(lat, lng, H3LITE_TARGET_RESOLUTION);
    if (h3 == 0) {
        return INVALID_REGION;
    }
    
    return h3ToRegion(h3);
}

/**
 * Get string name for a region ID
 */
const char* getRegionName(RegionId regionId) {
    // Use the region names from the generated lookup table
    // Make sure we don't go out of bounds of the array size
    int max_region_id = 12; // From the generated table
    if (regionId > max_region_id) {
        return regionNames[0];  // "Unknown"
    }
    return regionNames[regionId];
}

/**
 * Find up to MAX_NEAREST_REGIONS nearest regions by searching surrounding H3 rings
 * Returns regions in order of discovery (regions in same ring have similar distance)
 */
NearestRegionsInfo findNearestRegions(double lat, double lng, int maxRings) {
    NearestRegionsInfo result = {0};
    H3Index h3 = latLngToH3(lat, lng, H3LITE_TARGET_RESOLUTION);
    
    // Try current cell first
    RegionId region = h3ToRegion(h3);
    if (region != 0) {
        result.numRegions = 1;
        result.regions[0].regionId = region;
        result.regions[0].regionName = getRegionName(region);
        result.regions[0].ringDistance = 0;
        result.regions[0].distanceKm = 0.0;
        return result;
    }
    
    // PRODUCTION-SAFE: Limit maxRings to prevent buffer overflow
    #define MAX_SUPPORTED_RINGS 6
    #define RING_BUFFER_SIZE 42  // 6*6=36, +6 safety margin
    
    if (maxRings > MAX_SUPPORTED_RINGS) {
        maxRings = MAX_SUPPORTED_RINGS;  // Clamp to safe maximum
    }
    
    // Search rings 1 to maxRings
    for (int k = 1; k <= maxRings && result.numRegions < MAX_NEAREST_REGIONS; k++) {
        int ringSize = 6 * k;
        
        // BOUNDS CHECK: Verify ring size won't overflow buffer
        if (ringSize > RING_BUFFER_SIZE) {
            break;  // Stop if ring would overflow - prevents corruption
        }
        
        H3Index ringCells[RING_BUFFER_SIZE];
        
        // Get the ring of cells at distance k
        if (h3GetRing(h3, k, ringCells) == 0) {
            // Check all cells in this ring
            for (int i = 0; i < ringSize && result.numRegions < MAX_NEAREST_REGIONS; i++) {
                region = h3ToRegion(ringCells[i]);
                if (region != 0) {
                    // Check if we already have this region
                    bool alreadyFound = false;
                    for (int j = 0; j < result.numRegions; j++) {
                        if (result.regions[j].regionId == region) {
                            alreadyFound = true;
                            break;
                        }
                    }
                    
                    if (!alreadyFound) {
                        // Add new region
                        result.regions[result.numRegions].regionId = region;
                        result.regions[result.numRegions].regionName = getRegionName(region);
                        result.regions[result.numRegions].ringDistance = k;
                        // Approximate distance: ~65km per ring at resolution 3
                        result.regions[result.numRegions].distanceKm = k * 65.0;
                        result.numRegions++;
                    }
                }
            }
        }
        
        // If we found at least one region in this ring, stop searching
        // (all regions in same ring are roughly equidistant)
        if (result.numRegions > 0) break;
    }
    
    return result;
}

// Removed unused helper function
