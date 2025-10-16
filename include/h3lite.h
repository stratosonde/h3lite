/*
 * H3Lite - Minimalist version of Uber's H3 Hexagonal Hierarchical Geospatial Indexing System
 * Optimized for STM32 microcontrollers with limited resources
 * 
 * Based on Uber's H3 library (https://github.com/uber/h3)
 * Licensed under the Apache License, Version 2.0
 */

#ifndef H3LITE_H
#define H3LITE_H

#include <stdint.h>
#include <stdbool.h>

#ifdef H3LITE_DEBUG
#define H3LITE_DEBUG_PRINT(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define H3LITE_DEBUG_PRINT(fmt, ...) ((void)0)
#endif

/**
 * H3 index is a 64-bit integer containing geospatial information
 */
typedef uint64_t H3Index;

/**
 * Simple latitude/longitude point in radians
 */
typedef struct {
    double lat;  // Latitude in radians
    double lng;  // Longitude in radians
} LatLng;

/**
 * Region identifier for LoRaWAN regional configuration
 */
typedef uint8_t RegionId;

/**
 * Maximum number of nearest regions to return
 */
#define MAX_NEAREST_REGIONS 3

/**
 * Information about a single region result
 */
typedef struct {
    RegionId regionId;        // Region ID (0 if slot unused)
    const char* regionName;   // Region name
    int ringDistance;         // Ring where found (0=current, 1-3)
    double distanceKm;        // Approximate distance in kilometers
} RegionResult;

/**
 * Information about multiple nearest regions
 */
typedef struct {
    int numRegions;                          // Number of regions found
    RegionResult regions[MAX_NEAREST_REGIONS]; // Array of found regions
} NearestRegionsInfo;

/**
 * Converts a lat/lng point to an H3 index at the specified resolution
 * 
 * @param lat Latitude in degrees
 * @param lng Longitude in degrees
 * @param resolution Resolution (2-4 recommended for memory constraints)
 * @return The H3 index or 0 on failure
 */
H3Index latLngToH3(double lat, double lng, int resolution);

/**
 * Looks up which region contains the given H3 index
 * 
 * @param h3 The H3 index to look up
 * @return The region ID or 0 if not found
 */
RegionId h3ToRegion(H3Index h3);

/**
 * Directly convert lat/lng to region ID in one step 
 * (internally uses latLngToH3 and h3ToRegion)
 * 
 * @param lat Latitude in degrees
 * @param lng Longitude in degrees
 * @return The region ID or 0 if not found
 */
RegionId latLngToRegion(double lat, double lng);

/**
 * Returns the string name of a region ID
 * 
 * @param regionId The region ID to look up
 * @return Pointer to static string name like "US915", "EU868", etc.
 */
const char* getRegionName(RegionId regionId);

/**
 * Find up to MAX_NEAREST_REGIONS nearest regions by searching surrounding H3 rings
 * Returns regions in order of discovery (regions in same ring have similar distance)
 * 
 * @param lat Latitude in degrees
 * @param lng Longitude in degrees  
 * @param maxRings Maximum rings to search (1-3 recommended)
 * @return NearestRegionsInfo with up to 3 nearest regions
 */
NearestRegionsInfo findNearestRegions(double lat, double lng, int maxRings);

/**
 * Initialize the H3Lite library
 * Must be called before using any other H3Lite functions
 * 
 * @return true if initialization succeeded, false otherwise
 */
bool h3liteInit(void);

#endif /* H3LITE_H */
