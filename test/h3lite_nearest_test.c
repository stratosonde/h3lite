/*
 * H3Lite Nearest Regions Test
 * 
 * This program tests the findNearestRegions functionality by testing
 * points in known regions, at region boundaries, and over oceans.
 */

#include <stdio.h>
#include <stdlib.h>
#include "../include/h3lite.h"

typedef struct {
    double lat;
    double lng;
    const char* description;
} TestPoint;

int main(void) {
    printf("H3Lite Nearest Regions Test\n");
    printf("============================\n\n");
    
    // Initialize H3Lite
    if (!h3liteInit()) {
        printf("ERROR: Failed to initialize H3Lite\n");
        return 1;
    }
    
    // Test points
    TestPoint tests[] = {
        // Points inside regions
        {37.7749, -122.4194, "San Francisco (inside US915)"},
        {48.8566, 2.3522, "Paris (inside EU868)"},
        {-33.8688, 151.2093, "Sydney (inside AU915)"},
        
        // Points over ocean (should find nearest regions)
        {45.0, -30.0, "Mid-Atlantic (between US/EU)"},
        {0.0, -90.0, "Pacific Ocean (west of Americas)"},
        {0.0, 90.0, "Indian Ocean (east of Africa)"},
        
        // Edge cases
        {0.0, 0.0, "Gulf of Guinea"},
        {90.0, 0.0, "North Pole"},
        {-90.0, 0.0, "South Pole"},
    };
    
    int num_tests = sizeof(tests) / sizeof(TestPoint);
    
    printf("Testing %d locations with maxRings=3\n\n", num_tests);
    
    for (int i = 0; i < num_tests; i++) {
        printf("Location: %s\n", tests[i].description);
        printf("Coordinates: %.4f, %.4f\n", tests[i].lat, tests[i].lng);
        
        // Find nearest regions
        NearestRegionsInfo info = findNearestRegions(tests[i].lat, tests[i].lng, 3);
        
        if (info.numRegions == 0) {
            printf("Result: No regions found within 3 rings (~195km)\n");
        } else {
            printf("Found %d region(s):\n", info.numRegions);
            for (int j = 0; j < info.numRegions; j++) {
                printf("  %d. %s - Ring %d (~%.0f km)\n",
                       j + 1,
                       info.regions[j].regionName,
                       info.regions[j].ringDistance,
                       info.regions[j].distanceKm);
            }
        }
        printf("\n");
    }
    
    // Additional test: Compare with direct region lookup
    printf("=== Comparison Test ===\n");
    printf("Comparing findNearestRegions() with latLngToRegion() for known locations:\n\n");
    
    TestPoint comparison_tests[] = {
        {37.7749, -122.4194, "San Francisco"},
        {48.8566, 2.3522, "Paris"},
        {35.6762, 139.6503, "Tokyo"},
    };
    
    for (int i = 0; i < 3; i++) {
        RegionId directRegion = latLngToRegion(comparison_tests[i].lat, comparison_tests[i].lng);
        NearestRegionsInfo nearestInfo = findNearestRegions(comparison_tests[i].lat, comparison_tests[i].lng, 3);
        
        printf("%s:\n", comparison_tests[i].description);
        printf("  Direct lookup: %s\n", getRegionName(directRegion));
        if (nearestInfo.numRegions > 0) {
            printf("  Nearest regions: %s (ring %d)\n", 
                   nearestInfo.regions[0].regionName,
                   nearestInfo.regions[0].ringDistance);
            
            if (directRegion == nearestInfo.regions[0].regionId && 
                nearestInfo.regions[0].ringDistance == 0) {
                printf("  ✓ MATCH\n");
            } else {
                printf("  ✗ MISMATCH\n");
            }
        }
        printf("\n");
    }
    
    printf("Test complete!\n");
    return 0;
}
