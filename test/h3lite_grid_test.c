/*
 * H3Lite Grid Coverage Test
 * 
 * This program systematically tests H3Lite by converting a coarse grid
 * of lat/lng coordinates to H3 indexes and verifying region lookups.
 * This helps validate coverage across the entire globe and identify
 * any gaps in the region lookup table.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../include/h3lite.h"
#include "../include/h3lite_constants.h"

// Test configuration
#define GRID_STEP_DEGREES 10  // Test every 10 degrees
#define TEST_RESOLUTION 3     // H3 resolution to test

// Statistics structure
typedef struct {
    int total_tests;
    int valid_h3;
    int invalid_h3;
    int region_found;
    int region_unknown;
    int region_counts[16];  // Count for each region ID
} TestStats;

// Function prototypes
void test_global_grid(TestStats* stats);
void test_known_locations(TestStats* stats);
void test_boundary_cases(TestStats* stats);
void print_statistics(const TestStats* stats);
void print_h3_details(H3Index h3, double lat, double lng);

int main(void) {
    printf("H3Lite Grid Coverage Test\n");
    printf("=========================\n");
    printf("Grid step: %d degrees\n", GRID_STEP_DEGREES);
    printf("Test resolution: %d\n", TEST_RESOLUTION);
    printf("\n");
    
    // Initialize H3Lite
    if (!h3liteInit()) {
        printf("ERROR: Failed to initialize H3Lite\n");
        return 1;
    }
    printf("H3Lite initialized successfully\n\n");
    
    // Initialize statistics
    TestStats stats = {0};
    
    // Start timing
    clock_t start_time = clock();
    
    // Run test suites
    printf("=== Running Global Grid Coverage Test ===\n");
    test_global_grid(&stats);
    
    printf("\n=== Running Known Locations Test ===\n");
    test_known_locations(&stats);
    
    printf("\n=== Running Boundary Cases Test ===\n");
    test_boundary_cases(&stats);
    
    // End timing
    clock_t end_time = clock();
    double elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    
    // Print final statistics
    printf("\n=== Test Summary ===\n");
    print_statistics(&stats);
    printf("\nTotal execution time: %.3f seconds\n", elapsed);
    printf("Average time per conversion: %.3f ms\n", 
           (elapsed * 1000.0) / stats.total_tests);
    
    return 0;
}

void test_global_grid(TestStats* stats) {
    int points_tested = 0;
    
    printf("Testing lat/lng grid from -90° to +90° latitude, -180° to +180° longitude\n");
    printf("Step size: %d degrees\n\n", GRID_STEP_DEGREES);
    
    // Output format header
    printf("%-6s %-6s %-20s %-12s %-10s %-15s\n", 
           "Lat", "Lng", "H3 Index", "BaseCell", "PartIdx", "Region");
    printf("------------------------------------------------------------------------\n");
    
    for (int lat = -90; lat <= 90; lat += GRID_STEP_DEGREES) {
        for (int lng = -180; lng <= 180; lng += GRID_STEP_DEGREES) {
            stats->total_tests++;
            
            // Convert to H3
            H3Index h3 = latLngToH3((double)lat, (double)lng, TEST_RESOLUTION);
            
            if (h3 == 0) {
                stats->invalid_h3++;
                printf("%-6d %-6d %-20s %-12s %-10s %-15s\n",
                       lat, lng, "INVALID", "-", "-", "-");
                continue;
            }
            
            stats->valid_h3++;
            
            // Extract H3 components
            int baseCell = H3_GET_BASE_CELL(h3);
            
            // Calculate partial index
            uint16_t partialIndex = 0;
            for (int r = 1; r <= TEST_RESOLUTION; r++) {
                int digit = H3_GET_INDEX_DIGIT(h3, r);
                partialIndex = (partialIndex * 8) + digit;
            }
            
            // Look up region
            RegionId regionId = h3ToRegion(h3);
            const char* regionName = getRegionName(regionId);
            
            // Update statistics
            if (regionId > 0 && regionId < 16) {
                stats->region_counts[regionId]++;
                stats->region_found++;
            } else {
                stats->region_unknown++;
            }
            
            // Print result
            printf("%-6d %-6d 0x%016llx %-12d %-10d %-15s\n",
                   lat, lng, (unsigned long long)h3, baseCell, partialIndex, regionName);
            
            points_tested++;
        }
    }
    
    printf("\nGrid test completed: %d points tested\n", points_tested);
}

void test_known_locations(TestStats* stats) {
    // Known test locations
    typedef struct {
        double lat;
        double lng;
        const char* name;
        const char* expected_region;
    } TestLocation;
    
    TestLocation locations[] = {
        { 37.7749, -122.4194, "San Francisco, USA", "US915" },
        { 40.7128, -74.0060, "New York, USA", "US915" },
        { 34.0522, -118.2437, "Los Angeles, USA", "US915" },
        { 48.8566, 2.3522, "Paris, France", "EU868" },
        { 51.5074, -0.1278, "London, UK", "EU868" },
        { 52.5200, 13.4050, "Berlin, Germany", "EU868" },
        { -33.8688, 151.2093, "Sydney, Australia", "AU915" },
        { -37.8136, 144.9631, "Melbourne, Australia", "AU915" },
        { 35.6762, 139.6503, "Tokyo, Japan", "AS923-1" },
        { 37.5665, 126.9780, "Seoul, South Korea", "KR920" },
        { 28.6139, 77.2090, "New Delhi, India", "IN865" },
        { 55.7558, 37.6173, "Moscow, Russia", "RU864" },
        { 39.9042, 116.4074, "Beijing, China", "CN470" },
    };
    
    int num_locations = sizeof(locations) / sizeof(TestLocation);
    
    printf("Testing %d known locations\n\n", num_locations);
    printf("%-25s %-15s %-15s %-10s\n", "Location", "Region", "Expected", "Result");
    printf("--------------------------------------------------------------------\n");
    
    for (int i = 0; i < num_locations; i++) {
        stats->total_tests++;
        
        H3Index h3 = latLngToH3(locations[i].lat, locations[i].lng, TEST_RESOLUTION);
        
        if (h3 == 0) {
            stats->invalid_h3++;
            printf("%-25s %-15s %-15s %-10s\n",
                   locations[i].name, "INVALID", locations[i].expected_region, "FAIL");
            continue;
        }
        
        stats->valid_h3++;
        
        RegionId regionId = h3ToRegion(h3);
        const char* regionName = getRegionName(regionId);
        
        if (regionId > 0 && regionId < 16) {
            stats->region_counts[regionId]++;
            stats->region_found++;
        } else {
            stats->region_unknown++;
        }
        
        // Check if result matches expected
        const char* result = (strcmp(regionName, locations[i].expected_region) == 0) ? "PASS" : "FAIL";
        
        printf("%-25s %-15s %-15s %-10s\n",
               locations[i].name, regionName, locations[i].expected_region, result);
    }
}

void test_boundary_cases(TestStats* stats) {
    typedef struct {
        double lat;
        double lng;
        const char* description;
    } BoundaryTest;
    
    BoundaryTest tests[] = {
        { 0.0, 0.0, "Equator/Prime Meridian" },
        { 0.0, 180.0, "Equator/Date Line (East)" },
        { 0.0, -180.0, "Equator/Date Line (West)" },
        { 90.0, 0.0, "North Pole" },
        { -90.0, 0.0, "South Pole" },
        { 89.0, 0.0, "Near North Pole" },
        { -89.0, 0.0, "Near South Pole" },
        { 0.0, 179.9, "Near Date Line (East)" },
        { 0.0, -179.9, "Near Date Line (West)" },
        { 23.5, 0.0, "Tropic of Cancer" },
        { -23.5, 0.0, "Tropic of Capricorn" },
        { 66.5, 0.0, "Arctic Circle" },
        { -66.5, 0.0, "Antarctic Circle" },
    };
    
    int num_tests = sizeof(tests) / sizeof(BoundaryTest);
    
    printf("Testing %d boundary cases\n\n", num_tests);
    printf("%-30s %-20s %-15s\n", "Description", "H3 Index", "Region");
    printf("--------------------------------------------------------------------\n");
    
    for (int i = 0; i < num_tests; i++) {
        stats->total_tests++;
        
        H3Index h3 = latLngToH3(tests[i].lat, tests[i].lng, TEST_RESOLUTION);
        
        if (h3 == 0) {
            stats->invalid_h3++;
            printf("%-30s %-20s %-15s\n",
                   tests[i].description, "INVALID", "-");
            continue;
        }
        
        stats->valid_h3++;
        
        RegionId regionId = h3ToRegion(h3);
        const char* regionName = getRegionName(regionId);
        
        if (regionId > 0 && regionId < 16) {
            stats->region_counts[regionId]++;
            stats->region_found++;
        } else {
            stats->region_unknown++;
        }
        
        printf("%-30s 0x%016llx %-15s\n",
               tests[i].description, (unsigned long long)h3, regionName);
    }
}

void print_statistics(const TestStats* stats) {
    printf("Total tests run: %d\n", stats->total_tests);
    printf("Valid H3 indexes: %d (%.1f%%)\n", 
           stats->valid_h3, 
           100.0 * stats->valid_h3 / stats->total_tests);
    printf("Invalid H3 indexes: %d (%.1f%%)\n", 
           stats->invalid_h3,
           100.0 * stats->invalid_h3 / stats->total_tests);
    printf("\n");
    
    printf("Region lookup results:\n");
    printf("  Found region: %d (%.1f%%)\n", 
           stats->region_found,
           100.0 * stats->region_found / stats->valid_h3);
    printf("  Unknown region: %d (%.1f%%)\n", 
           stats->region_unknown,
           100.0 * stats->region_unknown / stats->valid_h3);
    printf("\n");
    
    printf("Regions detected:\n");
    const char* region_names[] = {
        "Unknown", "US915", "EU868", "AU915", "AS923-1", "AS923-2",
        "AS923-3", "AS923-4", "KR920", "IN865", "RU864", "CN470", "EU433"
    };
    
    for (int i = 0; i <= 12; i++) {
        if (stats->region_counts[i] > 0) {
            printf("  %s: %d points (%.1f%%)\n", 
                   region_names[i],
                   stats->region_counts[i],
                   100.0 * stats->region_counts[i] / stats->valid_h3);
        }
    }
}
