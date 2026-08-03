/*
 * H3Lite Grid Coverage Test
 * 
 * This program systematically tests H3Lite by converting a coarse grid
 * of lat/lng coordinates to H3 indexes and verifying region lookups.
 * This helps validate coverage across the entire globe and identify
 * any gaps in the region lookup table.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "h3lite.h"
#include "h3lite_constants.h"

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
    int known_failures;    // known-location / invalid-input mismatches
    int region_counts[16];  // Count for each region ID
} TestStats;

// Function prototypes
void test_global_grid(TestStats* stats);
void test_known_locations(TestStats* stats);
void test_boundary_cases(TestStats* stats);
void test_invalid_inputs(TestStats* stats);
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

    printf("\n=== Running Invalid Input Test ===\n");
    test_invalid_inputs(&stats);

    // End timing
    clock_t end_time = clock();
    double elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // Print final statistics
    printf("\n=== Test Summary ===\n");
    print_statistics(&stats);
    printf("\nTotal execution time: %.3f seconds\n", elapsed);
    printf("Average time per conversion: %.3f ms\n",
           (elapsed * 1000.0) / stats.total_tests);

    /* A test that cannot fail is not a test: report known-location and
     * invalid-input mismatches via the exit code (review F-18). */
    if (stats.known_failures > 0) {
        printf("\nRESULT: FAIL (%d known/invalid-input mismatches)\n",
               stats.known_failures);
        return 1;
    }
    printf("\nRESULT: PASS\n");
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
        int pass = (strcmp(regionName, locations[i].expected_region) == 0);
        if (!pass) {
            stats->known_failures++;
        }
        const char* result = pass ? "PASS" : "FAIL";

        printf("%-25s %-15s %-15s %-10s\n",
               locations[i].name, regionName, locations[i].expected_region, result);
    }
}

/* NaN / inf / out-of-domain coordinates must produce h3 == 0, never a
 * valid-looking index. On Cortex-M the NaN->int cast flushes to 0, which
 * previously walked straight through to a real cell and could select a
 * real LoRaWAN region from a lost GPS fix (review F-03). */
void test_invalid_inputs(TestStats* stats) {
    typedef struct {
        double lat;
        double lng;
        const char* description;
    } InvalidTest;

    InvalidTest tests[] = {
        { NAN, 0.0, "NaN latitude" },
        { 0.0, NAN, "NaN longitude" },
        { INFINITY, 0.0, "+inf latitude" },
        { 0.0, -INFINITY, "-inf longitude" },
        { 91.0, 0.0, "latitude > 90" },
        { -90.1, 0.0, "latitude < -90" },
        { 0.0, 180.1, "longitude > 180" },
        { 0.0, -181.0, "longitude < -180" },
    };

    int num_tests = sizeof(tests) / sizeof(InvalidTest);

    printf("Testing %d invalid inputs (all must yield h3 == 0)\n\n", num_tests);
    printf("%-25s %-20s %-10s\n", "Input", "H3 Index", "Result");
    printf("----------------------------------------------------------\n");

    for (int i = 0; i < num_tests; i++) {
        stats->total_tests++;
        H3Index h3 = latLngToH3(tests[i].lat, tests[i].lng, TEST_RESOLUTION);
        int pass = (h3 == 0);
        if (!pass) {
            stats->known_failures++;
        }
        printf("%-25s 0x%016llx %-10s\n",
               tests[i].description, (unsigned long long)h3,
               pass ? "PASS" : "FAIL");
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
    if (stats->total_tests == 0 || stats->valid_h3 == 0) {
        printf("No results to summarise.\n");
        return;
    }
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
    /* Use the library's own name table. A private copy here previously
     * skipped AS923-1B/1C, so every ID >= 5 was mislabelled, and the loop
     * stopped at 12 so CN470/EU433/CD900-1A never appeared at all. */
    for (int i = 1; i < (int)(sizeof(stats->region_counts) /
                              sizeof(stats->region_counts[0])); i++) {
        if (stats->region_counts[i] > 0) {
            printf("  %s: %d points (%.1f%%)\n",
                   getRegionName((RegionId)i),
                   stats->region_counts[i],
                   100.0 * stats->region_counts[i] / stats->valid_h3);
        }
    }
}
