/*
 * H3Lite - Core implementation
 * Minimalist version of Uber's H3 Hexagonal Hierarchical Geospatial Indexing System
 * Optimized for STM32 microcontrollers with limited resources
 */

#include "h3lite.h"
#include "h3lite_constants.h"
#include "h3lite_faceijk.h"
#include "h3lite_regions_table.h"
#include <math.h>
#include <string.h>

// Initialization flag
static bool initialized = false;

/* H3-9: baseCellTable[12][12] deleted — it was a simplified/fabricated table
 * that no code path actually consumed (only an unused extern in faceijk).
 * 576 B of misleading dead data removed. */

// The region lookup table is now imported from h3lite_regions_table.h and
// defined in h3lite_regions_table.c that was auto-generated

/**
 * Initialize the H3Lite library
 */
bool h3liteInit(void) {
  if (initialized) {
    return true; // Already initialized
  }

  /* F-013 (#210): this used to be an impossible-failure API - it set the
   * flag and returned true unconditionally, making the caller's
   * Error_Handler branch dead code. Now it performs a real integrity
   * self-check of the geo pipeline AND the region table before
   * reporting success:
   *  1. the coordinate->H3 pipeline produces a nonzero index for a known
   *     land coordinate;
   *  2. the region table maps that coordinate (Paris) to EU868 - catches
   *     a truncated/corrupt/regenerated-wrong table;
   *  3. a mid-Atlantic coordinate resolves to REGION_UNKNOWN - a positive
   *     probe for the no-region path. This does NOT prove global land
   *     coverage: one ocean probe cannot rule out unmapped land elsewhere.
   *     The safety contract is the F-006 (#208) / GEO-02 disposition:
   *     selected territories are marked REGION_RESTRICTED in the dataset
   *     (#257, 53 cells) and silenced; every other UNKNOWN cell (ocean
   *     AND any land not so marked) transmits on the held region by
   *     product decision, not by dataset claim.
   * Failure here means the geofence data cannot be trusted: the caller's
   * failure branch is now meaningful. */
  /* h3ToRegion refuses lookups while !initialized, so arm the flag first
   * and drop it again if any probe fails. */
  initialized = true;
  H3Index probe = latLngToH3(48.8566, 2.3522, H3LITE_TARGET_RESOLUTION);
  if (probe == 0) {
    initialized = false;
    return false;
  }
  if (strcmp(getRegionName(latLngToRegion(48.8566, 2.3522)), "EU868") != 0) {
    initialized = false;
    return false;
  }
  if (latLngToRegion(0.0, -30.0) != REGION_UNKNOWN) {
    initialized = false;
    return false;
  }

  /* GEO-04: the F-013 probes prove the table is present and well-formed,
   * but silently losing the restricted-enforcement dataset (GEO-01: Yemen
   * + North Korea) would pass all three - Paris and the mid-Atlantic are
   * unaffected by dropping those cells. Scan the table for at least one
   * REGION_RESTRICTED entry and resolve a known restricted coordinate
   * (Pyongyang) through the production lookup path: a table that lost its
   * enforcement set is boot-fatal here instead of silently permissive
   * over restricted airspace. */
  bool restricted_found = false;
  for (int i = 0; i < REGION_ENTRY_COUNT; i++) {
    if (RE_REGION(regionLookup[i]) == REGION_RESTRICTED) {
      restricted_found = true;
      break;
    }
  }
  if (!restricted_found) {
    initialized = false;
    return false;
  }
  if (latLngToRegion(39.0392, 125.7625) != REGION_RESTRICTED) {
    initialized = false;
    return false;
  }

  return true;
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
    return 0; // Invalid resolution
  }

  /* Reject non-finite and out-of-domain coordinates before they reach
   * the double->int casts in _hex2dToCoordIJK. Casting a NaN to int is
   * undefined behaviour and, critically, is NOT consistent across our
   * platforms: x86 SSE yields INT_MIN (caught by the later range check,
   * so latLngToH3 returns 0), while Cortex-M VCVT flushes NaN to 0,
   * which walks straight through to a *valid-looking* index near face 0
   * and would silently select a real LoRaWAN region. A GPS module that
   * reports a lost fix as NaN must not be able to pick a frequency
   * plan. Bounds are checked too so that a corrupt NMEA parse cannot
   * alias onto a legitimate cell. */
  if (!isfinite(lat) || !isfinite(lng) ||
      lat < -90.0 || lat > 90.0 || lng < -180.0 || lng > 180.0) {
    return 0;
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

  /* H3-1: the table now carries resolution in its key, so probe
   * most-specific-first: res 3, then res 2, then res 1 — first hit wins.
   * (Previously the device always built a 3-digit key, which made every
   * res-1/res-2 table entry unmatchable: 75.55% global Unknown rate.) */
  int d1 = H3_GET_INDEX_DIGIT(h3, 1);
  int d2 = (res >= 2) ? H3_GET_INDEX_DIGIT(h3, 2) : 0;
  int d3 = (res >= 3) ? H3_GET_INDEX_DIGIT(h3, 3) : 0;
  uint16_t k3 = (uint16_t)(d1 * 64 + d2 * 8 + d3);
  uint16_t k2 = (uint16_t)(d1 * 8 + d2);
  uint16_t k1 = (uint16_t)d1;

  RegionId region = INVALID_REGION;
  if (res >= 3)
    region = findRegion(baseCell, 3, k3);
  if (region == INVALID_REGION && res >= 2)
    region = findRegion(baseCell, 2, k2);
  if (region == INVALID_REGION && res >= 1)
    region = findRegion(baseCell, 1, k1);

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
const char *getRegionName(RegionId regionId) {
  /* H3-8: restricted territories have no entry in regionNames; name them
   * honestly. Enforcement keys on the ID (REGION_RESTRICTED), not the
   * name — REGION_UNKNOWN (0) is a different policy (keep current region
   * and transmit) and must never be conflated with restricted. */
  if (regionId == REGION_RESTRICTED) {
    return "RESTRICTED";
  }
  /* H3-7: clamp to the actual array size. Was hardcoded to 12, which
   * mislabeled CN470 (13), EU433 (14) and CD900-1A (15) as "Unknown"
   * and broke the region-name strcmp mapping in the firmware. */
  if (regionId >= (RegionId)(sizeof(regionNames) / sizeof(regionNames[0]))) {
    return regionNames[0]; // "Unknown"
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

  /* If the conversion failed there is no origin to ring-search from.
   * Without this, h3GetRing(0, k, ...) would walk baseCellNeighbors[0]
   * off an all-zero index and hand back meaningless cells. */
  if (h3 == 0) {
    return result;
  }

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
#define RING_BUFFER_SIZE 42 // 6*6=36, +6 safety margin

  if (maxRings > MAX_SUPPORTED_RINGS) {
    maxRings = MAX_SUPPORTED_RINGS; // Clamp to safe maximum
  }

  // Search rings 1 to maxRings
  for (int k = 1; k <= maxRings && result.numRegions < MAX_NEAREST_REGIONS; k++) {
    int ringSize = 6 * k;

    // BOUNDS CHECK: Verify ring size won't overflow buffer
    if (ringSize > RING_BUFFER_SIZE) {
      break; // Stop if ring would overflow - prevents corruption
    }

    H3Index ringCells[RING_BUFFER_SIZE];

    /* Get the ring of cells at distance k. H3-6(1): h3GetRing returns
     * the actual cell count (negative on pentagon failure) — iterate
     * that count, never the nominal 6*k. */
    int ringCount = h3GetRing(h3, k, ringCells);
    if (ringCount > 0) {
      // Check all cells in this ring
      for (int i = 0; i < ringCount && result.numRegions < MAX_NEAREST_REGIONS; i++) {
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
            /* Approximate distance: res-3 hexagons have ~69 km
             * mean edge length, so centre-to-centre ring spacing
             * is ~120 km (sqrt(3) * 69). The old 65 km/ring
             * figure underestimated by ~45% — a "3 rings
             * (~195 km)" claim was really ~360 km. */
            result.regions[result.numRegions].distanceKm = k * 120.0;
            result.numRegions++;
          }
        }
      }
    }

    // If we found at least one region in this ring, stop searching
    // (all regions in same ring are roughly equidistant)
    if (result.numRegions > 0)
      break;
  }

  return result;
}

// Removed unused helper function
