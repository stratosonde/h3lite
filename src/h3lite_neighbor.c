/*
 * H3Lite - Neighbor and Ring Functions
 * Ported from H3 reference implementation (algos.c)
 * Copyright 2016-2021 Uber Technologies, Inc.
 * Licensed under the Apache License, Version 2.0
 */

#include "../include/h3lite_faceijk.h"
#include "../include/h3lite_constants.h"
#include "../include/h3lite.h"
#include <stdlib.h>

// External declarations
extern const Direction DIRECTIONS[6];
extern const Direction NEXT_RING_DIRECTION;
extern const Direction NEW_DIGIT_II[7][7];
extern const Direction NEW_ADJUSTMENT_II[7][7];
extern const Direction NEW_DIGIT_III[7][7];
extern const Direction NEW_ADJUSTMENT_III[7][7];

// Need access to basecell data
extern const int baseCellNeighbors[NUM_BASE_CELLS][7];
extern const int baseCellNeighbor60CCWRots[NUM_BASE_CELLS][7];
extern const FaceIJK baseCellData[NUM_BASE_CELLS];

// Forward declarations
static Direction _rotate60ccw(Direction digit);
static Direction _rotate60cw(Direction digit);
static uint64_t _h3Rotate60ccw(uint64_t h);
static uint64_t _h3Rotate60cw(uint64_t h);
static uint64_t _h3RotatePent60ccw(uint64_t h);
static Direction _h3LeadingNonZeroDigit(uint64_t h);

/**
 * Rotate a direction 60 degrees counter-clockwise
 */
static Direction _rotate60ccw(Direction digit) {
    switch(digit) {
        case K_AXES_DIGIT: return IK_AXES_DIGIT;
        case IK_AXES_DIGIT: return I_AXES_DIGIT;
        case I_AXES_DIGIT: return IJ_AXES_DIGIT;
        case IJ_AXES_DIGIT: return J_AXES_DIGIT;
        case J_AXES_DIGIT: return JK_AXES_DIGIT;
        case JK_AXES_DIGIT: return K_AXES_DIGIT;
        default: return digit;
    }
}

/**
 * Rotate a direction 60 degrees clockwise
 */
static Direction _rotate60cw(Direction digit) {
    switch(digit) {
        case K_AXES_DIGIT: return JK_AXES_DIGIT;
        case JK_AXES_DIGIT: return J_AXES_DIGIT;
        case J_AXES_DIGIT: return IJ_AXES_DIGIT;
        case IJ_AXES_DIGIT: return I_AXES_DIGIT;
        case I_AXES_DIGIT: return IK_AXES_DIGIT;
        case IK_AXES_DIGIT: return K_AXES_DIGIT;
        default: return digit;
    }
}

/**
 * Get the leading non-zero digit of an H3 index
 */
static Direction _h3LeadingNonZeroDigit(uint64_t h) {
    int res = H3_GET_RESOLUTION(h);
    for (int r = 1; r <= res; r++) {
        Direction digit = H3_GET_INDEX_DIGIT(h, r);
        if (digit != 0) {
            return digit;
        }
    }
    return CENTER_DIGIT;
}

/**
 * Rotate an H3 index 60 degrees counter-clockwise
 */
static uint64_t _h3Rotate60ccw(uint64_t h) {
    int res = H3_GET_RESOLUTION(h);
    for (int r = 1; r <= res; r++) {
        H3_SET_INDEX_DIGIT(h, r, _rotate60ccw(H3_GET_INDEX_DIGIT(h, r)));
    }
    return h;
}

/**
 * Rotate an H3 index 60 degrees clockwise
 */
static uint64_t _h3Rotate60cw(uint64_t h) {
    int res = H3_GET_RESOLUTION(h);
    for (int r = 1; r <= res; r++) {
        H3_SET_INDEX_DIGIT(h, r, _rotate60cw(H3_GET_INDEX_DIGIT(h, r)));
    }
    return h;
}

/**
 * Rotate a pentagon H3 index 60 degrees counter-clockwise
 */
static uint64_t _h3RotatePent60ccw(uint64_t h) {
    int foundFirstNonZeroDigit = 0;
    int res = H3_GET_RESOLUTION(h);
    
    for (int r = 1; r <= res; r++) {
        Direction digit = H3_GET_INDEX_DIGIT(h, r);
        H3_SET_INDEX_DIGIT(h, r, _rotate60ccw(digit));
        
        if (!foundFirstNonZeroDigit && digit != 0) {
            foundFirstNonZeroDigit = 1;
            if (_h3LeadingNonZeroDigit(h) == K_AXES_DIGIT) {
                h = _h3Rotate60ccw(h);
            }
        }
    }
    return h;
}

/**
 * Check if an H3 index is a pentagon
 */
static bool _h3IsPentagon(H3Index h) {
    int baseCell = H3_GET_BASE_CELL(h);
    return _isBaseCellPentagon(baseCell);
}

/**
 * Returns the hexagon index neighboring the origin, in the direction dir.
 * Ported from H3 reference implementation
 *
 * @param origin Origin index
 * @param dir Direction to move in
 * @param rotations Number of ccw rotations to perform
 * @param out H3Index of the specified neighbor if successful
 * @return 0 on success, non-zero on failure
 */
int h3NeighborRotations(H3Index origin, Direction dir, int *rotations, H3Index *out) {
    H3Index current = origin;

    if (dir < CENTER_DIGIT || dir >= INVALID_DIGIT) {
        return -1;  // Invalid direction
    }
    
    // Ensure rotations is modulo'd by 6
    *rotations = *rotations % 6;
    
    for (int i = 0; i < *rotations; i++) {
        dir = _rotate60ccw(dir);
    }

    int newRotations = 0;
    int oldBaseCell = H3_GET_BASE_CELL(current);
    if (oldBaseCell < 0 || oldBaseCell >= NUM_BASE_CELLS) {
        return -1;  // Invalid base cell
    }
    
    Direction oldLeadingDigit = _h3LeadingNonZeroDigit(current);

    // Adjust the indexing digits and, if needed, the base cell.
    int r = H3_GET_RESOLUTION(current) - 1;
    
    while (true) {
        if (r == -1) {
            H3_SET_BASE_CELL(current, baseCellNeighbors[oldBaseCell][dir]);
            newRotations = baseCellNeighbor60CCWRots[oldBaseCell][dir];

            if (H3_GET_BASE_CELL(current) == INVALID_BASE_CELL) {
                // Adjust for the deleted k vertex at the base cell level
                H3_SET_BASE_CELL(current, baseCellNeighbors[oldBaseCell][IK_AXES_DIGIT]);
                newRotations = baseCellNeighbor60CCWRots[oldBaseCell][IK_AXES_DIGIT];

                // Perform the adjustment for the k-subsequence we're skipping over
                current = _h3Rotate60ccw(current);
                *rotations = *rotations + 1;
            }
            break;
        } else {
            Direction oldDigit = H3_GET_INDEX_DIGIT(current, r + 1);
            Direction nextDir;
            
            if (oldDigit == INVALID_DIGIT) {
                return -1;  // Invalid digit
            } else if (isResolutionClassIII(r + 1)) {
                H3_SET_INDEX_DIGIT(current, r + 1, NEW_DIGIT_II[oldDigit][dir]);
                nextDir = NEW_ADJUSTMENT_II[oldDigit][dir];
            } else {
                H3_SET_INDEX_DIGIT(current, r + 1, NEW_DIGIT_III[oldDigit][dir]);
                nextDir = NEW_ADJUSTMENT_III[oldDigit][dir];
            }

            if (nextDir != CENTER_DIGIT) {
                dir = nextDir;
                r--;
            } else {
                break;
            }
        }
    }

    int newBaseCell = H3_GET_BASE_CELL(current);
    
    if (_isBaseCellPentagon(newBaseCell)) {
        int alreadyAdjustedKSubsequence = 0;

        // Force rotation out of missing k-axes sub-sequence
        if (_h3LeadingNonZeroDigit(current) == K_AXES_DIGIT) {
            if (oldBaseCell != newBaseCell) {
                // Traversed into the deleted k subsequence of a pentagon base cell
                if (_baseCellIsCwOffset(newBaseCell, baseCellData[oldBaseCell].face)) {
                    current = _h3Rotate60cw(current);
                } else {
                    current = _h3Rotate60ccw(current);
                }
                alreadyAdjustedKSubsequence = 1;
            } else {
                // Traversed into deleted k subsequence from within same pentagon
                if (oldLeadingDigit == CENTER_DIGIT) {
                    return 1;  // Pentagon - k direction is deleted
                } else if (oldLeadingDigit == JK_AXES_DIGIT) {
                    current = _h3Rotate60ccw(current);
                    *rotations = *rotations + 1;
                } else if (oldLeadingDigit == IK_AXES_DIGIT) {
                    current = _h3Rotate60cw(current);
                    *rotations = *rotations + 5;
                } else {
                    return -1;  // Unexpected case
                }
            }
        }

        for (int i = 0; i < newRotations; i++) {
            current = _h3RotatePent60ccw(current);
        }

        // Account for differing orientation of base cells
        if (oldBaseCell != newBaseCell) {
            if (_isBaseCellPentagon(newBaseCell) && 
                (newBaseCell == 4 || newBaseCell == 117)) {  // Polar pentagons
                if (oldBaseCell != 118 && oldBaseCell != 8 &&
                    _h3LeadingNonZeroDigit(current) != JK_AXES_DIGIT) {
                    *rotations = *rotations + 1;
                }
            } else if (_h3LeadingNonZeroDigit(current) == IK_AXES_DIGIT &&
                       !alreadyAdjustedKSubsequence) {
                *rotations = *rotations + 1;
            }
        }
    } else {
        for (int i = 0; i < newRotations; i++) {
            current = _h3Rotate60ccw(current);
        }
    }

    *rotations = (*rotations + newRotations) % 6;
    *out = current;
    return 0;
}

/**
 * Returns the "hollow" ring of hexagons at exactly grid distance k from
 * the origin hexagon. Ported from H3 reference implementation.
 *
 * @param origin Origin location
 * @param k k >= 0
 * @param out Array which must be of size 6 * k (or 1 if k == 0)
 * @return 0 if successful; 1 if pentagon encountered
 */
int h3GetRing(H3Index origin, int k, H3Index *out) {
    // Short-circuit on 'identity' ring
    if (k == 0) {
        out[0] = origin;
        return 0;
    }
    
    int idx = 0;
    int rotations = 0;
    
    if (_h3IsPentagon(origin)) {
        return 1;  // Pentagon encountered
    }

    // Move to the start of the ring (k steps in NEXT_RING_DIRECTION)
    for (int ring = 0; ring < k; ring++) {
        int result = h3NeighborRotations(origin, NEXT_RING_DIRECTION, &rotations, &origin);
        if (result != 0) {
            return result;
        }

        if (_h3IsPentagon(origin)) {
            return 1;  // Pentagon encountered
        }
    }

    H3Index lastIndex = origin;
    out[idx++] = origin;

    // Traverse the ring
    for (int direction = 0; direction < 6; direction++) {
        for (int pos = 0; pos < k; pos++) {
            int result = h3NeighborRotations(origin, DIRECTIONS[direction], &rotations, &origin);
            if (result != 0) {
                return result;
            }

            // Skip the very last index (it was already added as first)
            if (pos != k - 1 || direction != 5) {
                out[idx++] = origin;

                if (_h3IsPentagon(origin)) {
                    return 1;  // Pentagon encountered
                }
            }
        }
    }

    // Check that we completed the ring properly
    if (lastIndex != origin) {
        return 1;  // Pentagon distortion
    }

    return 0;
}
