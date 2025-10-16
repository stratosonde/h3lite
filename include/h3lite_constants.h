/*
 * H3Lite - Constants
 * Minimized set of constants for H3 calculations
 */

#ifndef H3LITE_CONSTANTS_H
#define H3LITE_CONSTANTS_H

#include <math.h>
#include <stdint.h>

// Resolution-specific constants
#define H3LITE_MAX_RESOLUTION 4  // Only supporting up to resolution 4
#define H3LITE_TARGET_RESOLUTION 3  // Primary resolution to use for region lookup

// Icosahedron constants
#define NUM_ICOSA_FACES 20
#define NUM_PENTAGONS 12
#define NUM_BASE_CELLS 122
#define NUM_HEX_VERTS 6
#define NUM_PENT_VERTS 5

// Direction constants
#define IJ 1
#define JK 2
#define KI 3

// Overage types
#define NO_OVERAGE 0
#define FACE_EDGE 1
#define NEW_FACE 2

// Base cell constants
#define H3LITE_NUM_BASE_CELLS 122
#define INVALID_BASE_CELL 127  // Invalid base cell marker
#define INVALID_FACE -1
#define MAX_FACE_COORD 2  // Maximum i, j, k coordinate on an icosahedron face

// H3 index bit layout constants
#define H3_NUM_BITS 64
#define H3_MAX_OFFSET 63
#define H3_MODE_OFFSET 59
#define H3_BC_OFFSET 45
#define H3_RES_OFFSET 52
#define H3_PER_DIGIT_OFFSET 3

// H3 bit masks
#define H3_BC_MASK ((uint64_t)(127) << H3_BC_OFFSET)
#define H3_RES_MASK ((uint64_t)(15) << H3_RES_OFFSET)
#define H3_DIGIT_MASK ((uint64_t)(7))

// H3 bit mask negatives (for bit clearing)
#define H3_BC_MASK_NEGATIVE (~H3_BC_MASK)
#define H3_RES_MASK_NEGATIVE (~H3_RES_MASK)
#define H3_DIGIT_MASK_NEGATIVE (~H3_DIGIT_MASK)
#define H3_HIGH_BIT_MASK ((uint64_t)(1) << H3_MAX_OFFSET)
#define H3_HIGH_BIT_MASK_NEGATIVE (~H3_HIGH_BIT_MASK)
#define H3_MODE_MASK ((uint64_t)(15) << H3_MODE_OFFSET)
#define H3_MODE_MASK_NEGATIVE (~H3_MODE_MASK)

// H3 bit manipulation macros
#define H3_GET_HIGH_BIT(h3) ((int)((((h3)&H3_HIGH_BIT_MASK) >> H3_MAX_OFFSET)))
#define H3_SET_HIGH_BIT(h3, v) (h3) = (((h3)&H3_HIGH_BIT_MASK_NEGATIVE) | (((uint64_t)(v)) << H3_MAX_OFFSET))
#define H3_GET_MODE(h3) ((int)((((h3)&H3_MODE_MASK) >> H3_MODE_OFFSET)))
#define H3_SET_MODE(h3, v) (h3) = (((h3)&H3_MODE_MASK_NEGATIVE) | (((uint64_t)(v)) << H3_MODE_OFFSET))
#define H3_GET_BASE_CELL(h3) ((int)((((h3)&H3_BC_MASK) >> H3_BC_OFFSET)))
#define H3_SET_BASE_CELL(h3, bc) (h3) = (((h3)&H3_BC_MASK_NEGATIVE) | (((uint64_t)(bc)) << H3_BC_OFFSET))
#define H3_GET_RESOLUTION(h3) ((int)((((h3)&H3_RES_MASK) >> H3_RES_OFFSET)))
#define H3_SET_RESOLUTION(h3, res) (h3) = (((h3)&H3_RES_MASK_NEGATIVE) | (((uint64_t)(res)) << H3_RES_OFFSET))
#define H3_GET_INDEX_DIGIT(h3, res) ((int)((((h3) >> ((15 - (res)) * H3_PER_DIGIT_OFFSET)) & H3_DIGIT_MASK)))
#define H3_SET_INDEX_DIGIT(h3, res, digit) (h3) = (((h3) & ~((H3_DIGIT_MASK << ((15 - (res)) * H3_PER_DIGIT_OFFSET)))) | (((uint64_t)(digit)) << ((15 - (res)) * H3_PER_DIGIT_OFFSET)))

// H3 index modes
#define H3_CELL_MODE 1

// Direction enum type (matching original H3 library)
typedef int Direction;

// Direction digits (matching H3 reference library)
#define CENTER_DIGIT 0
#define K_AXES_DIGIT 1
#define J_AXES_DIGIT 2
#define JK_AXES_DIGIT 3  // J_AXES_DIGIT | K_AXES_DIGIT
#define I_AXES_DIGIT 4
#define IK_AXES_DIGIT 5  // I_AXES_DIGIT | K_AXES_DIGIT
#define IJ_AXES_DIGIT 6  // I_AXES_DIGIT | J_AXES_DIGIT

// Invalid digit marker
#define INVALID_DIGIT 7
#define NUM_DIGITS 7

/**
 * H3 index with mode 0, res 0, base cell 0, and 7 for all index digits.
 * Typically used to initialize the creation of an H3 cell index, which
 * expects all direction digits to be 7 beyond the cell's resolution.
 */
#define H3_INIT (UINT64_C(35184372088831))

// Math constants
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
#define M_PI_180 (M_PI / 180.0)
#define M_180_PI (180.0 / M_PI)
#define EPSILON 1e-10
#define M_AP7_ROT_RADS 0.333473172251832115336090755351601070065900389
#define M_SQRT3_2 0.8660254037844386467637231707
#define M_SQRT7 2.6457513110645905905016157536392604257102
#define M_RSQRT7 0.37796447300922722721451653623418006081576
#define M_ONESEVENTH 0.14285714285714285714285714285714
#define M_RSIN60 1.1547005383792515290182975610039149112953
#define INV_RES0_U_GNOMONIC 2.61803398874989588842
#define RES0_U_GNOMONIC 2.6181773447340
#define M_SIN60 0.8660254037844386 // sin(60 degrees)
#define M_COS60 0.5000000000000000 // cos(60 degrees)

// Earth radius in kilometers
#define EARTH_RADIUS_KM 6371.0

// Lookup specific constants
#define INVALID_REGION 0
#define MAX_REGIONS 12  // Based on the approximately 12 regions mentioned

#endif /* H3LITE_CONSTANTS_H */
