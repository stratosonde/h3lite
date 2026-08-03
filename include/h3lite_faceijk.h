/*
 * H3Lite - FaceIJK functions
 * Simplified version of coordinate transformations
 */

#ifndef H3LITE_FACEIJK_H
#define H3LITE_FACEIJK_H

#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include "h3lite.h"
#include "h3lite_constants.h"

#ifdef __cplusplus
extern "C" {
#endif

// LatLng is defined in h3lite.h

/**
 * 2D cartesian coordinate
 */
typedef struct {
    double x;
    double y;
} Vec2d;

/**
 * 3D cartesian coordinate
 */
typedef struct {
    double x;
    double y;
    double z;
} Vec3d;

/**
 * IJK hexagonal grid coordinates
 */
typedef struct {
    int i;      // i coordinate
    int j;      // j coordinate
    int k;      // k coordinate (where i + j + k = 0)
} CoordIJK;

/**
 * Face number and ijk coordinates on that face
 */
typedef struct {
    int face;         // Face number
    CoordIJK coord;   // IJK coordinates on that face
} FaceIJK;

/**
 * Narrow (int8_t) mirrors of CoordIJK/FaceIJK used only by the static
 * base-cell tables. Every value they hold fits in a signed byte
 * (face 0-19, ijk 0-2), so storing them as int is a 4x flash waste.
 * The brace layout is identical to CoordIJK/FaceIJK, so the generated
 * table initializers are unchanged.
 */
typedef struct {
    int8_t i;
    int8_t j;
    int8_t k;
} CoordIJK8;

typedef struct {
    int8_t face;
    CoordIJK8 coord;
} FaceIJK8;

/**
 * Base cell home face/coords, pentagon flag and cw-offset faces.
 * Declared here (not privately in each .c) so that every translation
 * unit agrees on the layout. Previously h3lite_neighbor.c declared this
 * as `extern const FaceIJK baseCellData[]`, a 16-byte stride against a
 * 28-byte definition: an ODR violation that silently read wrong data.
 */
typedef struct {
    FaceIJK8 homeFijk;      // Home face and IJK
    int8_t isPentagon;      // Is this a pentagon
    int8_t cwOffsetPent[2]; // CW offset faces for pentagons (-1 if none)
} BaseCellData;

/** @brief base cell at a given ijk and required rotations into its system */
typedef struct {
    uint8_t baseCell;  // base cell number (0-121)
    int8_t ccwRot60;   // number of ccw 60 degree rotations
} BaseCellRotation;

extern const BaseCellData baseCellData[NUM_BASE_CELLS];
extern const uint8_t baseCellNeighbors[NUM_BASE_CELLS][7];
extern const int8_t baseCellNeighbor60CCWRots[NUM_BASE_CELLS][7];

/**
 * Unit cube face IJK coordinate structure
 */
typedef struct {
    int baseCell;         // Base cell number
    FaceIJK faceIJK;      // Face and coordinates
} BaseCellOrient;

/**
 * Convert FaceIJK coordinates to H3 index
 * 
 * @param fijk FaceIJK coordinates
 * @param res Resolution
 * @return H3 index
 */
uint64_t faceIjkToH3(const FaceIJK *fijk, int res);

// Helper functions for IJK coordinate manipulation
void _setIJK(CoordIJK *ijk, int i, int j, int k);
void _ijkNormalize(CoordIJK *ijk);
void _ijkAdd(const CoordIJK *a, const CoordIJK *b, CoordIJK *sum);
void _ijkSub(const CoordIJK *a, const CoordIJK *b, CoordIJK *diff);
void _ijkScale(CoordIJK *ijk, int factor);
void _ijkRotate60ccw(CoordIJK *ijk);
void _ijkRotate60cw(CoordIJK *ijk);

// Aperture operations
void _upAp7(CoordIJK *ijk);
void _upAp7r(CoordIJK *ijk);
void _downAp7(CoordIJK *ijk);
void _downAp7r(CoordIJK *ijk);

// Vector operations
void _hex2dToCoordIJK(const Vec2d *v, CoordIJK *h);

// Geo conversion helpers
void _geoToVec3d(const LatLng *geo, Vec3d *vec);
double _pointSquareDist(const Vec3d *p1, const Vec3d *p2);
void _geoToClosestFace(const LatLng *g, int *face, double *sqd);
void _geoToHex2d(const LatLng *g, int res, int *face, Vec2d *v);
void _geoToFaceIjk(const LatLng *g, int res, FaceIJK *h);
double _geoAzimuthRads(const LatLng *p1, const LatLng *p2);
double _posAngleRads(double rads);

// Base cell functions
int _faceIjkToBaseCell(const FaceIJK *h);
int _faceIjkToBaseCellCCWrot60(const FaceIJK *h);
bool _isBaseCellPentagon(int baseCell);
bool _baseCellIsCwOffset(int baseCell, int testFace);
void _initBaseCellFaceOrient(int baseCell, int face, FaceIJK *result);

// Resolution class check
bool isResolutionClassIII(int res);

// H3 neighbor and ring functions
int h3NeighborRotations(H3Index origin, Direction dir, int *rotations, H3Index *out);
/* Returns the number of cells written to out (6*k, or 1 if k == 0) on
 * success; negative if a pentagon was encountered. Callers must iterate
 * the returned count, not the nominal 6*k. */
int h3GetRing(H3Index origin, int k, H3Index *out);

#ifdef __cplusplus
}
#endif

#endif /* H3LITE_FACEIJK_H */
