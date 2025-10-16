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
 * Face ijk orientation for transformations
 */
typedef struct {
    int face;             // Face number
    CoordIJK translate;   // Translation coordinates
    int ccwRot60;         // Number of 60° CCW rotations
} FaceOrientIJK;

/**
 * Unit cube face IJK coordinate structure
 */
typedef struct {
    int baseCell;         // Base cell number
    FaceIJK faceIJK;      // Face and coordinates
} BaseCellOrient;

/**
 * Overage type for indicating coordinate overage conditions
 */
typedef int Overage;

// Overage adjustment function
Overage _adjustOverageClassII(FaceIJK *fijk, int res, int pentLeading4, int substrate);

/**
 * Convert spherical coordinates (lat/lng) to face and ijk coordinates
 * 
 * @param lat Latitude in radians
 * @param lng Longitude in radians
 * @param res Resolution
 * @param fijk Output FaceIJK structure
 */
void geoToFaceIjk(double lat, double lng, int res, FaceIJK *fijk);

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
void _upAp3(CoordIJK *ijk);
void _downAp3(CoordIJK *ijk);
void _downAp3r(CoordIJK *ijk);
void _upAp7(CoordIJK *ijk);
void _upAp7r(CoordIJK *ijk);
void _downAp7(CoordIJK *ijk);
void _downAp7r(CoordIJK *ijk);

// Vector operations
void _ijkToHex2d(const CoordIJK *ijk, Vec2d *h);
void _hex2dToCoordIJK(const Vec2d *v, CoordIJK *h);
double _v2dMag(const Vec2d *v);
bool _v2dAlmostEquals(const Vec2d *v1, const Vec2d *v2);
void _v2dIntersect(const Vec2d *p1, const Vec2d *p2, const Vec2d *q1, const Vec2d *q2, Vec2d *r);

// Geo conversion helpers
void _geoToVec3d(const LatLng *geo, Vec3d *vec);
double _pointSquareDist(const Vec3d *p1, const Vec3d *p2);
void _geoToClosestFace(const LatLng *g, int *face, double *sqd);
void _geoToHex2d(const LatLng *g, int res, int *face, Vec2d *v);
void _geoToFaceIjk(const LatLng *g, int res, FaceIJK *h);
double _geoAzimuthRads(const LatLng *p1, const LatLng *p2);
void _geoAzDistanceRads(const LatLng *p1, double az, double distance, LatLng *p2);
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
int h3GetRing(H3Index origin, int k, H3Index *out);

#endif /* H3LITE_FACEIJK_H */
