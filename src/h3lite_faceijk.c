/*
 * H3Lite - FaceIJK implementation
 * Accurate implementation of H3 coordinate transformations
 */

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include "../include/h3lite_faceijk.h"
#include "../include/h3lite_constants.h"

// Forward declaration from h3lite.c
extern const int baseCellTable[12][12];

/**
 * Icosahedron face centers in x/y/z on the unit sphere (from H3)
 */
static const Vec3d faceCenterPoint[NUM_ICOSA_FACES] = {
    {0.2199307791404606, 0.6583691780274996, 0.7198475378926182},     // face  0
    {-0.2139234834501421, 0.1478171829550703, 0.9656017935214205},    // face  1
    {0.1092625278784797, -0.4811951572873210, 0.8697775121287253},    // face  2
    {0.7428567301586791, -0.3593941678278028, 0.5648005936517033},    // face  3
    {0.8112534709140969, 0.3448953237639384, 0.4721387736413930},     // face  4
    {-0.1055498149613921, 0.9794457296411413, 0.1718874610009365},    // face  5
    {-0.8075407579970092, 0.1533552485898818, 0.5695261994882688},    // face  6
    {-0.2846148069787907, -0.8644080972654206, 0.4144792552473539},   // face  7
    {0.7405621473854482, -0.6673299564565524, -0.0789837646326737},   // face  8
    {0.8512303986474293, 0.4722343788582681, -0.2289137388687808},    // face  9
    {-0.7405621473854481, 0.6673299564565524, 0.0789837646326737},    // face 10
    {-0.8512303986474292, -0.4722343788582682, 0.2289137388687808},   // face 11
    {0.1055498149613919, -0.9794457296411413, -0.1718874610009365},   // face 12
    {0.8075407579970092, -0.1533552485898819, -0.5695261994882688},   // face 13
    {0.2846148069787908, 0.8644080972654204, -0.4144792552473539},    // face 14
    {-0.7428567301586791, 0.3593941678278027, -0.5648005936517033},   // face 15
    {-0.8112534709140971, -0.3448953237639382, -0.4721387736413930},  // face 16
    {-0.2199307791404607, -0.6583691780274996, -0.7198475378926182},  // face 17
    {0.2139234834501420, -0.1478171829550704, -0.9656017935214205},   // face 18
    {-0.1092625278784796, 0.4811951572873210, -0.8697775121287253},   // face 19
};

/**
 * Icosahedron face ijk axes as azimuth in radians (from H3)
 */
static const double faceAxesAzRadsCII[NUM_ICOSA_FACES][3] = {
    {5.619958268523939882, 3.525563166130744542, 1.431168063737548730},  // face  0
    {5.760339081714187279, 3.665943979320991689, 1.571548876927796127},  // face  1
    {0.780213654393430055, 4.969003859179821079, 2.874608756786625655},  // face  2
    {0.430469363979999913, 4.619259568766391033, 2.524864466373195467},  // face  3
    {6.130269123335111400, 4.035874020941915804, 1.941478918548720291},  // face  4
    {2.692877706530642877, 0.598482604137447119, 4.787272808923838195},  // face  5
    {2.982963003477243874, 0.888567901084048369, 5.077358105870439581},  // face  6
    {3.532912002790141181, 1.438516900396945656, 5.627307105183336758},  // face  7
    {3.494305004259568154, 1.399909901866372864, 5.588700106652763840},  // face  8
    {3.003214169499538391, 0.908819067106342928, 5.097609271892733906},  // face  9
    {5.930472956509811562, 3.836077854116615875, 1.741682751723420374},  // face 10
    {0.138378484090254847, 4.327168688876645809, 2.232773586483450311},  // face 11
    {0.448714947059150361, 4.637505151845541521, 2.543110049452346120},  // face 12
    {0.158629650112549365, 4.347419854898940135, 2.253024752505744869},  // face 13
    {5.891865957979238535, 3.797470855586042958, 1.703075753192847583},  // face 14
    {2.711123289609793325, 0.616728187216597771, 4.805518392002988683},  // face 15
    {3.294508837434268316, 1.200113735041072948, 5.388903939827463911},  // face 16
    {3.804819692245439833, 1.710424589852244509, 5.899214794638635174},  // face 17
    {3.664438879055192436, 1.570043776661997111, 5.758833981448388027},  // face 18
    {2.361378999196363184, 0.266983896803167583, 4.455774101589558636},  // face 19
};

/**
 * Face neighbor data for overage handling
 */
static const FaceOrientIJK faceNeighbors[NUM_ICOSA_FACES][4] = {
    {{0, {0, 0, 0}, 0}, {4, {2, 0, 2}, 1}, {1, {2, 2, 0}, 5}, {5, {0, 2, 2}, 3}},
    {{1, {0, 0, 0}, 0}, {0, {2, 0, 2}, 1}, {2, {2, 2, 0}, 5}, {6, {0, 2, 2}, 3}},
    {{2, {0, 0, 0}, 0}, {1, {2, 0, 2}, 1}, {3, {2, 2, 0}, 5}, {7, {0, 2, 2}, 3}},
    {{3, {0, 0, 0}, 0}, {2, {2, 0, 2}, 1}, {4, {2, 2, 0}, 5}, {8, {0, 2, 2}, 3}},
    {{4, {0, 0, 0}, 0}, {3, {2, 0, 2}, 1}, {0, {2, 2, 0}, 5}, {9, {0, 2, 2}, 3}},
    {{5, {0, 0, 0}, 0}, {10, {2, 2, 0}, 3}, {14, {2, 0, 2}, 3}, {0, {0, 2, 2}, 3}},
    {{6, {0, 0, 0}, 0}, {11, {2, 2, 0}, 3}, {10, {2, 0, 2}, 3}, {1, {0, 2, 2}, 3}},
    {{7, {0, 0, 0}, 0}, {12, {2, 2, 0}, 3}, {11, {2, 0, 2}, 3}, {2, {0, 2, 2}, 3}},
    {{8, {0, 0, 0}, 0}, {13, {2, 2, 0}, 3}, {12, {2, 0, 2}, 3}, {3, {0, 2, 2}, 3}},
    {{9, {0, 0, 0}, 0}, {14, {2, 2, 0}, 3}, {13, {2, 0, 2}, 3}, {4, {0, 2, 2}, 3}},
    {{10, {0, 0, 0}, 0}, {5, {2, 2, 0}, 3}, {6, {2, 0, 2}, 3}, {15, {0, 2, 2}, 3}},
    {{11, {0, 0, 0}, 0}, {6, {2, 2, 0}, 3}, {7, {2, 0, 2}, 3}, {16, {0, 2, 2}, 3}},
    {{12, {0, 0, 0}, 0}, {7, {2, 2, 0}, 3}, {8, {2, 0, 2}, 3}, {17, {0, 2, 2}, 3}},
    {{13, {0, 0, 0}, 0}, {8, {2, 2, 0}, 3}, {9, {2, 0, 2}, 3}, {18, {0, 2, 2}, 3}},
    {{14, {0, 0, 0}, 0}, {9, {2, 2, 0}, 3}, {5, {2, 0, 2}, 3}, {19, {0, 2, 2}, 3}},
    {{15, {0, 0, 0}, 0}, {16, {2, 0, 2}, 1}, {19, {2, 2, 0}, 5}, {10, {0, 2, 2}, 3}},
    {{16, {0, 0, 0}, 0}, {17, {2, 0, 2}, 1}, {15, {2, 2, 0}, 5}, {11, {0, 2, 2}, 3}},
    {{17, {0, 0, 0}, 0}, {18, {2, 0, 2}, 1}, {16, {2, 2, 0}, 5}, {12, {0, 2, 2}, 3}},
    {{18, {0, 0, 0}, 0}, {19, {2, 0, 2}, 1}, {17, {2, 2, 0}, 5}, {13, {0, 2, 2}, 3}},
    {{19, {0, 0, 0}, 0}, {15, {2, 0, 2}, 1}, {18, {2, 2, 0}, 5}, {14, {0, 2, 2}, 3}}
};

/**
 * Maximum dimension for each resolution in Class II
 */
static const int maxDimByCIIres[] = {
    2, -1, 14, -1, 98, -1, 686, -1, 4802, -1, 33614, -1, 235298, -1, 1647086, -1, 11529602
};

/**
 * Unit scale for each resolution in Class II
 */
static const int unitScaleByCIIres[] = {
    1, -1, 7, -1, 49, -1, 343, -1, 2401, -1, 16807, -1, 117649, -1, 823543, -1, 5764801
};

/**
 * Directions used for traversing a hexagonal ring counterclockwise
 */
const Direction DIRECTIONS[6] = {J_AXES_DIGIT, JK_AXES_DIGIT,
                                        K_AXES_DIGIT, IK_AXES_DIGIT,
                                        I_AXES_DIGIT, IJ_AXES_DIGIT};

/**
 * Direction used for traversing to the next outward hexagonal ring
 */
const Direction NEXT_RING_DIRECTION = I_AXES_DIGIT;

/**
 * New digit when traversing along class II grids
 */
const Direction NEW_DIGIT_II[7][7] = {
    {CENTER_DIGIT, K_AXES_DIGIT, J_AXES_DIGIT, JK_AXES_DIGIT, I_AXES_DIGIT,
     IK_AXES_DIGIT, IJ_AXES_DIGIT},
    {K_AXES_DIGIT, I_AXES_DIGIT, JK_AXES_DIGIT, IJ_AXES_DIGIT, IK_AXES_DIGIT,
     J_AXES_DIGIT, CENTER_DIGIT},
    {J_AXES_DIGIT, JK_AXES_DIGIT, K_AXES_DIGIT, I_AXES_DIGIT, IJ_AXES_DIGIT,
     CENTER_DIGIT, IK_AXES_DIGIT},
    {JK_AXES_DIGIT, IJ_AXES_DIGIT, I_AXES_DIGIT, IK_AXES_DIGIT, CENTER_DIGIT,
     K_AXES_DIGIT, J_AXES_DIGIT},
    {I_AXES_DIGIT, IK_AXES_DIGIT, IJ_AXES_DIGIT, CENTER_DIGIT, J_AXES_DIGIT,
     JK_AXES_DIGIT, K_AXES_DIGIT},
    {IK_AXES_DIGIT, J_AXES_DIGIT, CENTER_DIGIT, K_AXES_DIGIT, JK_AXES_DIGIT,
     IJ_AXES_DIGIT, I_AXES_DIGIT},
    {IJ_AXES_DIGIT, CENTER_DIGIT, IK_AXES_DIGIT, J_AXES_DIGIT, K_AXES_DIGIT,
     I_AXES_DIGIT, JK_AXES_DIGIT}};

/**
 * New traversal direction when traversing along class II grids
 */
const Direction NEW_ADJUSTMENT_II[7][7] = {
    {CENTER_DIGIT, CENTER_DIGIT, CENTER_DIGIT, CENTER_DIGIT, CENTER_DIGIT,
     CENTER_DIGIT, CENTER_DIGIT},
    {CENTER_DIGIT, K_AXES_DIGIT, CENTER_DIGIT, K_AXES_DIGIT, CENTER_DIGIT,
     IK_AXES_DIGIT, CENTER_DIGIT},
    {CENTER_DIGIT, CENTER_DIGIT, J_AXES_DIGIT, JK_AXES_DIGIT, CENTER_DIGIT,
     CENTER_DIGIT, J_AXES_DIGIT},
    {CENTER_DIGIT, K_AXES_DIGIT, JK_AXES_DIGIT, JK_AXES_DIGIT, CENTER_DIGIT,
     CENTER_DIGIT, CENTER_DIGIT},
    {CENTER_DIGIT, CENTER_DIGIT, CENTER_DIGIT, CENTER_DIGIT, I_AXES_DIGIT,
     I_AXES_DIGIT, IJ_AXES_DIGIT},
    {CENTER_DIGIT, IK_AXES_DIGIT, CENTER_DIGIT, CENTER_DIGIT, I_AXES_DIGIT,
     IK_AXES_DIGIT, CENTER_DIGIT},
    {CENTER_DIGIT, CENTER_DIGIT, J_AXES_DIGIT, CENTER_DIGIT, IJ_AXES_DIGIT,
     CENTER_DIGIT, IJ_AXES_DIGIT}};

/**
 * New traversal direction when traversing along class III grids
 */
const Direction NEW_DIGIT_III[7][7] = {
    {CENTER_DIGIT, K_AXES_DIGIT, J_AXES_DIGIT, JK_AXES_DIGIT, I_AXES_DIGIT,
     IK_AXES_DIGIT, IJ_AXES_DIGIT},
    {K_AXES_DIGIT, J_AXES_DIGIT, JK_AXES_DIGIT, I_AXES_DIGIT, IK_AXES_DIGIT,
     IJ_AXES_DIGIT, CENTER_DIGIT},
    {J_AXES_DIGIT, JK_AXES_DIGIT, I_AXES_DIGIT, IK_AXES_DIGIT, IJ_AXES_DIGIT,
     CENTER_DIGIT, K_AXES_DIGIT},
    {JK_AXES_DIGIT, I_AXES_DIGIT, IK_AXES_DIGIT, IJ_AXES_DIGIT, CENTER_DIGIT,
     K_AXES_DIGIT, J_AXES_DIGIT},
    {I_AXES_DIGIT, IK_AXES_DIGIT, IJ_AXES_DIGIT, CENTER_DIGIT, K_AXES_DIGIT,
     J_AXES_DIGIT, JK_AXES_DIGIT},
    {IK_AXES_DIGIT, IJ_AXES_DIGIT, CENTER_DIGIT, K_AXES_DIGIT, J_AXES_DIGIT,
     JK_AXES_DIGIT, I_AXES_DIGIT},
    {IJ_AXES_DIGIT, CENTER_DIGIT, K_AXES_DIGIT, J_AXES_DIGIT, JK_AXES_DIGIT,
     I_AXES_DIGIT, IK_AXES_DIGIT}};

/**
 * New traversal direction when traversing along class III grids
 */
const Direction NEW_ADJUSTMENT_III[7][7] = {
    {CENTER_DIGIT, CENTER_DIGIT, CENTER_DIGIT, CENTER_DIGIT, CENTER_DIGIT,
     CENTER_DIGIT, CENTER_DIGIT},
    {CENTER_DIGIT, K_AXES_DIGIT, CENTER_DIGIT, JK_AXES_DIGIT, CENTER_DIGIT,
     K_AXES_DIGIT, CENTER_DIGIT},
    {CENTER_DIGIT, CENTER_DIGIT, J_AXES_DIGIT, J_AXES_DIGIT, CENTER_DIGIT,
     CENTER_DIGIT, IJ_AXES_DIGIT},
    {CENTER_DIGIT, JK_AXES_DIGIT, J_AXES_DIGIT, JK_AXES_DIGIT, CENTER_DIGIT,
     CENTER_DIGIT, CENTER_DIGIT},
    {CENTER_DIGIT, CENTER_DIGIT, CENTER_DIGIT, CENTER_DIGIT, I_AXES_DIGIT,
     IK_AXES_DIGIT, I_AXES_DIGIT},
    {CENTER_DIGIT, K_AXES_DIGIT, CENTER_DIGIT, CENTER_DIGIT, IK_AXES_DIGIT,
     IK_AXES_DIGIT, CENTER_DIGIT},
    {CENTER_DIGIT, CENTER_DIGIT, IJ_AXES_DIGIT, CENTER_DIGIT, I_AXES_DIGIT,
     CENTER_DIGIT, IJ_AXES_DIGIT}};

/**
 * Icosahedron face centers in lat/lng radians
 */
const LatLng faceCenterGeo[NUM_ICOSA_FACES] = {
    {0.803582649718989942, 1.248397419617396099},
    {1.307747883455638156, 2.536945009877921159},
    {1.054751253523952054, -1.347517358900396623},
    {0.600191595538186799, -0.450603909469755746},
    {0.491715428198773866, 0.401988202911306943},
    {0.172745327415618701, 1.678146885280433686},
    {0.605929321571350690, 2.953923329812411617},
    {0.427370518328979641, -1.888876200336285401},
    {-0.079066118549212831, -0.733429513380867741},
    {-0.230961644455383637, 0.506495587332349035},
    {0.079066118549212831, 2.408163140208925497},
    {0.230961644455383637, -2.635097066257444203},
    {-0.172745327415618701, -1.463445768309359553},
    {-0.605929321571350690, -0.187669323777381622},
    {-0.427370518328979641, 1.252716453253507838},
    {-0.600191595538186799, 2.690988744120037492},
    {-0.491715428198773866, -2.739604450678486295},
    {-0.803582649718989942, -1.893195233972397139},
    {-1.307747883455638156, -0.604647643711872080},
    {-1.054751253523952054, 1.794075294689396615}
};

/**
 * Set i, j, k values in a CoordIJK
 */
void _setIJK(CoordIJK *ijk, int i, int j, int k) {
    ijk->i = i;
    ijk->j = j;
    ijk->k = k;
}

/**
 * Normalize CoordIJK to maintain i + j + k = 0 (Official H3 algorithm)
 */
void _ijkNormalize(CoordIJK *c) {
    // remove any negative values
    if (c->i < 0) {
        c->j -= c->i;
        c->k -= c->i;
        c->i = 0;
    }

    if (c->j < 0) {
        c->i -= c->j;
        c->k -= c->j;
        c->j = 0;
    }

    if (c->k < 0) {
        c->i -= c->k;
        c->j -= c->k;
        c->k = 0;
    }

    // remove the min value if needed
    int min = c->i;
    if (c->j < min) min = c->j;
    if (c->k < min) min = c->k;
    if (min > 0) {
        c->i -= min;
        c->j -= min;
        c->k -= min;
    }
}

/**
 * Add two CoordIJK's
 */
void _ijkAdd(const CoordIJK *a, const CoordIJK *b, CoordIJK *sum) {
    sum->i = a->i + b->i;
    sum->j = a->j + b->j;
    sum->k = a->k + b->k;
}

/**
 * Subtract two CoordIJK's
 */
void _ijkSub(const CoordIJK *a, const CoordIJK *b, CoordIJK *diff) {
    diff->i = a->i - b->i;
    diff->j = a->j - b->j;
    diff->k = a->k - b->k;
}

/**
 * Scale a CoordIJK by a factor
 */
void _ijkScale(CoordIJK *ijk, int factor) {
    ijk->i *= factor;
    ijk->j *= factor;
    ijk->k *= factor;
}

/**
 * Check if resolution is Class III
 */
bool isResolutionClassIII(int res) {
    return res % 2 != 0;
}

/**
 * Aperture operations
 */
void _upAp3(CoordIJK *ijk) {
    _setIJK(ijk, 
            (ijk->i - ijk->k) / 3,
            (ijk->j - ijk->i) / 3,
            (ijk->k - ijk->j) / 3);
}

void _downAp3(CoordIJK *ijk) {
    // res r unit vectors in res r+1
    CoordIJK iVec = {2, 0, 1};
    CoordIJK jVec = {1, 2, 0};
    CoordIJK kVec = {0, 1, 2};
    _ijkScale(&iVec, ijk->i);
    _ijkScale(&jVec, ijk->j);
    _ijkScale(&kVec, ijk->k);
    _ijkAdd(&iVec, &jVec, ijk);
    _ijkAdd(ijk, &kVec, ijk);
    _ijkNormalize(ijk);
}

void _downAp3r(CoordIJK *ijk) {
    // res r unit vectors in res r+1
    CoordIJK iVec = {2, 1, 0};
    CoordIJK jVec = {0, 2, 1};
    CoordIJK kVec = {1, 0, 2};
    _ijkScale(&iVec, ijk->i);
    _ijkScale(&jVec, ijk->j);
    _ijkScale(&kVec, ijk->k);
    _ijkAdd(&iVec, &jVec, ijk);
    _ijkAdd(ijk, &kVec, ijk);
    _ijkNormalize(ijk);
}

void _upAp7(CoordIJK *ijk) {
    int i = ijk->i - ijk->k;
    int j = ijk->j - ijk->k;
    ijk->i = (int)lround((3 * i - j) * M_ONESEVENTH);
    ijk->j = (int)lround((i + 2 * j) * M_ONESEVENTH);
    ijk->k = 0;
    _ijkNormalize(ijk);
}

void _upAp7r(CoordIJK *ijk) {
    int i = ijk->i - ijk->k;
    int j = ijk->j - ijk->k;
    ijk->i = (int)lround((2 * i + j) * M_ONESEVENTH);
    ijk->j = (int)lround((3 * j - i) * M_ONESEVENTH);
    ijk->k = 0;
    _ijkNormalize(ijk);
}

void _downAp7(CoordIJK *ijk) {
    CoordIJK iVec = {3, 0, 1};
    CoordIJK jVec = {1, 3, 0};
    CoordIJK kVec = {0, 1, 3};
    _ijkScale(&iVec, ijk->i);
    _ijkScale(&jVec, ijk->j);
    _ijkScale(&kVec, ijk->k);
    _ijkAdd(&iVec, &jVec, ijk);
    _ijkAdd(ijk, &kVec, ijk);
    _ijkNormalize(ijk);
}

void _downAp7r(CoordIJK *ijk) {
    CoordIJK iVec = {3, 1, 0};
    CoordIJK jVec = {0, 3, 1};
    CoordIJK kVec = {1, 0, 3};
    _ijkScale(&iVec, ijk->i);
    _ijkScale(&jVec, ijk->j);
    _ijkScale(&kVec, ijk->k);
    _ijkAdd(&iVec, &jVec, ijk);
    _ijkAdd(ijk, &kVec, ijk);
    _ijkNormalize(ijk);
}

/**
 * Convert hexagonal coordinates to CoordIJK (Exact H3 implementation)
 */
void _hex2dToCoordIJK(const Vec2d *v, CoordIJK *h) {
    double a1, a2, x1, x2;
    int m1, m2;
    double r1, r2;

    h->k = 0;
    a1 = fabs(v->x);
    a2 = fabs(v->y);

    x2 = a2 * M_RSIN60;
    x1 = a1 + x2 / 2.0;

    m1 = (int)x1;
    m2 = (int)x2;

    r1 = x1 - m1;
    r2 = x2 - m2;

    if (r1 < 0.5) {
        if (r1 < 1.0 / 3.0) {
            if (r2 < (1.0 + r1) / 2.0) {
                h->i = m1;
                h->j = m2;
            } else {
                h->i = m1;
                h->j = m2 + 1;
            }
        } else {
            if (r2 < (1.0 - r1)) {
                h->j = m2;
            } else {
                h->j = m2 + 1;
            }
            if ((1.0 - r1) <= r2 && r2 < (2.0 * r1)) {
                h->i = m1 + 1;
            } else {
                h->i = m1;
            }
        }
    } else {
        if (r1 < 2.0 / 3.0) {
            if (r2 < (1.0 - r1)) {
                h->j = m2;
            } else {
                h->j = m2 + 1;
            }
            if ((2.0 * r1 - 1.0) < r2 && r2 < (1.0 - r1)) {
                h->i = m1;
            } else {
                h->i = m1 + 1;
            }
        } else {
            if (r2 < (r1 / 2.0)) {
                h->i = m1 + 1;
                h->j = m2;
            } else {
                h->i = m1 + 1;
                h->j = m2 + 1;
            }
        }
    }

    // Fold across axes if necessary
    if (v->x < 0.0) {
        if ((h->j % 2) == 0) {
            long long int axisi = h->j / 2;
            long long int diff = h->i - axisi;
            h->i = (int)(h->i - 2.0 * diff);
        } else {
            long long int axisi = (h->j + 1) / 2;
            long long int diff = h->i - axisi;
            h->i = (int)(h->i - (2.0 * diff + 1));
        }
    }

    if (v->y < 0.0) {
        h->i = h->i - (2 * h->j + 1) / 2;
        h->j = -1 * h->j;
    }

    _ijkNormalize(h);
}

static Direction _h3LeadingNonZeroDigit(uint64_t h) {
    int res = H3_GET_RESOLUTION(h);
    for (int r = 1; r <= res; r++) {
        if (H3_GET_INDEX_DIGIT(h, r)) {
            return H3_GET_INDEX_DIGIT(h, r);
        }
    }
    return CENTER_DIGIT;
}

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

static uint64_t _h3Rotate60ccw(uint64_t h) {
    int res = H3_GET_RESOLUTION(h);
    for (int r = 1; r <= res; r++) {
        H3_SET_INDEX_DIGIT(h, r, _rotate60ccw(H3_GET_INDEX_DIGIT(h, r)));
    }
    return h;
}

/**
 * Encodes a coordinate on the sphere to the corresponding icosahedral face and
 * containing 2D hex coordinates relative to that face center.
 * 
 * Copied exactly from H3 reference lib/faceijk.c
 */
void _geoToHex2d(const LatLng *g, int res, int *face, Vec2d *v) {
    // determine the icosahedron face
    double sqd;
    _geoToClosestFace(g, face, &sqd);

    // cos(r) = 1 - 2 * sin^2(r/2) = 1 - 2 * (sqd / 4) = 1 - sqd/2
    double r = acos(1 - sqd * 0.5);

    if (r < EPSILON) {
        v->x = v->y = 0.0;
        return;
    }

    // now have face and r, now find CCW theta from CII i-axis
    double theta =
        _posAngleRads(faceAxesAzRadsCII[*face][0] -
                      _posAngleRads(_geoAzimuthRads(&faceCenterGeo[*face], g)));

    // adjust theta for Class III (odd resolutions)
    if (isResolutionClassIII(res))
        theta = _posAngleRads(theta - M_AP7_ROT_RADS);

    // perform gnomonic scaling of r
    r = tan(r);

    // scale for current resolution length u
    r *= INV_RES0_U_GNOMONIC;
    for (int i = 0; i < res; i++) r *= M_SQRT7;

    // we now have (r, theta) in hex2d with theta ccw from x-axes

    // convert to local x,y
    v->x = r * cos(theta);
    v->y = r * sin(theta);
}

/**
 * Encodes a coordinate on the sphere to the corresponding icosahedral face and
 * the squared euclidean distance to that face center.
 */
void _geoToClosestFace(const LatLng *g, int *face, double *sqd) {
    Vec3d v3d;
    _geoToVec3d(g, &v3d);

    // determine the icosahedron face
    *face = 0;
    // The distance between two farthest points is 2.0, therefore the square of
    // the distance between two points should always be less or equal than 4.0 .
    *sqd = 5.0;
    for (int f = 0; f < NUM_ICOSA_FACES; ++f) {
        double sqdT = _pointSquareDist(&faceCenterPoint[f], &v3d);
        if (sqdT < *sqd) {
            *face = f;
            *sqd = sqdT;
        }
    }
}

/**
 * Adjusts a FaceIJK address in place so that the resulting cell address is
 * relative to the correct icosahedral face.
 */
Overage _adjustOverageClassII(FaceIJK *fijk, int res, int pentLeading4,
                              int substrate) {
    Overage overage = NO_OVERAGE;

    CoordIJK *ijk = &fijk->coord;

    // get the maximum dimension value; scale if a substrate grid
    int maxDim = maxDimByCIIres[res];
    if (substrate) maxDim *= 3;

    // check for overage
    if (substrate && ijk->i + ijk->j + ijk->k == maxDim)  // on edge
        overage = FACE_EDGE;
    else if (ijk->i + ijk->j + ijk->k > maxDim)  // overage
    {
        overage = NEW_FACE;

        const FaceOrientIJK *fijkOrient;
        if (ijk->k > 0) {
            if (ijk->j > 0)  // jk "quadrant"
                fijkOrient = &faceNeighbors[fijk->face][JK];
            else  // ik "quadrant"
            {
                fijkOrient = &faceNeighbors[fijk->face][KI];

                // adjust for the pentagonal missing sequence
                if (pentLeading4) {
                    // translate origin to center of pentagon
                    CoordIJK origin;
                    _setIJK(&origin, maxDim, 0, 0);
                    CoordIJK tmp;
                    _ijkSub(ijk, &origin, &tmp);
                    // rotate to adjust for the missing sequence
                    _ijkRotate60cw(&tmp);
                    // translate the origin back to the center of the triangle
                    _ijkAdd(&tmp, &origin, ijk);
                }
            }
        } else  // ij "quadrant"
            fijkOrient = &faceNeighbors[fijk->face][IJ];

        fijk->face = fijkOrient->face;

        // rotate and translate for adjacent face
        for (int i = 0; i < fijkOrient->ccwRot60; i++) _ijkRotate60ccw(ijk);

        CoordIJK transVec = fijkOrient->translate;
        int unitScale = unitScaleByCIIres[res];
        if (substrate) unitScale *= 3;
        _ijkScale(&transVec, unitScale);
        _ijkAdd(ijk, &transVec, ijk);
        _ijkNormalize(ijk);

        // overage points on pentagon boundaries can end up on edges
        if (substrate && ijk->i + ijk->j + ijk->k == maxDim)  // on edge
            overage = FACE_EDGE;
    }

    return overage;
}

static uint64_t _h3RotatePent60ccw(uint64_t h) {
    int foundFirstNonZeroDigit = 0;
    for (int r = 1, res = H3_GET_RESOLUTION(h); r <= res; r++) {
        H3_SET_INDEX_DIGIT(h, r, _rotate60ccw(H3_GET_INDEX_DIGIT(h, r)));
        if (!foundFirstNonZeroDigit && H3_GET_INDEX_DIGIT(h, r) != 0) {
            foundFirstNonZeroDigit = 1;
            if (_h3LeadingNonZeroDigit(h) == K_AXES_DIGIT)
                h = _h3Rotate60ccw(h);
        }
    }
    return h;
}

static uint64_t _h3Rotate60cw(uint64_t h) {
    int res = H3_GET_RESOLUTION(h);
    for (int r = 1; r <= res; r++) {
        H3_SET_INDEX_DIGIT(h, r, _rotate60cw(H3_GET_INDEX_DIGIT(h, r)));
    }
    return h;
}

static Direction _unitIjkToDigit(const CoordIJK *ijk) {
    CoordIJK c = *ijk;
    _ijkNormalize(&c);
    
    // Check against H3 UNIT_VECS
    if (c.i == 0 && c.j == 0 && c.k == 0) return CENTER_DIGIT;  // {0,0,0}
    if (c.i == 0 && c.j == 0 && c.k == 1) return K_AXES_DIGIT;  // {0,0,1}
    if (c.i == 0 && c.j == 1 && c.k == 0) return J_AXES_DIGIT;  // {0,1,0}
    if (c.i == 0 && c.j == 1 && c.k == 1) return JK_AXES_DIGIT; // {0,1,1}
    if (c.i == 1 && c.j == 0 && c.k == 0) return I_AXES_DIGIT;  // {1,0,0}
    if (c.i == 1 && c.j == 0 && c.k == 1) return IK_AXES_DIGIT; // {1,0,1}
    if (c.i == 1 && c.j == 1 && c.k == 0) return IJ_AXES_DIGIT; // {1,1,0}
    return INVALID_DIGIT;
}

void _geoToVec3d(const LatLng *geo, Vec3d *vec) {
    vec->x = cos(geo->lat) * cos(geo->lng);
    vec->y = cos(geo->lat) * sin(geo->lng);
    vec->z = sin(geo->lat);
}

double _pointSquareDist(const Vec3d *p1, const Vec3d *p2) {
    return (p1->x - p2->x) * (p1->x - p2->x) +
           (p1->y - p2->y) * (p1->y - p2->y) +
           (p1->z - p2->z) * (p1->z - p2->z);
}

double _posAngleRads(double rads) {
    return rads < 0 ? (rads + 2 * M_PI) : rads;
}

double _geoAzimuthRads(const LatLng *p1, const LatLng *p2) {
    return atan2(
        cos(p2->lat) * sin(p2->lng - p1->lng),
        cos(p1->lat) * sin(p2->lat) -
            sin(p1->lat) * cos(p2->lat) * cos(p2->lng - p1->lng));
}

void _ijkRotate60ccw(CoordIJK *ijk) {
    CoordIJK iVec = {1, 1, 0};
    CoordIJK jVec = {0, 1, 1};
    CoordIJK kVec = {1, 0, 1};
    _ijkScale(&iVec, ijk->i);
    _ijkScale(&jVec, ijk->j);
    _ijkScale(&kVec, ijk->k);
    _ijkAdd(&iVec, &jVec, ijk);
    _ijkAdd(ijk, &kVec, ijk);
    _ijkNormalize(ijk);
}

void _ijkRotate60cw(CoordIJK *ijk) {
    CoordIJK iVec = {1, 0, 1};
    CoordIJK jVec = {1, 1, 0};
    CoordIJK kVec = {0, 1, 1};
    _ijkScale(&iVec, ijk->i);
    _ijkScale(&jVec, ijk->j);
    _ijkScale(&kVec, ijk->k);
    _ijkAdd(&iVec, &jVec, ijk);
    _ijkAdd(ijk, &kVec, ijk);
    _ijkNormalize(ijk);
}

/**
 * Convert IJK coordinates to 2D hex coordinates
 */
void _ijkToHex2d(const CoordIJK *ijk, Vec2d *h) {
    int i = ijk->i - ijk->k;
    int j = ijk->j - ijk->k;
    h->x = i - 0.5 * j;
    h->y = j * M_SQRT3_2;
}

/**
 * Get magnitude of 2D vector
 */
double _v2dMag(const Vec2d *v) {
    return sqrt(v->x * v->x + v->y * v->y);
}

/**
 * Check if two 2D vectors are almost equal
 */
bool _v2dAlmostEquals(const Vec2d *v1, const Vec2d *v2) {
    return fabs(v1->x - v2->x) < EPSILON && fabs(v1->y - v2->y) < EPSILON;
}

/**
 * Find intersection of two line segments
 */
void _v2dIntersect(const Vec2d *p1, const Vec2d *p2, const Vec2d *q1, const Vec2d *q2, Vec2d *r) {
    double s1_x = p2->x - p1->x;
    double s1_y = p2->y - p1->y;
    double s2_x = q2->x - q1->x;
    double s2_y = q2->y - q1->y;

    double t = (s2_x * (p1->y - q1->y) - s2_y * (p1->x - q1->x)) /
               (s2_x * s1_y - s2_y * s1_x);

    r->x = p1->x + (t * s1_x);
    r->y = p1->y + (t * s1_y);
}

/**
 * Computes the point on the sphere at azimuth and distance from another point
 */
void _geoAzDistanceRads(const LatLng *p1, double az, double distance, LatLng *p2) {
    if (distance < EPSILON) {
        *p2 = *p1;
        return;
    }

    double sinlat = sin(p1->lat) * cos(distance) +
                    cos(p1->lat) * sin(distance) * cos(az);
    if (sinlat > 1.0) sinlat = 1.0;
    if (sinlat < -1.0) sinlat = -1.0;
    p2->lat = asin(sinlat);
    if (fabs(p2->lat - M_PI_2) < EPSILON) {
        p2->lng = 0.0;
    } else if (fabs(p2->lat + M_PI_2) < EPSILON) {
        p2->lng = 0.0;
    } else {
        double sinlng = sin(az) * sin(distance) / cos(p2->lat);
        double coslng = (cos(distance) - sin(p1->lat) * sin(p2->lat)) /
                        (cos(p1->lat) * cos(p2->lat));
        if (sinlng > 1.0) sinlng = 1.0;
        if (sinlng < -1.0) sinlng = -1.0;
        if (coslng > 1.0) coslng = 1.0;
        if (coslng < -1.0) coslng = -1.0;
        p2->lng = p1->lng + atan2(sinlng, coslng);
    }
}

/**
 * Encodes a coordinate on the sphere to the FaceIJK address
 */
void _geoToFaceIjk(const LatLng *g, int res, FaceIJK *h) {
    // first convert to hex2d
    Vec2d v;
    _geoToHex2d(g, res, &h->face, &v);

    // then convert to ijk+
    _hex2dToCoordIJK(&v, &h->coord);
}

/**
 * Public wrapper for geoToFaceIjk
 * Takes lat/lng in degrees and converts to radians
 */
void geoToFaceIjk(double lat, double lng, int res, FaceIJK *fijk) {
    LatLng g = {lat * M_PI_180, lng * M_PI_180};
    _geoToFaceIjk(&g, res, fijk);
}

/**
 * Convert FaceIJK to H3 index
 * This is adapted from H3 reference implementation _faceIjkToH3
 */
uint64_t faceIjkToH3(const FaceIJK *fijk, int res) {
    // initialize the index
    uint64_t h = H3_INIT;
    H3_SET_MODE(h, H3_CELL_MODE);
    H3_SET_RESOLUTION(h, res);

    // check for res 0/base cell
    if (res == 0) {
        if (fijk->coord.i > MAX_FACE_COORD || fijk->coord.j > MAX_FACE_COORD ||
            fijk->coord.k > MAX_FACE_COORD) {
            // out of range input
            return 0;
        }

        H3_SET_BASE_CELL(h, _faceIjkToBaseCell(fijk));
        return h;
    }

    // we need to find the correct base cell FaceIJK for this H3 index;
    // start with the passed in face and resolution res ijk coordinates
    // in that face's coordinate system
    FaceIJK fijkBC = *fijk;

    // build the H3Index from finest res up
    // adjust r for the fact that the res 0 base cell offsets the indexing
    // digits
    CoordIJK *ijk = &fijkBC.coord;
    for (int r = res - 1; r >= 0; r--) {
        CoordIJK lastIJK = *ijk;
        CoordIJK lastCenter;
        if (isResolutionClassIII(r + 1)) {
            // rotate ccw
            _upAp7(ijk);
            lastCenter = *ijk;
            _downAp7(&lastCenter);
        } else {
            // rotate cw
            _upAp7r(ijk);
            lastCenter = *ijk;
            _downAp7r(&lastCenter);
        }

        CoordIJK diff;
        _ijkSub(&lastIJK, &lastCenter, &diff);
        _ijkNormalize(&diff);

        H3_SET_INDEX_DIGIT(h, r + 1, _unitIjkToDigit(&diff));
    }

    // fijkBC should now hold the IJK of the base cell in the
    // coordinate system of the current face

    if (fijkBC.coord.i > MAX_FACE_COORD || fijkBC.coord.j > MAX_FACE_COORD ||
        fijkBC.coord.k > MAX_FACE_COORD) {
        // out of range input
        return 0;
    }

    // lookup the correct base cell
    int baseCell = _faceIjkToBaseCell(&fijkBC);
    H3_SET_BASE_CELL(h, baseCell);

    // rotate if necessary to get canonical base cell orientation
    // for this base cell
    int numRots = _faceIjkToBaseCellCCWrot60(&fijkBC);
    if (_isBaseCellPentagon(baseCell)) {
        // force rotation out of missing k-axes sub-sequence
        if (_h3LeadingNonZeroDigit(h) == K_AXES_DIGIT) {
            // check for a cw/ccw offset face; default is ccw
            if (_baseCellIsCwOffset(baseCell, fijkBC.face)) {
                h = _h3Rotate60cw(h);
            } else {
                h = _h3Rotate60ccw(h);
            }
        }

        for (int i = 0; i < numRots; i++) h = _h3RotatePent60ccw(h);
    } else {
        for (int i = 0; i < numRots; i++) {
            h = _h3Rotate60ccw(h);
        }
    }

    return h;
}
