/*
 * xval_pts.c — H3Lite cross-validation harness (handoff §1.1)
 *
 * Reads "lat lon" pairs on stdin, emits one line per point on stderr:
 *   lat lon baseCell partialIndex region h3hex
 *
 * Build (from h3lite repo root, after `make lib`):
 *   gcc -Wall -O2 -I./include test/xval_pts.c -L./bin -lh3lite -lm -o xval_pts
 */
#include <stdio.h>
#include "h3lite.h"
#include "h3lite_constants.h"

int main(void) {
    double lat, lon;
    h3liteInit();
    while (scanf("%lf %lf", &lat, &lon) == 2) {
        H3Index h = latLngToH3(lat, lon, 3);
        int bc = -1;
        unsigned pi = 0;
        if (h) {
            bc = H3_GET_BASE_CELL(h);
            int res = H3_GET_RESOLUTION(h);
            int nd = res < 3 ? res : 3;
            for (int r = 1; r <= nd; r++)
                pi = pi * 8 + H3_GET_INDEX_DIGIT(h, r);
        }
        RegionId reg = latLngToRegion(lat, lon);
        /* Print the 64-bit index as two 32-bit halves: %llx is not portable
         * (MSVCRT/MinGW.org lacks it), and zero-padding keeps the string
         * comparison in t_index.py trivial. */
        fprintf(stderr, "%.6f %.6f %d %u %d %08lx%08lx\n",
                lat, lon, bc, pi, (int)reg,
                (unsigned long)(h >> 32), (unsigned long)(h & 0xFFFFFFFFu));
    }
    return 0;
}
