/*
 * xval_ring.c — ring cross-validation harness (review F-01 regression)
 *
 * Reads "h3hex k" pairs on stdin, emits one line per request on stderr:
 *   <count> <cellhex> <cellhex> ...      on success
 *   FAIL                                 when h3GetRing reports failure
 *
 * Build (from h3lite repo root, after `make lib`):
 *   gcc -Wall -O2 -I./include test/xval_ring.c -L./bin -lh3lite -lm -o xval_ring
 */
#include <stdio.h>
#include "h3lite.h"

#define MAX_RING 42 /* k <= 6 -> 36 cells + margin, matches h3lite.c */

int main(void) {
    unsigned long hi, lo;
    int k;
    h3liteInit();
    while (scanf("%8lx%8lx %d", &hi, &lo, &k) == 3) {
        H3Index origin = ((H3Index)hi << 32) | (H3Index)lo;
        H3Index ring[MAX_RING];
        int n = h3GetRing(origin, k, ring);
        if (n < 0) {
            fprintf(stderr, "FAIL\n");
            continue;
        }
        fprintf(stderr, "%d", n);
        for (int i = 0; i < n; i++) {
            fprintf(stderr, " %08lx%08lx",
                    (unsigned long)(ring[i] >> 32),
                    (unsigned long)(ring[i] & 0xFFFFFFFFu));
        }
        fprintf(stderr, "\n");
    }
    return 0;
}
