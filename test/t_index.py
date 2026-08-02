"""
T-INDEX — global grid cross-validation of the H3Lite geometry engine.

Compares h3lite's full 64-bit index at resolution 3 against the reference
`h3` Python package (validated against h3==4.5.0) on a global 0.5-degree
grid. Baseline (handoff Appendix A): 0 mismatches. Must stay 0.

Usage:  python test/t_index.py
Requires: pip install h3==4.5.0
"""
import sys

import harness_common as hc


def grid_points(step=0.5):
    """Global grid, cell centers: lat -89.75..89.75, lon -179.75..179.75."""
    n_lat = int(180.0 / step)   # 360
    n_lon = int(360.0 / step)   # 720
    pts = []
    for i in range(n_lat):
        lat = -90.0 + step / 2 + i * step
        for j in range(n_lon):
            lon = -180.0 + step / 2 + j * step
            pts.append((lat, lon))
    return pts


def main():
    try:
        import h3
    except ImportError:
        sys.exit("error: pip install h3==4.5.0 is required for T-INDEX")

    pts = grid_points()
    print(f"T-INDEX: {len(pts)} grid points at res 3", file=sys.stderr)

    rows = hc.run_points(pts)

    mismatches = 0
    first = []
    for (lat, lon), row in zip(pts, rows):
        ref = h3.latlng_to_cell(lat, lon, 3)
        if int(row["h3"], 16) != int(ref, 16):
            mismatches += 1
            if len(first) < 10:
                first.append((lat, lon, row["h3"], ref))

    print(f"T-INDEX: {mismatches} / {len(pts)} mismatches")
    for lat, lon, got, want in first:
        print(f"  {lat:+.2f} {lon:+.2f}  h3lite={got}  ref={want}")
    sys.exit(0 if mismatches == 0 else 1)


if __name__ == "__main__":
    main()
