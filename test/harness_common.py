"""
harness_common.py — shared helpers for the H3Lite cross-validation harness.

Locates and (if needed) builds the xval_pts binary, feeds it lat/lon pairs,
and parses its stderr output lines:
    lat lon baseCell partialIndex region h3hex
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)  # h3lite repo root

EXE = os.path.join(ROOT, "xval_pts.exe" if os.name == "nt" else "xval_pts")

SRCS = [os.path.join(ROOT, "src", f) for f in (
    "h3lite.c", "h3lite_faceijk.c", "h3lite_basecells.c",
    "h3lite_neighbor.c", "h3lite_regions_table.c")]

GCC_CANDIDATES = ["gcc", r"C:\MinGW\bin\gcc.exe"]


def _find_gcc():
    for cand in GCC_CANDIDATES:
        try:
            subprocess.run([cand, "--version"], capture_output=True, check=True)
            return cand
        except (OSError, subprocess.CalledProcessError):
            continue
    return None


def build():
    """Build xval_pts if it is missing or older than the sources."""
    if os.path.exists(EXE):
        exe_mt = os.path.getmtime(EXE)
        srcs = SRCS + [os.path.join(HERE, "xval_pts.c")]
        if all(os.path.getmtime(s) < exe_mt for s in srcs if os.path.exists(s)):
            return EXE
    gcc = _find_gcc()
    if gcc is None:
        sys.exit("error: no host gcc found (looked for gcc, C:\\MinGW\\bin\\gcc.exe)")
    cmd = [gcc, "-Wall", "-O2",
           "-I", os.path.join(ROOT, "include"),
           os.path.join(HERE, "xval_pts.c")] + SRCS + \
          ["-lm", "-o", EXE]
    print("build:", " ".join(cmd), file=sys.stderr)
    # cc1 resolves its support DLLs (libgmp/mpfr/mpc) via PATH; ensure the
    # compiler's own bin directory is present (MinGW is not on the user PATH).
    env = dict(os.environ)
    gcc_dir = os.path.dirname(os.path.abspath(gcc))
    if gcc_dir:
        env["PATH"] = gcc_dir + os.pathsep + env.get("PATH", "")
    subprocess.run(cmd, check=True, env=env)
    return EXE


def run_points(points):
    """
    Feed an iterable of (lat, lon) floats to xval_pts.
    Returns a list of dicts in the same order:
        {lat, lon, base_cell, partial_index, region, h3 (hex string)}
    """
    exe = build()
    data = "".join(f"{lat:.6f} {lon:.6f}\n" for lat, lon in points)
    proc = subprocess.run([exe], input=data, capture_output=True,
                          text=True, check=True)
    out = []
    for line in proc.stderr.splitlines():
        parts = line.split()
        if len(parts) != 6:
            continue
        out.append({
            "lat": float(parts[0]),
            "lon": float(parts[1]),
            "base_cell": int(parts[2]),
            "partial_index": int(parts[3]),
            "region": int(parts[4]),
            "h3": parts[5],
        })
    if len(out) != len(list(points)):
        raise RuntimeError(f"xval_pts emitted {len(out)} rows for "
                           f"{len(list(points))} input points")
    return out


def region_names():
    """Parse regionNames[] from the generated table source: {id: name}."""
    import re
    path = os.path.join(ROOT, "src", "h3lite_regions_table.c")
    txt = open(path, encoding="utf-8").read()
    m = re.search(r'regionNames\[\d+\]\s*=\s*\{(.*?)\};', txt, re.S)
    if not m:
        return {0: "Unknown"}
    names = re.findall(r'"([^"]*)"', m.group(1))
    return {i: n for i, n in enumerate(names)}
