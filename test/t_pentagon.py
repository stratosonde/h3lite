#!/usr/bin/env python3
"""
t_pentagon.py — pentagon ring differential test (review F-01 regression)

For every res-0 pentagon and a sample of res-3 cells within ~3 rings of
each pentagon, compare h3lite's h3GetRing(k=0..6) against the reference
h3 library's grid_ring as SETS.

Classification (per handoff §5.1):
  wrong   : C returns cells but the set differs from reference — FATAL
  clean   : C returns FAIL while reference has a ring — allowed ONLY for
            k=1 with a pentagon origin (legitimate H3 pentagon distortion)
  exact   : sets match

Exit 0 iff wrong == 0 and every clean failure is the allowed k=1 pentagon
case.

Requires: gcc, the h3 python package, and `make lib` artifacts.
Run from the h3lite repo root:  python3 test/t_pentagon.py
"""
import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "xval_ring")

try:
    import h3
except ImportError:
    print("SKIP: h3 python package not installed")
    sys.exit(0)


def build():
    if os.path.exists(BIN):
        return True
    gcc = shutil.which("gcc")
    if not gcc:
        return False
    lib = os.path.join(ROOT, "bin", "libh3lite.a")
    if not os.path.exists(lib):
        r = subprocess.run(["make", "lib"], cwd=ROOT)
        if r.returncode != 0:
            return False
    r = subprocess.run(
        [gcc, "-Wall", "-O2", "-I./include", "test/xval_ring.c",
         "-L./bin", "-lh3lite", "-lm", "-o", "xval_ring"], cwd=ROOT)
    return r.returncode == 0


def main():
    if not build():
        print("SKIP: cannot build xval_ring (gcc or libh3lite.a missing)")
        sys.exit(0)

    # Origins: the 12 res-0 pentagons, their immediate res-3 children,
    # and a sample of cells at distance 2-3 (where the F-01 predicate bug
    # silently killed every ring search).
    origins = set()
    pentagons0 = sorted(h3.get_pentagons(0))
    for p in pentagons0:
        origins.add(p)
        kids = sorted(h3.cell_to_children(p, 3))
        origins.update(kids[::7])           # sample of pentagon descendants
        for ring_cell in h3.grid_ring(p, 2):
            origins.add(ring_cell)
            origins.update(h3.grid_ring(ring_cell, 1)[::3])
    origins = sorted(origins)
    print(f"Testing {len(origins)} origins around {len(pentagons0)} "
          f"pentagons, k=0..6")

    # Batch all requests through one xval_ring process
    requests = [(o, k) for o in origins for k in range(0, 7)]
    stdin = "".join(
        f"{int(o, 16) >> 32:08x}{int(o, 16) & 0xFFFFFFFF:08x} {k}\n"
        for o, k in requests)
    proc = subprocess.run([BIN], input=stdin, capture_output=True,
                          text=True, cwd=ROOT)
    lines = proc.stderr.strip().splitlines()
    if len(lines) != len(requests):
        print(f"FAIL: expected {len(requests)} results, got {len(lines)}")
        sys.exit(1)

    wrong = clean = exact = 0
    bad_examples = []
    pentagon_set = set(pentagons0)
    for (origin, k), line in zip(requests, lines):
        ref = set(h3.grid_ring(origin, k))
        if line.strip() == "FAIL":
            if k == 1 and origin in pentagon_set:
                clean += 1  # legitimate pentagon distortion failure
            else:
                wrong += 1  # clean-looking but should have produced cells
                bad_examples.append((origin, k, "FAIL-vs-%d" % len(ref)))
            continue
        parts = line.split()
        got = set(p.lower().lstrip("0") or "0" for p in parts[1:])
        ref_n = set(r.lower().lstrip("0") or "0" for r in ref)
        if got == ref_n:
            exact += 1
        else:
            wrong += 1
            if len(bad_examples) < 5:
                bad_examples.append((origin, k,
                                     f"got {len(got)} want {len(ref_n)}"))

    print(f"exact={exact}  clean-fail(k=1 pentagon)={clean}  wrong={wrong}")
    for o, k, why in bad_examples:
        print(f"  WRONG: origin={o} k={k} {why}")

    if wrong:
        print("RESULT: FAIL")
        sys.exit(1)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
