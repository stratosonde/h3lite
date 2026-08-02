"""
T-TABLE — static audit of the generated region lookup table (handoff Appendix C).

Run against src/h3lite_regions_table.c on every regeneration.
Must report zero duplicates and zero conflicts.

Supports both table formats:
  * legacy struct entries  { baseCell, partialIndex, regionId }   (pre-H3-1)
    — resolution is not carried; for the duplicate/conflict audit the key is
      (baseCell, partialIndex). Baseline: 269 duplicates / 35 conflicting.
  * packed uint32 entries  (post-H3-1, bits [31:25] bc, [24:23] res,
      [22:14] partialIndex, [13:10] regionId). Required: 0 duplicates,
      0 conflicts, table sorted by (baseCell, resolution, partialIndex).

Usage:  python test/t_table.py
"""
import os
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
TABLE = os.path.join(os.path.dirname(HERE), "src", "h3lite_regions_table.c")


def decode_packed(e):
    bc = (e >> 25) & 0x7F
    res = (e >> 23) & 0x03
    pi = (e >> 14) & 0x1FF
    region = (e >> 10) & 0x0F
    return bc, res, pi, region


def main():
    txt = open(TABLE, encoding="utf-8").read()
    m = re.search(r'regionLookup\[[^\]]*\]\s*=\s*\{(.*?)\n\};', txt, re.S)
    if not m:
        sys.exit("error: could not locate regionLookup[] body")
    body = m.group(1)

    legacy = re.findall(r'\{\s*(\d+),\s*(\d+),\s*(\d+)\s*\}', body)
    if legacy:
        fmt = "legacy struct"
        # (bc, res=3 implied, pi, region) — res unknown, use key (bc, pi)
        ents = [(int(bc), None, int(pi), int(reg)) for bc, pi, reg in legacy]
        key = lambda e: (e[0], e[2])
    else:
        fmt = "packed uint32"
        vals = [int(v, 0) for v in re.findall(r'0[xX][0-9a-fA-F]+|\b\d{5,}\b', body)]
        if not vals:
            sys.exit("error: no entries found in regionLookup body")
        ents = [decode_packed(v) for v in vals]
        key = lambda e: (e[0], e[1], e[2])

    print(f"T-TABLE: format = {fmt}, entries = {len(ents)}")

    problems = 0

    # sorted-order assertion (packed format must be binary-searchable)
    if fmt == "packed uint32":
        if ents != sorted(ents, key=lambda e: (e[0], e[1], e[2])):
            print("  FAIL: table not sorted by (baseCell, resolution, partialIndex)")
            problems += 1

    keys = Counter(key(e) for e in ents)
    dups = {k: v for k, v in keys.items() if v > 1}
    print(f"  duplicate keys:   {len(dups)}")
    if dups:
        problems += 1

    regions_by_key = defaultdict(set)
    for e in ents:
        regions_by_key[key(e)].add(e[3])
    conflicts = {k: v for k, v in regions_by_key.items() if len(v) > 1}
    print(f"  conflicting keys: {len(conflicts)}")
    if conflicts:
        problems += 1

    # range sanity (meaningful for both formats)
    if max(e[2] for e in ents) >= 512:
        print("  FAIL: partialIndex overflows 9 bits")
        problems += 1
    if max(e[0] for e in ents) > 121:
        print("  FAIL: baseCell out of range")
        problems += 1
    if fmt == "packed uint32" and not all(1 <= e[1] <= 3 for e in ents):
        print("  FAIL: resolution out of range")
        problems += 1

    if problems:
        print(f"T-TABLE: FAIL ({problems} problem classes)")
        sys.exit(1)
    print(f"T-TABLE: OK ({len(ents)} entries, no duplicates, no conflicts)")


if __name__ == "__main__":
    main()
