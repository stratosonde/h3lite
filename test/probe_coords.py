#!/usr/bin/env python3
"""Device-equivalent region probe for arbitrary coordinates.

Usage:  python test/probe_coords.py '<label>' <lat> <lng> '<label>' <lat> <lng> ...
Replicates latLngToH3 -> h3ToRegion probing res 3 -> 2 -> 1 and prints the
region id plus name. Exit code 0 on success even for Unknown (0).
"""
import re, sys, os
import h3

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

names = {}
entries = {}
with open(os.path.join(ROOT, "src", "h3lite_regions_table.c")) as f:
    text = f.read()
for m in re.finditer(r'"(Unknown|[A-Z0-9\-]+)", // ID (\d+)', text):
    names[int(m.group(2))] = m.group(1)
for m in re.finditer(r'0x([0-9A-F]{8})u, // bc=(\d+) res=(\d+) pi=(\d+) (\S+)',
                     text):
    packed = int(m.group(1), 16)
    bc, res, pi = int(m.group(2)), int(m.group(3)), int(m.group(4))
    entries[(bc, res, pi)] = (packed >> 10) & 0x0F


def lookup(lat, lng):
    h_int = int(h3.latlng_to_cell(lat, lng, 3), 16)
    base_res = (h_int >> 52) & 0xF
    bc = (h_int >> 45) & 0x7F
    digits = [0] + [(h_int >> (45 - 3 * r)) & 0x7 for r in (1, 2, 3)]
    for probe_res in (3, 2, 1):
        if base_res < probe_res:
            continue
        pi = 0
        for r in range(1, probe_res + 1):
            pi = pi * 8 + digits[r]
        rid = entries.get((bc, probe_res, pi))
        if rid:
            return rid
    return 0


def main():
    args = sys.argv[1:]
    if len(args) < 3 or len(args) % 3 != 0:
        sys.exit("usage: probe_coords.py '<label>' <lat> <lng> ...")
    for i in range(0, len(args), 3):
        label, lat, lng = args[i], float(args[i + 1]), float(args[i + 2])
        rid = lookup(lat, lng)
        print(f"{label} -> {names.get(rid, f'ID={rid}')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
