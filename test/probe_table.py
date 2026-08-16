#!/usr/bin/env python3
"""Probe the generated packed region table exactly as the device does.

Replicates: latLngToH3 (via reference h3) -> h3ToRegion probing
res 3 -> res 2 -> res 1 with the packed (baseCell, res, partialIndex) key.
Also checks which hplans polygons contain each probe point (F-04 analysis).
"""
import re, sys, json, os
import h3
from shapely.geometry import shape, Point

H3LITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HPLANS = r"c:\working\sonde\hplans"

# ---- parse generated table ----
entries = {}  # (baseCell, res, partialIndex) -> regionId
names = {}
with open(os.path.join(H3LITE, "src", "h3lite_regions_table.c")) as f:
    text = f.read()
for m in re.finditer(r'"(Unknown|[A-Z0-9\-]+)", // ID (\d+)', text):
    names[int(m.group(2))] = m.group(1)
for m in re.finditer(r'0x([0-9A-F]{8})u, // bc=(\d+) res=(\d+) pi=(\d+) (\S+)', text):
    packed = int(m.group(1), 16)
    bc, res, pi = int(m.group(2)), int(m.group(3)), int(m.group(4))
    rid = (packed >> 10) & 0x0F
    entries[(bc, res, pi)] = rid

def h3_to_region(h):
    h_int = int(h, 16) if isinstance(h, str) else h
    res = (h_int >> 52) & 0xF
    bc = (h_int >> 45) & 0x7F
    d = [0] + [(h_int >> (45 - 3 * r)) & 0x7 for r in (1, 2, 3)]
    for probe_res in (3, 2, 1):
        if res < probe_res:
            continue
        pi = 0
        for r in range(1, probe_res + 1):
            pi = pi * 8 + d[r]
        rid = entries.get((bc, probe_res, pi))
        if rid:
            return rid
    return 0

def lookup(lat, lng):
    cell = h3.latlng_to_cell(lat, lng, 3)
    return h3_to_region(cell)

PROBES = [
    ("Kinshasa, DRC",        -4.32,   15.31),
    ("DRC center",           -2.50,   23.00),
    ("Hong Kong",            22.32,  114.17),
    ("Macau",                22.17,  113.55),
    ("Ulaanbaatar, Mongolia", 47.92, 106.92),
    ("Suva, Fiji",          -17.71,  178.07),
    ("Funafuti, Tuvalu",     -8.52,  179.20),
    ("Tarawa, Kiribati",      1.45,  173.03),
    ("Majuro, Marshall Is.",  7.09,  171.38),
    ("Palikir, FSM",          6.92,  158.16),
    ("Yaren, Nauru",         -0.55,  166.92),
    ("Ngerulmud, Palau",      7.50,  134.62),
    ("Male, Maldives",        4.18,   73.51),
    ("NYC, USA",             40.71,  -74.01),
    ("Paris, France",        48.86,    2.35),
    ("Tokyo, Japan",         35.68,  139.69),
    ("Sydney, Australia",   -33.87,  151.21),
    ("New Delhi, India",     28.61,   77.21),
    ("Seoul, Korea",         37.57,  126.98),
    ("Wellington, NZ",      -41.29,  174.78),
    ("Auckland, NZ",        -36.85,  174.76),
    ("Beijing, China",       39.90,  116.40),
    ("Moscow, Russia",       55.76,   37.62),
    ("Sao Paulo, Brazil",   -23.55,  -46.63),
    ("Mid-Atlantic",         30.00,  -40.00),
]

print("== Device-equivalent lookup probes ==")
results = {}
for name, lat, lng in PROBES:
    rid = lookup(lat, lng)
    results[name] = names.get(rid, f"?{rid}")
    print(f"  {name:26s} -> {results[name]}")

# ---- F-04: which hplans polygons contain the probe points? ----
print("\n== hplans polygon containment (F-04 data-gap analysis) ==")
geoms = {}
for fn in sorted(os.listdir(HPLANS)):
    if not fn.endswith(".geojson"):
        continue
    with open(os.path.join(HPLANS, fn), encoding="utf-8") as f:
        gj = json.load(f)
    gs = []
    if gj.get("type") == "FeatureCollection":
        gs = [shape(ft["geometry"]) for ft in gj.get("features", [])]
    elif gj.get("type") == "Feature":
        gs = [shape(gj["geometry"])]
    else:
        gs = [shape(gj)]
    geoms[os.path.splitext(fn)[0]] = gs

for name, lat, lng in PROBES:
    pt = Point(lng, lat)
    hits = [rn for rn, gs in geoms.items() if any(g.contains(pt) for g in gs)]
    print(f"  {name:26s} polygons: {hits if hits else 'NONE'}")
