"""
T-CITY — 52 known-land cities vs expected LoRaWAN region (handoff Appendix B).

Treats any AS923-* result as matching expected AS923. New Zealand and Hong
Kong expectations are PROVISIONAL pending H3-3 (RP002 verification).

Baseline (handoff Appendix A): 36 correct / 15 Unknown / 1 wrong.
Required after Phase 2: >= 98% correct (>= 51/52), 0 wrong.

Usage:  python test/t_city.py
"""
import sys

import harness_common as hc

# lat, lon, expected_region, name   (NZ/HK provisional per H3-3)
CITIES = [
    (40.71, -74.01, "US915", "NewYork"),      (41.88, -87.63, "US915", "Chicago"),
    (39.74, -104.99, "US915", "Denver"),      (29.76, -95.37, "US915", "Houston"),
    (51.05, -114.07, "US915", "Calgary"),     (45.42, -75.70, "US915", "Ottawa"),
    (49.28, -123.12, "US915", "Vancouver"),   (19.43, -99.13, "US915", "MexicoCity"),
    (51.51, -0.13, "EU868", "London"),        (48.86, 2.35, "EU868", "Paris"),
    (52.52, 13.40, "EU868", "Berlin"),        (41.90, 12.50, "EU868", "Rome"),
    (40.42, -3.70, "EU868", "Madrid"),        (59.33, 18.07, "EU868", "Stockholm"),
    (64.13, -21.90, "EU868", "Reykjavik"),    (37.98, 23.73, "EU868", "Athens"),
    (53.35, -6.26, "EU868", "Dublin"),        (-33.87, 151.21, "AU915", "Sydney"),
    (-37.81, 144.96, "AU915", "Melbourne"),   (-31.95, 115.86, "AU915", "Perth"),
    (-27.47, 153.03, "AU915", "Brisbane"),    (-41.29, 174.78, "AU915", "Wellington"),
    (-36.85, 174.76, "AU915", "Auckland"),    (-23.55, -46.63, "AU915", "SaoPaulo"),
    (-34.60, -58.38, "AU915", "BuenosAires"), (-33.45, -70.67, "AU915", "Santiago"),
    (4.71, -74.07, "AU915", "Bogota"),        (35.68, 139.69, "AS923", "Tokyo"),
    (34.69, 135.50, "AS923", "Osaka"),        (1.35, 103.82, "AS923", "Singapore"),
    (13.76, 100.50, "AS923", "Bangkok"),      (21.03, 105.85, "AS923", "Hanoi"),
    (10.82, 106.63, "AS923", "HoChiMinh"),    (-6.21, 106.85, "AS923", "Jakarta"),
    (14.60, 120.98, "AS923", "Manila"),       (3.14, 101.69, "AS923", "KualaLumpur"),
    (25.03, 121.57, "AS923", "Taipei"),       (22.32, 114.17, "AS923", "HongKong"),
    (37.57, 126.98, "KR920", "Seoul"),        (39.90, 116.41, "CN470", "Beijing"),
    (31.23, 121.47, "CN470", "Shanghai"),     (28.61, 77.21, "IN865", "Delhi"),
    (19.08, 72.88, "IN865", "Mumbai"),        (13.08, 80.27, "IN865", "Chennai"),
    (55.75, 37.62, "RU864", "Moscow"),        (59.93, 30.34, "RU864", "StPetersburg"),
    (-1.29, 36.82, "EU868", "Nairobi"),       (-26.20, 28.05, "EU868", "Johannesburg"),
    (6.52, 3.38, "EU868", "Lagos"),           (30.04, 31.24, "EU868", "Cairo"),
    (33.57, -7.59, "EU868", "Casablanca"),    (25.20, 55.27, "EU868", "Dubai"),
]


def matches(expected, got_name):
    if expected == "AS923":
        return got_name.startswith("AS923")
    return got_name == expected


def main():
    rows = hc.run_points([(lat, lon) for lat, lon, _, _ in CITIES])
    names = hc.region_names()

    correct = unknown = wrong = 0
    failures = []
    for (lat, lon, expected, city), row in zip(CITIES, rows):
        got = names.get(row["region"], f"<id {row['region']}>")
        if matches(expected, got):
            correct += 1
        elif row["region"] == 0:
            unknown += 1
            failures.append((city, expected, "Unknown"))
        else:
            wrong += 1
            failures.append((city, expected, got))

    total = len(CITIES)
    print(f"T-CITY: {correct} correct / {unknown} Unknown / {wrong} wrong "
          f"({100.0 * correct / total:.1f}% correct)")
    for city, expected, got in failures:
        print(f"  FAIL {city:15s} expected={expected:8s} got={got}")

    ok = correct >= 51 and wrong == 0
    print("T-CITY:", "PASS" if ok else "FAIL (target: >=51 correct, 0 wrong)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
