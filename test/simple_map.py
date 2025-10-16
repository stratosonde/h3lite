#!/usr/bin/env python3
"""
Simple map visualization of h3lite region coverage
Reads from grid test output and creates a visual map
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import subprocess

# Run grid test and capture output
result = subprocess.run(
    ['./bin/h3lite_grid_test'],
    cwd='h3lite',
    capture_output=True,
    text=True
)

# Parse output
grid_data = {}
for line in result.stdout.split('\n'):
    if line.strip() and not line.startswith('=') and not line.startswith('Testing'):
        parts = line.split()
        if len(parts) >= 6 and parts[0].lstrip('-').isdigit():
            try:
                lat = int(parts[0])
                lng = int(parts[1])
                region = parts[5]
                grid_data[(lat, lng)] = region
            except:
                pass

# Create grid arrays
lats = sorted(set(lat for lat, lng in grid_data.keys()))
lngs = sorted(set(lng for lat, lng in grid_data.keys()))

print(f"Grid size: {len(lats)} x {len(lngs)} = {len(grid_data)} points")

# Create numeric grid
region_map = {
    'Unknown': 0, 'US915': 1, 'EU868': 2, 'AU915': 3, 'AS923-1': 4,
    'AS923-2': 5, 'AS923-3': 6, 'AS923-4': 7, 'KR920': 8, 'IN865': 9,
    'RU864': 10, 'CN470': 11, 'EU433': 12
}

colors = {
    'Unknown': '#CCCCCC', 'US915': '#FF5733', 'EU868': '#33FF57',
    'AU915': '#3357FF', 'AS923-1': '#FF33E9', 'AS923-2': '#E933FF',
    'AS923-3': '#33E9FF', 'AS923-4': '#E9FF33', 'KR920': '#33FFF6',
    'IN865': '#F6FF33', 'RU864': '#7F33FF', 'CN470': '#FF8333', 'EU433': '#83FF33'
}

grid = np.zeros((len(lats), len(lngs)))
for i, lat in enumerate(lats):
    for j, lng in enumerate(lngs):
        region = grid_data.get((lat, lng), 'Unknown')
        grid[i, j] = region_map.get(region, 0)

# Create plot
fig, ax = plt.subplots(figsize=(16, 10))

from matplotlib.colors import ListedColormap
cmap = ListedColormap([colors[r] for r in region_map.keys()])

im = ax.imshow(grid, cmap=cmap, extent=[-180, 180, -90, 90],
               aspect='auto', interpolation='nearest', origin='lower')

ax.set_title('H3Lite Region Coverage (10° Grid)', fontsize=16, weight='bold')
ax.set_xlabel('Longitude', fontsize=12)
ax.set_ylabel('Latitude', fontsize=12)
ax.grid(True, linestyle='--', alpha=0.3, color='white', linewidth=0.5)

# Add legend
patches = [mpatches.Patch(color=colors[name], label=name) 
           for name in region_map.keys() if name in colors]
ax.legend(handles=patches, loc='lower left', fontsize=9, ncol=2, 
          framealpha=0.9, edgecolor='black')

# Add some reference lines
ax.axhline(y=0, color='white', linestyle='-', linewidth=1, alpha=0.5)
ax.axvline(x=0, color='white', linestyle='-', linewidth=1, alpha=0.5)

plt.tight_layout()
plt.savefig('h3lite_region_map.png', dpi=200, bbox_inches='tight')
print("Saved map to h3lite_region_map.png")

# Statistics
total = len(grid_data)
by_region = {}
for region in grid_data.values():
    by_region[region] = by_region.get(region, 0) + 1

print("\n=== Coverage Statistics ===")
for region in sorted(by_region.keys()):
    count = by_region[region]
    pct = 100.0 * count / total
    print(f"{region:12s}: {count:3d} points ({pct:5.1f}%)")

known_regions = sum(c for r, c in by_region.items() if r != 'Unknown')
print(f"\nTotal coverage: {known_regions}/{total} ({100.0*known_regions/total:.1f}%)")
