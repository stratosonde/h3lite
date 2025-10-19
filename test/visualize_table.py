#!/usr/bin/env python3
"""
Visualize what's actually in the h3lite regions lookup table
"""

import re
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

# Parse the regions table
def parse_regions_table():
    with open('src/h3lite_regions_table.c', 'r') as f:
        content = f.read()
    
    # Extract entries
    entries = []
    for line in content.split('\n'):
        if '{ ' in line and ', ' in line and '}, //' in line:
            # Parse: { baseCell, partialIndex, regionId }, // RegionName
            match = re.search(r'\{\s*(\d+),\s*(\d+),\s*(\d+)\s*\},\s*//\s*(\w+[-\w]*)', line)
            if match:
                entries.append({
                    'baseCell': int(match.group(1)),
                    'partialIndex': int(match.group(2)),
                    'regionId': int(match.group(3)),
                    'regionName': match.group(4)
                })
    
    return entries

# Analyze the table
entries = parse_regions_table()
print(f"Total entries in lookup table: {len(entries)}")

# Group by region
by_region = defaultdict(list)
for entry in entries:
    by_region[entry['regionName']].append(entry)

print("\nEntries per region:")
for region, items in sorted(by_region.items(), key=lambda x: len(x[1]), reverse=True):
    print(f"  {region}: {len(items)} entries")

# Group by base cell
by_basecell = defaultdict(list)
for entry in entries:
    by_basecell[entry['baseCell']].append(entry)

print(f"\nNumber of base cells with entries: {len(by_basecell)}")
print(f"Base cells used: {sorted(by_basecell.keys())}")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Plot 1: Entries per region
ax1 = axes[0, 0]
regions = list(by_region.keys())
counts = [len(by_region[r]) for r in regions]
ax1.barh(regions, counts)
ax1.set_xlabel('Number of Entries')
ax1.set_title('Lookup Table Entries per Region')
ax1.grid(axis='x', alpha=0.3)

# Plot 2: Entries per base cell
ax2 = axes[0, 1]
basecells = sorted(by_basecell.keys())
bc_counts = [len(by_basecell[bc]) for bc in basecells]
ax2.bar(range(len(basecells)), bc_counts)
ax2.set_xlabel('Base Cell')
ax2.set_ylabel('Number of Entries')
ax2.set_title(f'Entries per Base Cell ({len(basecells)} base cells used)')
ax2.grid(axis='y', alpha=0.3)

# Plot 3: Partial index distribution
ax3 = axes[1, 0]
partial_indices = [e['partialIndex'] for e in entries]
ax3.hist(partial_indices, bins=50, edgecolor='black')
ax3.set_xlabel('Partial Index Value')
ax3.set_ylabel('Count')
ax3.set_title('Distribution of Partial Index Values')
ax3.grid(axis='y', alpha=0.3)

# Plot 4: Sample of actual entries
ax4 = axes[1, 1]
ax4.axis('off')
sample_text = "Sample Lookup Table Entries:\n\n"
for i, entry in enumerate(entries[:20]):
    sample_text += f"BC:{entry['baseCell']:3d} PI:{entry['partialIndex']:4d} -> {entry['regionName']}\n"
sample_text += f"\n... and {len(entries)-20} more entries"
ax4.text(0.1, 0.9, sample_text, fontsize=9, family='monospace',
         verticalalignment='top', transform=ax4.transAxes)
ax4.set_title('Sample Entries')

plt.tight_layout()
plt.savefig('h3lite_table_analysis.png', dpi=150, bbox_inches='tight')
print("\nSaved analysis to h3lite_table_analysis.png")

# Now let's check a specific point
print("\n=== Debug Specific Point ===")
print("San Francisco: 37.7749, -122.4194")
print("From DEBUG output: baseCell=20, partialIndex=42")
print("\nLooking for entries with baseCell=20:")
bc20_entries = [e for e in entries if e['baseCell'] == 20]
print(f"Found {len(bc20_entries)} entries with baseCell=20")
for e in sorted(bc20_entries, key=lambda x: x['partialIndex'])[:10]:
    print(f"  PI:{e['partialIndex']:4d} -> {e['regionName']}")

print("\nLooking for entries near partialIndex=42:")
near_entries = [e for e in entries if 30 <= e['partialIndex'] <= 60]
print(f"Found {len(near_entries)} entries with PI between 30-60")
for e in sorted(near_entries, key=lambda x: (x['baseCell'], x['partialIndex']))[:20]:
    print(f"  BC:{e['baseCell']:3d} PI:{e['partialIndex']:4d} -> {e['regionName']}")
