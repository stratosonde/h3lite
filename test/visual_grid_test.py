#!/usr/bin/env python3
"""
Visual Grid Test for H3Lite Region Coverage

This script creates a world map visualization showing where h3lite
successfully finds regions vs where it doesn't (using the actual C library).
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Region colors
REGION_COLORS = {
    'Unknown': '#808080',  # Gray
    'US915': '#FF5733',    # Red
    'EU868': '#33FF57',    # Green
    'AU915': '#3357FF',    # Blue
    'AS923-1': '#FF33E9',  # Pink
    'AS923-2': '#E933FF',  # Purple
    'AS923-3': '#33E9FF',  # Cyan
    'AS923-4': '#E9FF33',  # Yellow
    'KR920': '#33FFF6',    # Light cyan
    'IN865': '#F6FF33',    # Light yellow
    'RU864': '#7F33FF',    # Purple
    'CN470': '#FF8333',    # Orange
    'EU433': '#83FF33',    # Lime
}

def get_region_from_h3lite(lat, lng):
    """Call h3lite to get the region for a lat/lng point"""
    # Create a simple C test program on the fly
    code = f"""
#include "../include/h3lite.h"
#include <stdio.h>
int main() {{
    h3liteInit();
    NearestRegionsInfo info = findNearestRegions({lat}, {lng}, 3);
    if (info.numRegions > 0) {{
        printf("%s,%d\\n", info.regions[0].regionName, info.regions[0].ringDistance);
    }} else {{
        printf("Unknown,999\\n");
    }}
    return 0;
}}
"""
    
    # Write, compile, and run
    with open('/tmp/test_point.c', 'w') as f:
        f.write(code)
    
    try:
        subprocess.run([
            'gcc', '-O3', '-I./include', '/tmp/test_point.c',
            'obj/h3lite.o', 'obj/h3lite_faceijk.o', 'obj/h3lite_basecells.o',
            'obj/h3lite_neighbor.o', 'obj/h3lite_regions_table.o',
            '-o', '/tmp/test_point', '-lm'
        ], cwd='h3lite', check=True, capture_output=True)
        
        result = subprocess.run(['/tmp/test_point'], capture_output=True, text=True)
        region, ring = result.stdout.strip().split(',')
        return region, int(ring)
    except:
        return 'Unknown', 999

def create_visual_grid(step=10):
    """Create a visual grid showing region coverage"""
    print(f"Creating visual grid with {step}° step...")
    
    lats = range(-90, 91, step)
    lngs = range(-180, 181, step)
    
    # Create grids for region and ring distance
    region_grid = []
    ring_grid = []
    
    total_points = len(lats) * len(lngs)
    tested = 0
    
    for lat in lats:
        region_row = []
        ring_row = []
        for lng in lngs:
            region, ring = get_region_from_h3lite(lat, lng)
            region_row.append(region)
            ring_row.append(ring)
            
            tested += 1
            if tested % 50 == 0:
                print(f"Tested {tested}/{total_points} points ({100*tested/total_points:.1f}%)")
        
        region_grid.append(region_row)
        ring_grid.append(ring_row)
    
    # Create the visualization
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12))
    
    # Plot 1: Region map
    region_numeric = np.zeros((len(lats), len(lngs)))
    region_names = list(REGION_COLORS.keys())
    
    for i, lat_row in enumerate(region_grid):
        for j, region in enumerate(lat_row):
            if region in region_names:
                region_numeric[i, j] = region_names.index(region)
            else:
                region_numeric[i, j] = 0  # Unknown
    
    # Create custom colormap
    colors = [REGION_COLORS[name] for name in region_names]
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(colors)
    
    im1 = ax1.imshow(region_numeric, cmap=cmap, extent=[-180, 180, -90, 90], 
                     aspect='auto', interpolation='nearest')
    ax1.set_title(f'H3Lite Region Coverage Map (Step: {step}°)', fontsize=14, weight='bold')
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')
    ax1.grid(True, linestyle='--', alpha=0.3)
    
    # Add legend
    patches = [mpatches.Patch(color=REGION_COLORS[name], label=name) for name in region_names]
    ax1.legend(handles=patches, loc='lower left', fontsize=8, ncol=3)
    
    # Plot 2: Ring distance map
    ring_numeric = np.array([[ring if ring < 900 else -1 for ring in row] for row in ring_grid])
    
    im2 = ax2.imshow(ring_numeric, cmap='RdYlGn_r', extent=[-180, 180, -90, 90],
                     aspect='auto', vmin=0, vmax=3, interpolation='nearest')
    ax2.set_title('Ring Distance to Nearest Region (0=direct, 1-3=rings away)', fontsize=14, weight='bold')
    ax2.set_xlabel('Longitude')
    ax2.set_ylabel('Latitude')
    ax2.grid(True, linestyle='--', alpha=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(im2, ax=ax2)
    cbar.set_label('Ring Distance', rotation=270, labelpad=15)
    
    plt.tight_layout()
    plt.savefig('h3lite_coverage_visual.png', dpi=150, bbox_inches='tight')
    print(f"\nSaved visualization to h3lite_coverage_visual.png")
    
    # Print statistics
    print("\n=== Coverage Statistics ===")
    total = len(lats) * len(lngs)
    ring0 = sum(1 for row in ring_grid for ring in row if ring == 0)
    ring1 = sum(1 for row in ring_grid for ring in row if ring == 1)
    ring2 = sum(1 for row in ring_grid for ring in row if ring == 2)
    ring3 = sum(1 for row in ring_grid for ring in row if ring == 3)
    not_found = sum(1 for row in ring_grid for ring in row if ring > 3)
    
    print(f"Total points tested: {total}")
    print(f"Ring 0 (direct match): {ring0} ({100*ring0/total:.1f}%)")
    print(f"Ring 1 (~65km): {ring1} ({100*ring1/total:.1f}%)")
    print(f"Ring 2 (~130km): {ring2} ({100*ring2/total:.1f}%)")
    print(f"Ring 3 (~195km): {ring3} ({100*ring3/total:.1f}%)")
    print(f"Not found (>195km): {not_found} ({100*not_found/total:.1f}%)")
    
    # Region breakdown
    print("\n=== Regions Found ===")
    for region_name in sorted(set(r for row in region_grid for r in row)):
        count = sum(1 for row in region_grid for r in row if r == region_name)
        print(f"{region_name}: {count} points ({100*count/total:.1f}%)")

if __name__ == "__main__":
    print("H3Lite Visual Grid Coverage Test")
    print("=================================\n")
    
    # Use a coarser grid (10 degrees) for faster execution
    create_visual_grid(step=10)
    
    print("\nDone!")
