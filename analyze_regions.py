#!/usr/bin/env python3
"""
H3Lite Region Analysis Tool

This script helps analyze the H3Lite region lookup table and visualize
the boundaries of each region. It can help identify issues with the
region detection, such as missing Sydney/AU915 mapping.
"""

import sys
import json
import os
import h3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import pandas as pd

# Define regions with colors for visualization
REGIONS = {
    'US915': '#FF5733',
    'EU868': '#33FF57',
    'AU915': '#3357FF',
    'AS923-1': '#FF33E9',
    'KR920': '#33FFF6',
    'IN865': '#F6FF33',
    'RU864': '#7F33FF',
    'CN470': '#FF8333',
    'Unknown': '#CCCCCC'
}

# Test points from h3lite_test.c
TEST_POINTS = [
    {'name': 'San Francisco', 'lat': 37.7749, 'lng': -122.4194, 'expected_region': 'US915'},
    {'name': 'Paris', 'lat': 48.8566, 'lng': 2.3522, 'expected_region': 'EU868'},
    {'name': 'Sydney', 'lat': -33.8688, 'lng': 151.2093, 'expected_region': 'AU915'},
    {'name': 'Tokyo', 'lat': 35.6762, 'lng': 139.6503, 'expected_region': 'AS923-1'},
    {'name': 'Seoul', 'lat': 37.5665, 'lng': 126.9780, 'expected_region': 'KR920'},
    {'name': 'New Delhi', 'lat': 28.6139, 'lng': 77.2090, 'expected_region': 'IN865'},
    {'name': 'Moscow', 'lat': 55.7558, 'lng': 37.6173, 'expected_region': 'RU864'},
    {'name': 'Beijing', 'lat': 39.9042, 'lng': 116.4074, 'expected_region': 'CN470'}
]

def lat_lng_to_h3(lat, lng, resolution=3):
    """Convert latitude/longitude to H3 index"""
    # The h3 Python library seems to use latlng_to_cell instead of geo_to_h3
    try:
        # Try newer version API
        h3_index = h3.latlng_to_cell(lat, lng, resolution)
    except AttributeError:
        # Try older version API
        try:
            h3_index = h3.geo_to_h3(lat, lng, resolution)
        except AttributeError:
            # Older version of h3py used lat_lng_to_h3_address
            h3_index = h3.lat_lng_to_h3_address(lat, lng, resolution)
    return h3_index

def h3_to_components(h3_index):
    """Extract base cell and partial index from an H3 index"""
    h3_int = int(h3_index, 16)
    
    # Extract resolution
    res = (h3_int >> 52) & 0xF
    
    # Extract base cell (bits 45-52)
    base_cell = (h3_int >> 45) & 0x7F
    
    # Extract first 3 resolution digits
    partial_index = 0
    for r in range(1, min(4, res+1)):
        shift = 45 - (3 * r)
        digit = (h3_int >> shift) & 0x7
        partial_index = (partial_index * 8) + digit
    
    return {
        'h3_index': h3_index,
        'resolution': res,
        'base_cell': base_cell,
        'partial_index': partial_index
    }

def analyze_test_points():
    """Analyze test points and their H3 indices"""
    print("\nTest Point Analysis (Resolution 3):")
    print("-" * 80)
    print(f"{'Location':<15} {'H3 Index':<20} {'Base Cell':<10} {'Partial Index':<15} {'Expected Region':<15}")
    print("-" * 80)
    
    for point in TEST_POINTS:
        h3_index = lat_lng_to_h3(point['lat'], point['lng'], 3)
        components = h3_to_components(h3_index)
        
        print(f"{point['name']:<15} {h3_index:<20} {components['base_cell']:<10} {components['partial_index']:<15} {point['expected_region']:<15}")
        
        # Special focus on Sydney/AU915
        if point['name'] == 'Sydney':
            print("\nDetailed Analysis for Sydney:")
            print(f"Lat/Lng: {point['lat']}, {point['lng']}")
            print(f"H3 Index: {h3_index}")
            print(f"Resolution: {components['resolution']}")
            print(f"Base Cell: {components['base_cell']}")
            print(f"Partial Index: {components['partial_index']}")
            
            # Generate nearby indices
            print("\nNearby H3 cells around Sydney:")
            try:
                # Try newer version API
                neighbors = list(h3.grid_ring(h3_index, 1))
            except AttributeError:
                # Try older version API
                try:
                    neighbors = list(h3.k_ring(h3_index, 1))
                except AttributeError:
                    # Another possible API name
                    neighbors = list(h3.hex_ring(h3_index, 1))
                
            for i, neighbor in enumerate(neighbors):
                n_comp = h3_to_components(neighbor)
                print(f"  Neighbor {i}: {neighbor} (BaseCell: {n_comp['base_cell']}, PartialIndex: {n_comp['partial_index']})")

def visualize_regions(resolution=3):
    """Create a map visualization of the regions"""
    print("\nGenerating region visualization map...")
    
    # Create a world grid at the specified resolution
    lats = np.linspace(-90, 90, 180)
    lngs = np.linspace(-180, 180, 360)
    
    grid = np.zeros((len(lats), len(lngs)), dtype=int)
    
    # Sample points from the grid and assign regions
    for i, lat in enumerate(lats):
        for j, lng in enumerate(lngs):
            # Skip some points to speed up visualization
            if i % 5 != 0 or j % 5 != 0:
                continue
                
            # Find the nearest test point as a simple region assignment
            min_dist = float('inf')
            region_idx = 8  # Default to "Unknown" (index 8 in our color list)
            
            for k, point in enumerate(TEST_POINTS):
                dist = (lat - point['lat'])**2 + (lng - point['lng'])**2
                if dist < min_dist:
                    min_dist = dist
                    region_idx = k
            
            grid[i, j] = region_idx
    
    # Create the plot
    plt.figure(figsize=(15, 10))
    
    # Create color map from region colors
    colors = list(REGIONS.values())
    cmap = ListedColormap(colors)
    
    # Plot the grid
    plt.imshow(grid, cmap=cmap, extent=[-180, 180, -90, 90], aspect='auto')
    
    # Add test points
    for point in TEST_POINTS:
        plt.plot(point['lng'], point['lat'], 'o', markersize=10, 
                color='black', markerfacecolor='white')
        plt.text(point['lng'], point['lat'] + 3, point['name'], 
                 ha='center', va='bottom', fontsize=10, weight='bold')
    
    # Add legend
    patches = []
    for i, (region, color) in enumerate(REGIONS.items()):
        patches.append(mpatches.Patch(color=color, label=region))
    
    plt.legend(handles=patches, loc='lower left', framealpha=0.7)
    
    # Set plot attributes
    plt.title(f'LoRaWAN Regions (H3 Resolution {resolution})')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Save the plot
    plt.savefig('region_visualization.png', dpi=300, bbox_inches='tight')
    print("Visualization saved as 'region_visualization.png'")

def analyze_lookup_table():
    """
    Generate information about what would be in the lookup table
    for the AU915 region to help diagnose why Sydney isn't found
    """
    print("\nLookup Table Analysis for AU915 Region:")
    print("-" * 80)
    
    # First, let's confirm the information for Sydney
    sydney = next(p for p in TEST_POINTS if p['name'] == 'Sydney')
    sydney_h3 = lat_lng_to_h3(sydney['lat'], sydney['lng'], 3)
    sydney_comp = h3_to_components(sydney_h3)
    
    print(f"Sydney H3 Index: {sydney_h3}")
    print(f"Base Cell: {sydney_comp['base_cell']}")
    print(f"Partial Index: {sydney_comp['partial_index']}")
    
    # Create a recommendation for the lookup table
    print("\nRecommended entry for the region lookup table:")
    print("```c")
    print(f"{{ {sydney_comp['base_cell']}, {sydney_comp['partial_index']}, REGION_AU915 }}, // Sydney")
    print("```")
    
    # Generate a sample patch for the h3lite_regions_table.c file
    print("\nSample patch for h3lite_regions_table.c:")
    print("```c")
    # Find an appropriate place to insert (somewhere in the sorted table)
    print(f"    // Insert near other entries with baseCell={sydney_comp['base_cell']}")
    print(f"    {{ {sydney_comp['base_cell']}, {sydney_comp['partial_index']}, 3 }}, // AU915 (Sydney)")
    print("```")

if __name__ == "__main__":
    print("H3Lite Region Analysis Tool")
    print("=========================")
    
    analyze_test_points()
    analyze_lookup_table()
    visualize_regions()
    
    print("\nAnalysis complete!")
