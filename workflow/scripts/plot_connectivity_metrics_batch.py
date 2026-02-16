#!/usr/bin/env python3
"""
Create visualizations for batch-processed connectivity metrics.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# Snakemake inputs
summary_csv = Path(snakemake.input.summary)
stats_json = Path(snakemake.input.stats)

# Outputs
cc_dist_pdf = Path(snakemake.output.cc_dist)
size_dist_pdf = Path(snakemake.output.size_dist)
od_conn_pdf = Path(snakemake.output.od_conn)

# Ensure output directory exists
cc_dist_pdf.parent.mkdir(parents=True, exist_ok=True)

# Load data
df = pd.read_csv(summary_csv)

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

print(f"Creating plots from {len(df)} images...")

# ============================================================================
# Plot 1: Connected Components Distribution
# ============================================================================

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Histogram
ax = axes[0]
ax.hist(df['num_components'], bins=30, edgecolor='black', alpha=0.7)
ax.axvline(df['num_components'].mean(), color='red', linestyle='--', 
           label=f'Mean: {df["num_components"].mean():.1f}')
ax.axvline(df['num_components'].median(), color='blue', linestyle='--',
           label=f'Median: {df["num_components"].median():.1f}')
ax.set_xlabel('Number of Connected Components')
ax.set_ylabel('Frequency')
ax.set_title('Distribution of Connected Components')
ax.legend()
ax.grid(True, alpha=0.3)

# Box plot
ax = axes[1]
ax.boxplot([df['num_components']], labels=['All Images'])
ax.set_ylabel('Number of Connected Components')
ax.set_title('Connected Components Summary')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(cc_dist_pdf, dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Saved: {cc_dist_pdf}")

# ============================================================================
# Plot 2: Component Size Distribution
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Largest component (nodes) distribution
ax = axes[0, 0]
ax.hist(df['largest_component_proportion_nodes'], bins=30, edgecolor='black', alpha=0.7)
ax.axvline(df['largest_component_proportion_nodes'].mean(), color='red', linestyle='--',
           label=f'Mean: {df["largest_component_proportion_nodes"].mean():.2%}')
ax.set_xlabel('Largest Component (Proportion of Nodes)')
ax.set_ylabel('Frequency')
ax.set_title('Largest Component Size (Nodes)')
ax.legend()
ax.grid(True, alpha=0.3)

# Largest component (length) distribution
ax = axes[0, 1]
ax.hist(df['largest_component_proportion_length'], bins=30, edgecolor='black', alpha=0.7)
ax.axvline(df['largest_component_proportion_length'].mean(), color='red', linestyle='--',
           label=f'Mean: {df["largest_component_proportion_length"].mean():.2%}')
ax.set_xlabel('Largest Component (Proportion of Length)')
ax.set_ylabel('Frequency')
ax.set_title('Largest Component Size (Length)')
ax.legend()
ax.grid(True, alpha=0.3)

# Total nodes vs num components
ax = axes[1, 0]
ax.scatter(df['total_nodes'], df['num_components'], alpha=0.5, s=20)
ax.set_xlabel('Total Nodes')
ax.set_ylabel('Number of Connected Components')
ax.set_title('Graph Size vs Fragmentation')
ax.grid(True, alpha=0.3)

# Total length vs num components  
ax = axes[1, 1]
ax.scatter(df['total_length'], df['num_components'], alpha=0.5, s=20)
ax.set_xlabel('Total Length (pixels)')
ax.set_ylabel('Number of Connected Components')
ax.set_title('Total Length vs Fragmentation')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(size_dist_pdf, dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Saved: {size_dist_pdf}")

# ============================================================================
# Plot 3: Optic Disc Connectivity
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# OD connectivity (nodes) distribution
ax = axes[0, 0]
ax.hist(df['proportion_nodes_connected'], bins=30, edgecolor='black', alpha=0.7)
ax.axvline(df['proportion_nodes_connected'].mean(), color='red', linestyle='--',
           label=f'Mean: {df["proportion_nodes_connected"].mean():.2%}')
ax.set_xlabel('Proportion of Nodes Connected to OD')
ax.set_ylabel('Frequency')
ax.set_title('OD Connectivity (Nodes)')
ax.legend()
ax.grid(True, alpha=0.3)

# OD connectivity (length) distribution
ax = axes[0, 1]
ax.hist(df['proportion_length_connected'], bins=30, edgecolor='black', alpha=0.7)
ax.axvline(df['proportion_length_connected'].mean(), color='red', linestyle='--',
           label=f'Mean: {df["proportion_length_connected"].mean():.2%}')
ax.set_xlabel('Proportion of Length Connected to OD')
ax.set_ylabel('Frequency')
ax.set_title('OD Connectivity (Length)')
ax.legend()
ax.grid(True, alpha=0.3)

# Number of CCs touching OD
ax = axes[1, 0]
ax.hist(df['num_components_touching_od'], bins=20, edgecolor='black', alpha=0.7)
ax.axvline(df['num_components_touching_od'].mean(), color='red', linestyle='--',
           label=f'Mean: {df["num_components_touching_od"].mean():.1f}')
ax.set_xlabel('Number of Components Touching OD')
ax.set_ylabel('Frequency')
ax.set_title('Components Connected to OD')
ax.legend()
ax.grid(True, alpha=0.3)

# Nodes vs Length connectivity comparison
ax = axes[1, 1]
ax.scatter(df['proportion_nodes_connected'], df['proportion_length_connected'], 
           alpha=0.5, s=20)
ax.plot([0, 1], [0, 1], 'r--', alpha=0.5, label='1:1 line')
ax.set_xlabel('Proportion Nodes Connected')
ax.set_ylabel('Proportion Length Connected')
ax.set_title('Node vs Length Connectivity')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

plt.tight_layout()
plt.savefig(od_conn_pdf, dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Saved: {od_conn_pdf}")

print(f"\n✓ All plots created successfully!")
