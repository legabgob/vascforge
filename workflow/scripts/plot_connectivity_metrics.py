#!/usr/bin/env python3
"""
Create visualization plots for connectivity metrics.

Generates distribution plots for:
- Connected component counts
- Component sizes (nodes and length)  
- Optic disc connectivity proportions
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


def plot_cc_distribution(summary_df: pd.DataFrame, output_path: Path):
    """
    Plot distribution of connected component counts across images.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Histogram
    ax1 = axes[0]
    summary_df['num_components'].hist(bins=30, ax=ax1, edgecolor='black', alpha=0.7)
    ax1.axvline(summary_df['num_components'].mean(), color='r', linestyle='--', 
               linewidth=2, label=f'Mean: {summary_df["num_components"].mean():.1f}')
    ax1.axvline(summary_df['num_components'].median(), color='g', linestyle='--',
               linewidth=2, label=f'Median: {summary_df["num_components"].median():.1f}')
    ax1.set_xlabel('Number of Connected Components', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('Distribution of Connected Component Counts', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Box plot
    ax2 = axes[1]
    box_data = summary_df['num_components']
    bp = ax2.boxplot([box_data], labels=['All Images'], patch_artist=True)
    bp['boxes'][0].set_facecolor('lightblue')
    bp['boxes'][0].set_alpha(0.7)
    
    # Add individual points
    y = np.random.normal(1, 0.04, size=len(box_data))
    ax2.scatter(y, box_data, alpha=0.3, s=20)
    
    ax2.set_ylabel('Number of Connected Components', fontsize=12)
    ax2.set_title('Connected Component Count Distribution', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Add statistics text
    stats_text = f'n = {len(summary_df)}\n'
    stats_text += f'Mean ± SD: {summary_df["num_components"].mean():.1f} ± {summary_df["num_components"].std():.1f}\n'
    stats_text += f'Median [IQR]: {summary_df["num_components"].median():.1f} '
    stats_text += f'[{summary_df["num_components"].quantile(0.25):.1f}-{summary_df["num_components"].quantile(0.75):.1f}]'
    ax2.text(0.5, 0.98, stats_text, transform=ax2.transAxes,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
            fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved CC distribution plot to {output_path}")


def plot_component_sizes(components_df: pd.DataFrame, output_path: Path):
    """
    Plot distribution of component sizes (both nodes and length).
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Rank-size plot for nodes
    ax1 = axes[0, 0]
    for image_id in components_df['image_id'].unique()[:50]:  # Plot first 50 images
        img_data = components_df[components_df['image_id'] == image_id].sort_values(
            'num_nodes', ascending=False
        )
        ax1.plot(range(1, len(img_data)+1), img_data['num_nodes'], 
                alpha=0.3, linewidth=1)
    
    ax1.set_xlabel('Component Rank', fontsize=12)
    ax1.set_ylabel('Number of Nodes', fontsize=12)
    ax1.set_title('Component Size Distribution (Nodes)', fontsize=13, fontweight='bold')
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3)
    
    # Rank-size plot for length
    ax2 = axes[0, 1]
    for image_id in components_df['image_id'].unique()[:50]:
        img_data = components_df[components_df['image_id'] == image_id].sort_values(
            'total_length', ascending=False
        )
        ax2.plot(range(1, len(img_data)+1), img_data['total_length'],
                alpha=0.3, linewidth=1)
    
    ax2.set_xlabel('Component Rank', fontsize=12)
    ax2.set_ylabel('Total Segment Length', fontsize=12)
    ax2.set_title('Component Size Distribution (Length)', fontsize=13, fontweight='bold')
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3)
    
    # Relative size distribution (nodes)
    ax3 = axes[1, 0]
    if 'relative_nodes' in components_df.columns:
        # Only plot non-largest components
        non_largest = components_df[components_df['rank_by_nodes'] > 1]
        non_largest['relative_nodes'].hist(bins=50, ax=ax3, edgecolor='black', alpha=0.7)
        ax3.set_xlabel('Relative Size (Nodes)', fontsize=12)
        ax3.set_ylabel('Frequency', fontsize=12)
        ax3.set_title('Distribution of Relative Component Sizes (Nodes)\nExcluding Largest Component',
                     fontsize=13, fontweight='bold')
        ax3.set_xlim(0, 0.2)  # Focus on smaller components
        ax3.grid(True, alpha=0.3, axis='y')
    
    # Relative size distribution (length)
    ax4 = axes[1, 1]
    if 'relative_length' in components_df.columns:
        non_largest = components_df[components_df['rank_by_length'] > 1]
        non_largest['relative_length'].hist(bins=50, ax=ax4, edgecolor='black', alpha=0.7)
        ax4.set_xlabel('Relative Size (Length)', fontsize=12)
        ax4.set_ylabel('Frequency', fontsize=12)
        ax4.set_title('Distribution of Relative Component Sizes (Length)\nExcluding Largest Component',
                     fontsize=13, fontweight='bold')
        ax4.set_xlim(0, 0.2)
        ax4.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved component size distribution plot to {output_path}")


def plot_optic_disc_connectivity(summary_df: pd.DataFrame, output_path: Path):
    """
    Plot distribution of optic disc connectivity proportions.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Proportion of nodes connected - histogram
    ax1 = axes[0, 0]
    summary_df['proportion_nodes_connected'].hist(bins=30, ax=ax1, 
                                                   edgecolor='black', alpha=0.7, color='steelblue')
    ax1.axvline(summary_df['proportion_nodes_connected'].mean(), color='r', 
               linestyle='--', linewidth=2, 
               label=f'Mean: {summary_df["proportion_nodes_connected"].mean():.3f}')
    ax1.set_xlabel('Proportion of Nodes Connected to Optic Disc', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('Distribution of Optic Disc Connectivity (Nodes)', 
                 fontsize=13, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Proportion of length connected - histogram
    ax2 = axes[0, 1]
    summary_df['proportion_length_connected'].hist(bins=30, ax=ax2,
                                                    edgecolor='black', alpha=0.7, color='coral')
    ax2.axvline(summary_df['proportion_length_connected'].mean(), color='r',
               linestyle='--', linewidth=2,
               label=f'Mean: {summary_df["proportion_length_connected"].mean():.3f}')
    ax2.set_xlabel('Proportion of Length Connected to Optic Disc', fontsize=12)
    ax2.set_ylabel('Frequency', fontsize=12)
    ax2.set_title('Distribution of Optic Disc Connectivity (Length)',
                 fontsize=13, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Scatter plot: nodes vs length connectivity
    ax3 = axes[1, 0]
    ax3.scatter(summary_df['proportion_nodes_connected'],
               summary_df['proportion_length_connected'],
               alpha=0.5, s=50, color='purple')
    
    # Add diagonal line
    ax3.plot([0, 1], [0, 1], 'r--', alpha=0.5, linewidth=2, label='y=x')
    
    ax3.set_xlabel('Proportion of Nodes Connected', fontsize=12)
    ax3.set_ylabel('Proportion of Length Connected', fontsize=12)
    ax3.set_title('Nodes vs Length Connectivity Comparison',
                 fontsize=13, fontweight='bold')
    ax3.legend(fontsize=11)
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, 1)
    
    # Box plots comparing both metrics
    ax4 = axes[1, 1]
    data_to_plot = [
        summary_df['proportion_nodes_connected'],
        summary_df['proportion_length_connected']
    ]
    bp = ax4.boxplot(data_to_plot, labels=['Nodes', 'Length'],
                    patch_artist=True)
    
    colors = ['steelblue', 'coral']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Add individual points
    for i, data in enumerate(data_to_plot):
        y = np.random.normal(i+1, 0.04, size=len(data))
        ax4.scatter(y, data, alpha=0.3, s=20)
    
    ax4.set_ylabel('Proportion Connected to Optic Disc', fontsize=12)
    ax4.set_title('Comparison of Connectivity Metrics',
                 fontsize=13, fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='y')
    ax4.set_ylim(0, 1)
    
    # Add statistics
    stats_text = f'Nodes: {summary_df["proportion_nodes_connected"].mean():.3f} ± {summary_df["proportion_nodes_connected"].std():.3f}\n'
    stats_text += f'Length: {summary_df["proportion_length_connected"].mean():.3f} ± {summary_df["proportion_length_connected"].std():.3f}'
    ax4.text(0.5, 0.02, stats_text, transform=ax4.transAxes,
            verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
            fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved optic disc connectivity plot to {output_path}")


def main(snakemake):
    """Main function for Snakemake."""
    
    # Load data
    print("Loading data...")
    summary_df = pd.read_csv(snakemake.input.summary)
    components_df = pd.read_csv(snakemake.input.components)
    
    print(f"Loaded data: {len(summary_df)} images, {len(components_df)} components")
    
    # Create output directory
    Path(snakemake.output.cc_dist).parent.mkdir(parents=True, exist_ok=True)
    
    # Generate plots
    print("Creating connected component distribution plot...")
    plot_cc_distribution(summary_df, Path(snakemake.output.cc_dist))
    
    print("Creating component size distribution plot...")
    plot_component_sizes(components_df, Path(snakemake.output.size_dist))
    
    print("Creating optic disc connectivity plot...")
    plot_optic_disc_connectivity(summary_df, Path(snakemake.output.od_conn))
    
    print("All plots created successfully!")


if __name__ == '__main__':
    main(snakemake)
