#!/usr/bin/env python3
"""
Aggregate connectivity metrics across multiple images.

Combines individual image metrics into dataset-level summaries.
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List


def load_metrics_json(json_path: Path) -> dict:
    """Load metrics from JSON file."""
    with open(json_path, 'r') as f:
        return json.load(f)


def aggregate_metrics(json_files: List[Path]) -> pd.DataFrame:
    """
    Aggregate metrics from multiple JSON files into a DataFrame.
    
    Returns DataFrame with one row per image.
    """
    records = []
    
    for json_file in json_files:
        metrics = load_metrics_json(json_file)
        
        # Extract image ID from filename
        image_id = json_file.stem.replace('_metrics', '')
        
        # Flatten nested metrics into single record
        record = {
            'image_id': image_id,
            'total_nodes': metrics['graph_info']['total_nodes'],
            'total_edges': metrics['graph_info']['total_edges'],
            'num_components': metrics['connected_components']['num_components'],
            'largest_component_nodes': metrics['connected_components']['largest_component_nodes'],
            'largest_component_length': metrics['connected_components']['largest_component_length'],
            'total_length': metrics['connected_components']['total_length'],
            'optic_disc_node': metrics['optic_disc_connectivity']['optic_disc_node'],
            'od_connected_nodes': metrics['optic_disc_connectivity']['connected_nodes'],
            'od_connected_length': metrics['optic_disc_connectivity']['connected_length'],
            'proportion_nodes_connected': metrics['optic_disc_connectivity']['proportion_nodes'],
            'proportion_length_connected': metrics['optic_disc_connectivity']['proportion_length'],
            'mean_component_size_nodes': metrics['summary_stats']['mean_component_size_nodes'],
            'median_component_size_nodes': metrics['summary_stats']['median_component_size_nodes'],
            'std_component_size_nodes': metrics['summary_stats']['std_component_size_nodes'],
            'mean_component_size_length': metrics['summary_stats']['mean_component_size_length'],
            'median_component_size_length': metrics['summary_stats']['median_component_size_length'],
            'std_component_size_length': metrics['summary_stats']['std_component_size_length'],
        }
        
        # Add relative size of largest component
        if record['total_nodes'] > 0:
            record['largest_component_proportion_nodes'] = (
                record['largest_component_nodes'] / record['total_nodes']
            )
        else:
            record['largest_component_proportion_nodes'] = 0.0
        
        if record['total_length'] > 0:
            record['largest_component_proportion_length'] = (
                record['largest_component_length'] / record['total_length']
            )
        else:
            record['largest_component_proportion_length'] = 0.0
        
        records.append(record)
    
    return pd.DataFrame(records)


def aggregate_components(csv_files: List[Path]) -> pd.DataFrame:
    """
    Aggregate component-level data from multiple CSV files.
    
    Returns DataFrame with all components from all images.
    """
    dfs = []
    
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        
        # Extract image ID from filename
        image_id = csv_file.stem.replace('_components', '')
        df['image_id'] = image_id
        
        dfs.append(df)
    
    if dfs:
        return pd.concat(dfs, ignore_index=True)
    else:
        return pd.DataFrame()


def calculate_summary_statistics(df: pd.DataFrame) -> dict:
    """
    Calculate summary statistics across all images.
    """
    stats = {}
    
    # Metrics to summarize
    metrics = [
        'num_components',
        'proportion_nodes_connected',
        'proportion_length_connected',
        'largest_component_proportion_nodes',
        'largest_component_proportion_length',
    ]
    
    for metric in metrics:
        if metric in df.columns:
            stats[metric] = {
                'mean': float(df[metric].mean()),
                'median': float(df[metric].median()),
                'std': float(df[metric].std()),
                'min': float(df[metric].min()),
                'max': float(df[metric].max()),
                'q25': float(df[metric].quantile(0.25)),
                'q75': float(df[metric].quantile(0.75)),
            }
    
    # Additional statistics
    stats['dataset_summary'] = {
        'num_images': len(df),
        'total_nodes_all_images': int(df['total_nodes'].sum()),
        'total_edges_all_images': int(df['total_edges'].sum()),
        'total_components_all_images': int(df['num_components'].sum()),
    }
    
    return stats


def main(snakemake):
    """Main function for Snakemake."""
    
    # Load all JSON files
    json_files = [Path(f) for f in snakemake.input.json_files]
    csv_files = [Path(f) for f in snakemake.input.csv_files]
    
    print(f"Aggregating metrics from {len(json_files)} images...")
    
    # Aggregate image-level metrics
    summary_df = aggregate_metrics(json_files)
    summary_df.to_csv(snakemake.output.summary, index=False)
    print(f"Saved summary to {snakemake.output.summary}")
    
    # Aggregate component-level data
    components_df = aggregate_components(csv_files)
    components_df.to_csv(snakemake.output.components, index=False)
    print(f"Saved components to {snakemake.output.components}")
    
    # Calculate summary statistics
    stats = calculate_summary_statistics(summary_df)
    with open(snakemake.output.stats, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"Saved statistics to {snakemake.output.stats}")
    
    # Print summary
    print("\n" + "="*60)
    print("DATASET SUMMARY")
    print("="*60)
    print(f"Number of images: {len(summary_df)}")
    print(f"\nConnected Components (mean ± std): {stats['num_components']['mean']:.1f} ± {stats['num_components']['std']:.1f}")
    print(f"  Range: {stats['num_components']['min']:.0f} - {stats['num_components']['max']:.0f}")
    print(f"\nProportion connected to optic disc (nodes): {stats['proportion_nodes_connected']['mean']*100:.1f}% ± {stats['proportion_nodes_connected']['std']*100:.1f}%")
    print(f"Proportion connected to optic disc (length): {stats['proportion_length_connected']['mean']*100:.1f}% ± {stats['proportion_length_connected']['std']*100:.1f}%")
    print("="*60)


if __name__ == '__main__':
    main(snakemake)
