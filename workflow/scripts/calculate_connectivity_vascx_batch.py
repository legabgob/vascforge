#!/usr/bin/env python3
"""
Calculate connectivity metrics for ENTIRE VascX dataset (batch processing).
Avoids per-image job explosion by processing all images at once.
"""

import argparse
import errno
import gc
import json
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
from typing import Dict, List
import warnings
warnings.filterwarnings('ignore')

# VascX imports
from vascx.fundus.loader import RetinaLoader
from vascx.fundus.layer import VesselTreeLayer


# ============================================================================
# Graph utilities (copied from before)
# ============================================================================

def get_seg_ids_from_cc(graph):
    """Get list of segment IDs from a connected component."""
    list_seg_ids = []
    for s, e in graph.edges():
        list_seg_ids.append(frozenset([s, e]))
    return list_seg_ids


def get_segments_from_cc(graph, dict_seg: dict):
    """Get list of segments from a connected component."""
    seg_list = []
    list_seg_ids = get_seg_ids_from_cc(graph)
    
    for seg_id in list_seg_ids:
        if seg_id in dict_seg:
            seg_list.append(dict_seg[seg_id])
    
    return seg_list


def get_graph_total_length(graph, dict_seg: dict):
    """Calculate total length of segments in a connected component."""
    seg_list = get_segments_from_cc(graph, dict_seg)
    return sum([seg.length for seg in seg_list])


def get_connected_components(layer: VesselTreeLayer, sort: bool = True, dict_seg: dict = None):
    """Get connected components from a vessel layer.

    Uses subgraph *views* (no copying) to avoid duplicating graph data in
    memory.  Pass a pre-built ``dict_seg`` mapping to avoid rebuilding it
    when sorting is required.
    """
    graph = layer.digraph
    
    if nx.is_directed(graph):
        components = nx.weakly_connected_components(graph)
    else:
        components = nx.connected_components(graph)
    
    # Use views instead of copies to avoid duplicating graph data in memory
    cc_list = [graph.subgraph(c) for c in components]
    
    if sort:
        if dict_seg is None:
            dict_seg = {seg.id: seg for seg in layer.segments}
        return sorted(cc_list, key=lambda x: get_graph_total_length(x, dict_seg), reverse=True)
    else:
        return cc_list


def cc_out_of_optic_disc(layer: VesselTreeLayer, cc):
    """Check if connected component touches the optic disc."""
    pts = []
    for node_idx in cc.nodes:
        point = (cc.nodes[node_idx]['o'][0], cc.nodes[node_idx]['o'][1])
        pts.append(point)
    
    return layer.disc.intersects(pts, with_boundary=True)


def get_cc_subset_touching_od(layer: VesselTreeLayer, cc_list):
    """Get connected components that touch the optic disc."""
    return [cc for cc in cc_list if cc_out_of_optic_disc(layer, cc)]


# ============================================================================
# Per-image metrics calculation
# ============================================================================

def calculate_image_metrics(layer: VesselTreeLayer):
    """
    Calculate connectivity metrics for a single image.
    Returns dict with key metrics.
    """
    # Build the segment lookup dict once per image to avoid rebuilding it for
    # every connected component during sorting and length calculations.
    dict_seg = {seg.id: seg for seg in layer.segments}

    # Get connected components (subgraph views, no copies)
    cc_list = get_connected_components(layer, sort=True, dict_seg=dict_seg)
    num_components = len(cc_list)
    
    if num_components == 0:
        return {
            'num_components': 0,
            'total_nodes': 0,
            'total_edges': 0,
            'total_length': 0,
            'largest_component_nodes': 0,
            'largest_component_length': 0,
            'num_components_touching_od': 0,
            'proportion_nodes_connected': 0.0,
            'proportion_length_connected': 0.0,
        }
    
    # Calculate component sizes
    component_sizes_nodes = [cc.number_of_nodes() for cc in cc_list]
    component_sizes_length = [get_graph_total_length(cc, dict_seg) for cc in cc_list]
    
    total_nodes = layer.graph.number_of_nodes()
    total_edges = layer.graph.number_of_edges()
    total_length = sum(component_sizes_length)
    
    # Optic disc connectivity
    cc_touching_od = get_cc_subset_touching_od(layer, cc_list)
    connected_nodes = sum(cc.number_of_nodes() for cc in cc_touching_od)
    connected_length = sum(get_graph_total_length(cc, dict_seg) for cc in cc_touching_od)
    
    return {
        'num_components': num_components,
        'total_nodes': total_nodes,
        'total_edges': total_edges,
        'total_length': total_length,
        'largest_component_nodes': component_sizes_nodes[0] if component_sizes_nodes else 0,
        'largest_component_length': component_sizes_length[0] if component_sizes_length else 0,
        'largest_component_proportion_nodes': component_sizes_nodes[0] / total_nodes if total_nodes > 0 else 0,
        'largest_component_proportion_length': component_sizes_length[0] / total_length if total_length > 0 else 0,
        'num_components_touching_od': len(cc_touching_od),
        'connected_nodes': connected_nodes,
        'connected_length': connected_length,
        'proportion_nodes_connected': connected_nodes / total_nodes if total_nodes > 0 else 0,
        'proportion_length_connected': connected_length / total_length if total_length > 0 else 0,
    }


# ============================================================================
# Batch processing for entire dataset
# ============================================================================

def process_dataset_batch(
    dataset_path: Path,
    output_dir: Path,
    av_subfolder: str = 'av',
    fundus_subfolder: str = 'rgb',
    batch_size: int = 250
):
    """
    Process all images in a VascX dataset and output summary statistics.

    Images are processed in fixed-size batches (controlled by ``batch_size``)
    to bound peak memory usage.  Partial results for each batch are written to
    a temporary CSV inside ``output_dir`` and then concatenated at the end, so
    that the in-memory ``results`` list never grows beyond ``batch_size`` rows.
    """
    print(f"="*70)
    print(f"BATCH CONNECTIVITY ANALYSIS")
    print(f"="*70)
    print(f"Dataset: {dataset_path}")
    print(f"Batch size: {batch_size}")

    # Load dataset
    loader = RetinaLoader.from_folder(
        dataset_path,
        av_subfolder=av_subfolder,
        fundus_subfolder=fundus_subfolder
    )

    n_total = len(loader)
    print(f"Found {n_total} images to process")

    output_dir.mkdir(parents=True, exist_ok=True)
    partial_dir = output_dir / "_partial_batches"
    partial_dir.mkdir(parents=True, exist_ok=True)

    failed = []
    n_processed = 0

    # Split indices into fixed-size batches
    batch_indices = [
        range(start, min(start + batch_size, n_total))
        for start in range(0, n_total, batch_size)
    ]
    n_batches = len(batch_indices)

    for batch_num, indices in enumerate(batch_indices, start=1):
        batch_results = []
        print(f"\n--- Batch {batch_num}/{n_batches} (images {indices.start + 1}–{indices.stop}) ---")

        # Use index-based iteration so loader crashes (e.g. NaN in disc/fovea)
        # are caught per-image and don't kill the entire batch job.
        for i in indices:
            image_id = f"index_{i}"  # fallback ID if loading fails before we get retina.id
            retina = None
            try:
                retina = loader[i]  # ← inside try/except so VascX init errors are caught
                image_id = retina.id
                print(f"[{i+1}/{n_total}] Processing {image_id}...", end=' ')

                # Use arteries layer (where combined vessel masks are loaded by VascX)
                layer = retina.arteries
                if layer is None or layer.graph is None or layer.graph.number_of_nodes() == 0:
                    layer = retina.veins  # fallback
                    if layer is None or layer.graph is None or layer.graph.number_of_nodes() == 0:
                        raise ValueError("No vessel data found in arteries or veins layer")

                # Calculate metrics
                metrics = calculate_image_metrics(layer)
                metrics['image_id'] = image_id

                batch_results.append(metrics)
                print(f"✓ ({metrics['num_components']} CCs)")

            except Exception as e:
                print(f"✗ FAILED: {e}")
                failed.append({
                    'image_id': image_id,
                    'error': str(e)
                })
            finally:
                # Explicitly release the retina object and its graph data so
                # memory from each image is freed before loading the next one.
                del retina
                gc.collect()

        # Write this batch's results to a partial CSV, then free memory
        if batch_results:
            n_processed += len(batch_results)
            partial_csv = partial_dir / f"batch_{batch_num:04d}.csv"
            pd.DataFrame(batch_results).to_csv(partial_csv, index=False)
            print(f"  → Saved {len(batch_results)} rows to {partial_csv}")

        del batch_results
        gc.collect()

    print(f"\nProcessed {n_processed} images successfully")
    if failed:
        print(f"Failed: {len(failed)} images")

    # Concatenate all partial CSVs into the final DataFrame
    partial_csvs = sorted(partial_dir.glob("batch_[0-9]*.csv"))
    if partial_csvs:
        df = pd.concat([pd.read_csv(p) for p in partial_csvs], ignore_index=True)
    else:
        df = pd.DataFrame()

    # Clean up partial batch files now that they have been merged
    for p in partial_csvs:
        p.unlink()
    try:
        partial_dir.rmdir()
    except OSError as exc:
        if exc.errno != errno.ENOTEMPTY:
            raise
    
    # Calculate summary statistics
    summary_stats = {
        'dataset_info': {
            'dataset_path': str(dataset_path),
            'n_images': n_processed,
            'n_failed': len(failed),
        },
        'connectivity_stats': {
            'num_components': {
                'mean': float(df['num_components'].mean()),
                'median': float(df['num_components'].median()),
                'std': float(df['num_components'].std()),
                'min': float(df['num_components'].min()),
                'max': float(df['num_components'].max()),
                'q25': float(df['num_components'].quantile(0.25)),
                'q75': float(df['num_components'].quantile(0.75)),
            },
            'proportion_nodes_connected': {
                'mean': float(df['proportion_nodes_connected'].mean()),
                'median': float(df['proportion_nodes_connected'].median()),
                'std': float(df['proportion_nodes_connected'].std()),
                'min': float(df['proportion_nodes_connected'].min()),
                'max': float(df['proportion_nodes_connected'].max()),
            },
            'proportion_length_connected': {
                'mean': float(df['proportion_length_connected'].mean()),
                'median': float(df['proportion_length_connected'].median()),
                'std': float(df['proportion_length_connected'].std()),
                'min': float(df['proportion_length_connected'].min()),
                'max': float(df['proportion_length_connected'].max()),
            },
            'largest_component_proportion_nodes': {
                'mean': float(df['largest_component_proportion_nodes'].mean()),
                'median': float(df['largest_component_proportion_nodes'].median()),
                'std': float(df['largest_component_proportion_nodes'].std()),
            },
            'largest_component_proportion_length': {
                'mean': float(df['largest_component_proportion_length'].mean()),
                'median': float(df['largest_component_proportion_length'].median()),
                'std': float(df['largest_component_proportion_length'].std()),
            },
        }
    }
    
    # Save outputs
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save per-image results
    summary_csv = output_dir / "summary.csv"
    df.to_csv(summary_csv, index=False)
    print(f"\n✓ Saved per-image results: {summary_csv}")
    
    # Save summary statistics
    stats_json = output_dir / "statistics.json"
    with open(stats_json, 'w') as f:
        json.dump(summary_stats, f, indent=2)
    print(f"✓ Saved summary statistics: {stats_json}")
    
    # Save failed images log if any
    if failed:
        failed_csv = output_dir / "failed.csv"
        pd.DataFrame(failed).to_csv(failed_csv, index=False)
        print(f"✓ Saved failed images log: {failed_csv}")
    
    # Print summary
    print(f"\n{'='*70}")
    print("SUMMARY STATISTICS")
    print("="*70)
    print(f"Connected Components:")
    print(f"  Mean: {summary_stats['connectivity_stats']['num_components']['mean']:.2f}")
    print(f"  Median: {summary_stats['connectivity_stats']['num_components']['median']:.2f}")
    print(f"  Std: {summary_stats['connectivity_stats']['num_components']['std']:.2f}")
    print(f"\nOptic Disc Connectivity (proportion of nodes):")
    print(f"  Mean: {summary_stats['connectivity_stats']['proportion_nodes_connected']['mean']:.2%}")
    print(f"  Median: {summary_stats['connectivity_stats']['proportion_nodes_connected']['median']:.2%}")
    print(f"\nOptic Disc Connectivity (proportion of length):")
    print(f"  Mean: {summary_stats['connectivity_stats']['proportion_length_connected']['mean']:.2%}")
    print(f"  Median: {summary_stats['connectivity_stats']['proportion_length_connected']['median']:.2%}")
    print("="*70)
    
    return summary_stats


def main():
    parser = argparse.ArgumentParser(
        description='Batch calculate connectivity metrics for entire VascX dataset'
    )
    
    parser.add_argument('dataset_path', type=Path, help='Path to dataset root')
    parser.add_argument('output_dir', type=Path, help='Output directory for results')
    parser.add_argument('--av-subfolder', default='av', help='Vessel masks subfolder')
    parser.add_argument('--fundus-subfolder', default='rgb', help='Fundus images subfolder')
    parser.add_argument('--batch-size', type=int, default=250,
                        help='Number of images to process per batch (default: 250)')
    
    args = parser.parse_args()
    
    process_dataset_batch(
        args.dataset_path,
        args.output_dir,
        args.av_subfolder,
        args.fundus_subfolder,
        args.batch_size
    )


if __name__ == '__main__':
    main()
