#!/usr/bin/env python3
"""
Calculate connectivity metrics from binary segmentation masks using VascX.

This script:
1. Loads binary .png masks using VascX RetinaLoader
2. Accesses the automatically created vascular graphs
3. Calculates connectivity metrics:
   - Number of connected components
   - Component sizes (nodes and length)
   - Proportion connected to optic disc
"""

import argparse
import json
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import warnings
warnings.filterwarnings('ignore')

# VascX imports
from vascx.fundus.loader import RetinaLoader
from vascx.fundus.layer import VesselTreeLayer


# ============================================================================
# Graph utilities adapted from VascX tutorial
# ============================================================================

def get_seg_ids_from_cc(graph: Union[nx.Graph, nx.DiGraph]) -> List[frozenset]:
    """Get list of segment IDs from a connected component."""
    list_seg_ids = []
    for s, e in graph.edges():
        list_seg_ids.append(frozenset([s, e]))
    return list_seg_ids


def get_segments_from_cc(graph: Union[nx.Graph, nx.DiGraph], layer: VesselTreeLayer) -> List:
    """Get list of segments from a connected component."""
    seg_list = []
    list_seg_ids = get_seg_ids_from_cc(graph)
    
    ids = [seg.id for seg in layer.segments]
    segments = layer.segments
    dict_seg = dict(zip(ids, segments))
    
    for seg_id in list_seg_ids:
        if seg_id in dict_seg:
            seg_list.append(dict_seg[seg_id])
    
    return seg_list


def get_graph_total_length(graph: Union[nx.Graph, nx.DiGraph], layer: VesselTreeLayer) -> float:
    """Calculate total length of segments in a connected component."""
    seg_list = get_segments_from_cc(graph, layer)
    return sum([seg.length for seg in seg_list])


def get_connected_components(layer: VesselTreeLayer, sort: bool = True) -> List[Union[nx.Graph, nx.DiGraph]]:
    """Get connected components from a vessel layer."""
    graph = layer.digraph  # Use directed graph
    
    if nx.is_directed(graph):
        components = nx.weakly_connected_components(graph)
    else:
        components = nx.connected_components(graph)
    
    cc_list = [graph.subgraph(c).copy() for c in components]
    
    if sort:
        # Sort by total length (descending)
        return sorted(cc_list, key=lambda x: get_graph_total_length(x, layer), reverse=True)
    else:
        return cc_list


def cc_out_of_optic_disc(layer: VesselTreeLayer, cc: Union[nx.Graph, nx.DiGraph]) -> bool:
    """Check if connected component touches the optic disc."""
    pts = []
    for node_idx in cc.nodes:
        point = (cc.nodes[node_idx]['o'][0], cc.nodes[node_idx]['o'][1])
        pts.append(point)
    
    return layer.disc.intersects(pts, with_boundary=True)


def get_cc_subset_touching_od(layer: VesselTreeLayer, cc_list: List[Union[nx.Graph, nx.DiGraph]]) -> List[Union[nx.Graph, nx.DiGraph]]:
    """Get connected components that touch the optic disc."""
    return [cc for cc in cc_list if cc_out_of_optic_disc(layer, cc)]


# ============================================================================
# Connectivity metrics calculation
# ============================================================================

def calculate_connected_components_metrics(layer: VesselTreeLayer) -> Dict:
    """
    Calculate detailed metrics about connected components.
    
    Returns dictionary with:
    - num_components: Total number of connected components
    - component_sizes_nodes: List of component sizes (number of nodes)
    - component_sizes_length: List of component sizes (total segment length)
    - largest_component_nodes: Size of largest component (nodes)
    - largest_component_length: Size of largest component (length)
    """
    # Get connected components
    cc_list = get_connected_components(layer, sort=True)
    num_components = len(cc_list)
    
    component_sizes_nodes = []
    component_sizes_length = []
    
    for cc in cc_list:
        # Number of nodes
        num_nodes = cc.number_of_nodes()
        component_sizes_nodes.append(num_nodes)
        
        # Total segment length
        total_length = get_graph_total_length(cc, layer)
        component_sizes_length.append(total_length)
    
    # Calculate total graph stats
    total_nodes = layer.graph.number_of_nodes()
    total_length = sum(component_sizes_length)
    
    return {
        'num_components': num_components,
        'component_sizes_nodes': component_sizes_nodes,
        'component_sizes_length': component_sizes_length,
        'largest_component_nodes': component_sizes_nodes[0] if component_sizes_nodes else 0,
        'largest_component_length': component_sizes_length[0] if component_sizes_length else 0,
        'total_nodes': total_nodes,
        'total_edges': layer.graph.number_of_edges(),
        'total_length': total_length,
    }


def calculate_optic_disc_connectivity(layer: VesselTreeLayer) -> Dict:
    """
    Calculate proportion of graph connected to the optic disc.
    
    Returns dictionary with:
    - num_components_touching_od: Number of components touching optic disc
    - connected_nodes: Number of nodes in components touching optic disc
    - connected_length: Total length in components touching optic disc
    - proportion_nodes: Proportion of nodes connected to optic disc
    - proportion_length: Proportion of total length connected to optic disc
    """
    # Get all connected components
    cc_list = get_connected_components(layer, sort=True)
    
    # Get components touching optic disc
    cc_touching_od = get_cc_subset_touching_od(layer, cc_list)
    
    # Calculate metrics for components touching OD
    connected_nodes = sum(cc.number_of_nodes() for cc in cc_touching_od)
    connected_length = sum(get_graph_total_length(cc, layer) for cc in cc_touching_od)
    
    # Total graph stats
    total_nodes = layer.graph.number_of_nodes()
    total_length = sum(get_graph_total_length(cc, layer) for cc in cc_list)
    
    return {
        'num_components_touching_od': len(cc_touching_od),
        'connected_nodes': connected_nodes,
        'connected_length': connected_length,
        'proportion_nodes': connected_nodes / total_nodes if total_nodes > 0 else 0,
        'proportion_length': connected_length / total_length if total_length > 0 else 0,
        'total_nodes': total_nodes,
        'total_edges': layer.graph.number_of_edges(),
        'total_length': total_length,
    }


def calculate_component_size_distribution(layer: VesselTreeLayer) -> pd.DataFrame:
    """
    Calculate detailed size distribution of connected components.
    
    Returns DataFrame with one row per component.
    """
    cc_list = get_connected_components(layer, sort=True)
    
    component_data = []
    for i, cc in enumerate(cc_list, 1):
        num_nodes = cc.number_of_nodes()
        num_edges = cc.number_of_edges()
        total_length = get_graph_total_length(cc, layer)
        
        # Check if touches optic disc
        touches_od = cc_out_of_optic_disc(layer, cc)
        
        component_data.append({
            'component_id': i,
            'num_nodes': num_nodes,
            'num_edges': num_edges,
            'total_length': total_length,
            'touches_optic_disc': touches_od,
        })
    
    df = pd.DataFrame(component_data)
    
    # Add relative sizes
    if len(df) > 0:
        df['relative_nodes'] = df['num_nodes'] / df['num_nodes'].sum()
        df['relative_length'] = df['total_length'] / df['total_length'].sum()
        
        # Rank by size
        df['rank_by_length'] = range(1, len(df) + 1)
        df['rank_by_nodes'] = df['num_nodes'].rank(ascending=False, method='dense').astype(int)
    
    return df


# ============================================================================
# Main processing
# ============================================================================

def process_single_image(
    mask_path: Path,
    vessel_type: str,
    output_prefix: Path
) -> Dict:
    """
    Process a single binary mask using VascX.
    
    Args:
        mask_path: Path to binary .png mask
        vessel_type: 'arteries' or 'veins'
        output_prefix: Prefix for output files
    
    Returns:
        Dictionary with all calculated metrics
    """
    print(f"Processing {mask_path} as {vessel_type}...")
    
    # Load mask using VascX
    # Note: VascX expects a folder structure, so we need to adapt
    # For now, we'll load the mask directly and create a minimal Retina-like structure
    # This is a simplified version - you may need to adjust based on your exact setup
    
    import cv2
    from vascx.fundus.retina import Retina
    
    # Load binary mask
    mask = cv2.imread(str(mask_path), cv2.IMREAD_GRAYSCALE)
    if mask is None:
        raise ValueError(f"Could not load mask from {mask_path}")
    
    # Convert to binary
    mask = (mask > 127).astype(np.uint8) * 255
    
    # Create a minimal Retina object
    # This is simplified - in practice, VascX handles this automatically
    retina = Retina.from_segmentation(mask, image_id=mask_path.stem)
    
    # Get the appropriate vessel layer
    if vessel_type == 'arteries':
        layer = retina.arteries
    elif vessel_type == 'veins':
        layer = retina.veins
    else:
        raise ValueError(f"Unknown vessel type: {vessel_type}")
    
    # Calculate metrics
    print("Calculating connected components metrics...")
    cc_metrics = calculate_connected_components_metrics(layer)
    
    print("Calculating optic disc connectivity...")
    od_metrics = calculate_optic_disc_connectivity(layer)
    
    print("Calculating component size distribution...")
    component_df = calculate_component_size_distribution(layer)
    
    # Combine metrics
    all_metrics = {
        'image_info': {
            'image_id': mask_path.stem,
            'vessel_type': vessel_type,
            'mask_path': str(mask_path),
        },
        'graph_info': {
            'total_nodes': cc_metrics['total_nodes'],
            'total_edges': cc_metrics['total_edges'],
            'total_length': cc_metrics['total_length'],
        },
        'connected_components': {
            'num_components': cc_metrics['num_components'],
            'largest_component_nodes': cc_metrics['largest_component_nodes'],
            'largest_component_length': cc_metrics['largest_component_length'],
            'component_sizes_nodes': cc_metrics['component_sizes_nodes'],
            'component_sizes_length': cc_metrics['component_sizes_length'],
        },
        'optic_disc_connectivity': {
            'num_components_touching_od': od_metrics['num_components_touching_od'],
            'connected_nodes': od_metrics['connected_nodes'],
            'connected_length': od_metrics['connected_length'],
            'proportion_nodes': od_metrics['proportion_nodes'],
            'proportion_length': od_metrics['proportion_length'],
        },
        'summary_stats': {
            'mean_component_size_nodes': float(np.mean(cc_metrics['component_sizes_nodes'])),
            'median_component_size_nodes': float(np.median(cc_metrics['component_sizes_nodes'])),
            'std_component_size_nodes': float(np.std(cc_metrics['component_sizes_nodes'])),
            'mean_component_size_length': float(np.mean(cc_metrics['component_sizes_length'])),
            'median_component_size_length': float(np.median(cc_metrics['component_sizes_length'])),
            'std_component_size_length': float(np.std(cc_metrics['component_sizes_length'])),
        }
    }
    
    # Save outputs
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    
    # Save JSON metrics
    json_output = output_prefix.with_suffix('.json')
    with open(json_output, 'w') as f:
        json.dump(all_metrics, f, indent=2)
    print(f"Saved metrics to {json_output}")
    
    # Save CSV with component details
    csv_output = output_prefix.with_suffix('.csv')
    component_df.to_csv(csv_output, index=False)
    print(f"Saved component distribution to {csv_output}")
    
    return all_metrics


def process_from_vascx_loader(
    dataset_path: Path,
    image_id: str,
    output_prefix: Path,
    av_subfolder: str = 'av',
    fundus_subfolder: str = 'rgb'
) -> Dict:
    """
    Process using VascX RetinaLoader (for properly structured datasets).
    Works with COMBINED vessel masks.
    
    Args:
        dataset_path: Path to dataset root
        image_id: Image identifier
        output_prefix: Prefix for output files
        av_subfolder: Subfolder containing vessel masks
        fundus_subfolder: Subfolder containing fundus images
    
    Returns:
        Dictionary with all calculated metrics
    """
    print(f"Loading dataset from {dataset_path}...")
    
    # Load using VascX
    loader = RetinaLoader.from_folder(
        dataset_path,
        av_subfolder=av_subfolder,
        fundus_subfolder=fundus_subfolder
    )
    
    # Find the retina by image_id
    retina = None
    for r in loader:
        if r.id == image_id:
            retina = r
            break
    
    if retina is None:
        raise ValueError(f"Image {image_id} not found in dataset")
    
    # Use the COMBINED vessels layer (arteries + veins together)
    # VascX creates a combined layer automatically
    layer = retina.vessels  # Use .vessels for combined, not .arteries or .veins
    
    # Calculate metrics
    cc_metrics = calculate_connected_components_metrics(layer)
    od_metrics = calculate_optic_disc_connectivity(layer)
    component_df = calculate_component_size_distribution(layer)
    
    # Combine and save
    all_metrics = {
        'image_info': {
            'image_id': image_id,
            'dataset_path': str(dataset_path),
        },
        'graph_info': {
            'total_nodes': cc_metrics['total_nodes'],
            'total_edges': cc_metrics['total_edges'],
            'total_length': cc_metrics['total_length'],
        },
        'connected_components': {
            'num_components': cc_metrics['num_components'],
            'largest_component_nodes': cc_metrics['largest_component_nodes'],
            'largest_component_length': cc_metrics['largest_component_length'],
            'component_sizes_nodes': cc_metrics['component_sizes_nodes'],
            'component_sizes_length': cc_metrics['component_sizes_length'],
        },
        'optic_disc_connectivity': {
            'num_components_touching_od': od_metrics['num_components_touching_od'],
            'connected_nodes': od_metrics['connected_nodes'],
            'connected_length': od_metrics['connected_length'],
            'proportion_nodes': od_metrics['proportion_nodes'],
            'proportion_length': od_metrics['proportion_length'],
        },
        'summary_stats': {
            'mean_component_size_nodes': float(np.mean(cc_metrics['component_sizes_nodes'])) if cc_metrics['component_sizes_nodes'] else 0,
            'median_component_size_nodes': float(np.median(cc_metrics['component_sizes_nodes'])) if cc_metrics['component_sizes_nodes'] else 0,
            'std_component_size_nodes': float(np.std(cc_metrics['component_sizes_nodes'])) if cc_metrics['component_sizes_nodes'] else 0,
            'mean_component_size_length': float(np.mean(cc_metrics['component_sizes_length'])) if cc_metrics['component_sizes_length'] else 0,
            'median_component_size_length': float(np.median(cc_metrics['component_sizes_length'])) if cc_metrics['component_sizes_length'] else 0,
            'std_component_size_length': float(np.std(cc_metrics['component_sizes_length'])) if cc_metrics['component_sizes_length'] else 0,
        }
    }
    
    # Save outputs
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    
    json_output = output_prefix.with_suffix('.json')
    with open(json_output, 'w') as f:
        json.dump(all_metrics, f, indent=2)
    
    csv_output = output_prefix.with_suffix('.csv')
    component_df.to_csv(csv_output, index=False)
    
    print(f"\n{'='*60}")
    print("CONNECTIVITY METRICS SUMMARY")
    print("="*60)
    print(f"Image: {image_id}")
    print(f"Total nodes: {all_metrics['graph_info']['total_nodes']}")
    print(f"Total edges: {all_metrics['graph_info']['total_edges']}")
    print(f"\nConnected Components: {cc_metrics['num_components']}")
    if cc_metrics['total_nodes'] > 0:
        print(f"  Largest (nodes): {cc_metrics['largest_component_nodes']} ({cc_metrics['largest_component_nodes']/cc_metrics['total_nodes']*100:.1f}%)")
    if cc_metrics['total_length'] > 0:
        print(f"  Largest (length): {cc_metrics['largest_component_length']:.1f} ({cc_metrics['largest_component_length']/cc_metrics['total_length']*100:.1f}%)")
    print(f"\nOptic Disc Connectivity:")
    print(f"  Components touching OD: {od_metrics['num_components_touching_od']}")
    print(f"  Nodes connected: {od_metrics['connected_nodes']} ({od_metrics['proportion_nodes']*100:.1f}%)")
    print(f"  Length connected: {od_metrics['connected_length']:.1f} ({od_metrics['proportion_length']*100:.1f}%)")
    print("="*60)
    
    return all_metrics


def main():
    parser = argparse.ArgumentParser(
        description='Calculate connectivity metrics from binary masks using VascX'
    )
    
    subparsers = parser.add_subparsers(dest='mode', help='Processing mode')
    
    # Dataset mode (using VascX loader)
    dataset_parser = subparsers.add_parser('dataset', help='Process from VascX-structured dataset')
    dataset_parser.add_argument('dataset_path', type=Path, help='Path to dataset root')
    dataset_parser.add_argument('image_id', type=str, help='Image identifier')
    dataset_parser.add_argument('output_prefix', type=Path, help='Prefix for output files')
    dataset_parser.add_argument('--av-subfolder', default='av', help='Vessel masks subfolder')
    dataset_parser.add_argument('--fundus-subfolder', default='rgb', help='Fundus images subfolder')
    
    args = parser.parse_args()
    
    if args.mode == 'dataset':
        process_from_vascx_loader(
            args.dataset_path,
            args.image_id,
            args.output_prefix,
            args.av_subfolder,
            args.fundus_subfolder
        )
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
