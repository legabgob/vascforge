# workflow/rules/connectivity_metrics_vascx.smk
"""
Calculate connectivity metrics from VascX datasets.
Works with COMBINED vessel masks (not separate arteries/veins).
No conda environment needed - uses existing environment.
"""

import re
from pathlib import Path

# Read config values (matching your existing structure)
VASCX_VIEW_ROOT = config.get("vascx", {}).get("dataset_view_root", "results/vascx_datasets")
VASCX_UNREFINED_ROOT = config.get("vascx", {}).get("dataset_unrefined_root", "results/vascx_datasets_unrefined")
CONNECTIVITY_OUT = config.get("vascx", {}).get("connectivity_out", "results/connectivity_metrics")
N_JOBS = int(config.get("vascx", {}).get("n_jobs", 64))

# Expects these from your dataset discovery rules:
# SIMPLE_DATASETS, OTHERDIR_DATASETS, RESOLUTIONS, K_VALUES, OTHERDIRS


# ==============================================================================
# Helper function to get image IDs from VascX dataset
# ==============================================================================

def get_image_ids_from_vascx(ds_path):
    """Get list of image IDs from a VascX dataset directory."""
    ds_path = Path(ds_path)
    av_dir = ds_path / "av"
    
    if not av_dir.exists():
        return []
    
    # Get all .png files and extract IDs (without extension)
    image_ids = sorted([f.stem for f in av_dir.glob("*.png")])
    return image_ids


# ==============================================================================
# REFINED DATASETS - Calculate connectivity metrics
# ==============================================================================

rule connectivity_refined_simple:
    """Calculate connectivity for simple refined datasets (no intermediate dir)."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, SIMPLE_DATASETS)) if SIMPLE_DATASETS else "NO_MATCH"
    input:
        ds_dir = f"{VASCX_VIEW_ROOT}" + "/{dataset}/k{k}/downsampled/{res}px"
    output:
        json = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/{image_id}_metrics.json",
        csv = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/{image_id}_components.csv"
    params:
        output_prefix = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/{image_id}"
    log:
        "logs/connectivity/refined/{dataset}_k{k}_{res}px_{image_id}.log"
    shell:
        """
        python workflow/scripts/calculate_connectivity_vascx.py dataset \
            {input.ds_dir} \
            {wildcards.image_id} \
            {params.output_prefix} \
            --av-subfolder av \
            --fundus-subfolder rgb \
            2>&1 | tee {log}
        """


rule connectivity_refined_otherdir:
    """Calculate connectivity for refined datasets with intermediate directory."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, OTHERDIR_DATASETS)) if OTHERDIR_DATASETS else "NO_MATCH"
    input:
        ds_dir = f"{VASCX_VIEW_ROOT}" + "/{dataset}/{other_dir}/k{k}/downsampled/{res}px"
    output:
        json = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/{image_id}_metrics.json",
        csv = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/{image_id}_components.csv"
    params:
        output_prefix = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/{image_id}"
    log:
        "logs/connectivity/refined/{dataset}_{other_dir}_k{k}_{res}px_{image_id}.log"
    shell:
        """
        python workflow/scripts/calculate_connectivity_vascx.py dataset \
            {input.ds_dir} \
            {wildcards.image_id} \
            {params.output_prefix} \
            --av-subfolder av \
            --fundus-subfolder rgb \
            2>&1 | tee {log}
        """


# ==============================================================================
# UNREFINED DATASETS - Calculate connectivity metrics
# ==============================================================================

rule connectivity_unrefined_simple:
    """Calculate connectivity for simple unrefined datasets (no intermediate dir)."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, SIMPLE_DATASETS)) if SIMPLE_DATASETS else "NO_MATCH"
    input:
        ds_dir = f"{VASCX_UNREFINED_ROOT}" + "/{dataset}/downsampled/{res}px"
    output:
        json = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/downsampled/{res}px/{image_id}_metrics.json",
        csv = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/downsampled/{res}px/{image_id}_components.csv"
    params:
        output_prefix = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/downsampled/{res}px/{image_id}"
    log:
        "logs/connectivity/unrefined/{dataset}_{res}px_{image_id}.log"
    shell:
        """
        python workflow/scripts/calculate_connectivity_vascx.py dataset \
            {input.ds_dir} \
            {wildcards.image_id} \
            {params.output_prefix} \
            --av-subfolder av \
            --fundus-subfolder rgb \
            2>&1 | tee {log}
        """


rule connectivity_unrefined_otherdir:
    """Calculate connectivity for unrefined datasets with intermediate directory."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, OTHERDIR_DATASETS)) if OTHERDIR_DATASETS else "NO_MATCH"
    input:
        ds_dir = f"{VASCX_UNREFINED_ROOT}" + "/{dataset}/{other_dir}/downsampled/{res}px"
    output:
        json = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/{other_dir}/downsampled/{res}px/{image_id}_metrics.json",
        csv = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/{other_dir}/downsampled/{res}px/{image_id}_components.csv"
    params:
        output_prefix = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/{other_dir}/downsampled/{res}px/{image_id}"
    log:
        "logs/connectivity/unrefined/{dataset}_{other_dir}_{res}px_{image_id}.log"
    shell:
        """
        python workflow/scripts/calculate_connectivity_vascx.py dataset \
            {input.ds_dir} \
            {wildcards.image_id} \
            {params.output_prefix} \
            --av-subfolder av \
            --fundus-subfolder rgb \
            2>&1 | tee {log}
        """


# ==============================================================================
# AGGREGATE - Combine metrics across all images
# ==============================================================================

rule aggregate_connectivity_refined_simple:
    """Aggregate connectivity metrics for simple refined datasets."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, SIMPLE_DATASETS)) if SIMPLE_DATASETS else "NO_MATCH"
    input:
        json_files = lambda wc: expand(
            f"{CONNECTIVITY_OUT}/refined/{{dataset}}/k{{k}}/downsampled/{{res}}px/{{image_id}}_metrics.json",
            dataset=wc.dataset,
            k=wc.k,
            res=wc.res,
            image_id=get_image_ids_from_vascx(f"{VASCX_VIEW_ROOT}/{wc.dataset}/k{wc.k}/downsampled/{wc.res}px")
        ),
        csv_files = lambda wc: expand(
            f"{CONNECTIVITY_OUT}/refined/{{dataset}}/k{{k}}/downsampled/{{res}}px/{{image_id}}_components.csv",
            dataset=wc.dataset,
            k=wc.k,
            res=wc.res,
            image_id=get_image_ids_from_vascx(f"{VASCX_VIEW_ROOT}/{wc.dataset}/k{wc.k}/downsampled/{wc.res}px")
        )
    output:
        summary = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/summary.csv",
        components = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/all_components.csv",
        stats = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/statistics.json"
    log:
        "logs/connectivity/aggregate_refined_{dataset}_k{k}_{res}px.log"
    script:
        "../scripts/aggregate_connectivity_metrics.py"


rule aggregate_connectivity_refined_otherdir:
    """Aggregate connectivity metrics for refined datasets with intermediate directory."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, OTHERDIR_DATASETS)) if OTHERDIR_DATASETS else "NO_MATCH"
    input:
        json_files = lambda wc: expand(
            f"{CONNECTIVITY_OUT}/refined/{{dataset}}/{{other_dir}}/k{{k}}/downsampled/{{res}}px/{{image_id}}_metrics.json",
            dataset=wc.dataset,
            other_dir=wc.other_dir,
            k=wc.k,
            res=wc.res,
            image_id=get_image_ids_from_vascx(f"{VASCX_VIEW_ROOT}/{wc.dataset}/{wc.other_dir}/k{wc.k}/downsampled/{wc.res}px")
        ),
        csv_files = lambda wc: expand(
            f"{CONNECTIVITY_OUT}/refined/{{dataset}}/{{other_dir}}/k{{k}}/downsampled/{{res}}px/{{image_id}}_components.csv",
            dataset=wc.dataset,
            other_dir=wc.other_dir,
            k=wc.k,
            res=wc.res,
            image_id=get_image_ids_from_vascx(f"{VASCX_VIEW_ROOT}/{wc.dataset}/{wc.other_dir}/k{wc.k}/downsampled/{wc.res}px")
        )
    output:
        summary = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/summary.csv",
        components = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/all_components.csv",
        stats = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/statistics.json"
    log:
        "logs/connectivity/aggregate_refined_{dataset}_{other_dir}_k{k}_{res}px.log"
    script:
        "../scripts/aggregate_connectivity_metrics.py"


rule aggregate_connectivity_unrefined_simple:
    """Aggregate connectivity metrics for simple unrefined datasets."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, SIMPLE_DATASETS)) if SIMPLE_DATASETS else "NO_MATCH"
    input:
        json_files = lambda wc: expand(
            f"{CONNECTIVITY_OUT}/unrefined/{{dataset}}/downsampled/{{res}}px/{{image_id}}_metrics.json",
            dataset=wc.dataset,
            res=wc.res,
            image_id=get_image_ids_from_vascx(f"{VASCX_UNREFINED_ROOT}/{wc.dataset}/downsampled/{wc.res}px")
        ),
        csv_files = lambda wc: expand(
            f"{CONNECTIVITY_OUT}/unrefined/{{dataset}}/downsampled/{{res}}px/{{image_id}}_components.csv",
            dataset=wc.dataset,
            res=wc.res,
            image_id=get_image_ids_from_vascx(f"{VASCX_UNREFINED_ROOT}/{wc.dataset}/downsampled/{wc.res}px")
        )
    output:
        summary = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/downsampled/{res}px/summary.csv",
        components = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/downsampled/{res}px/all_components.csv",
        stats = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/downsampled/{res}px/statistics.json"
    log:
        "logs/connectivity/aggregate_unrefined_{dataset}_{res}px.log"
    script:
        "../scripts/aggregate_connectivity_metrics.py"


rule aggregate_connectivity_unrefined_otherdir:
    """Aggregate connectivity metrics for unrefined datasets with intermediate directory."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, OTHERDIR_DATASETS)) if OTHERDIR_DATASETS else "NO_MATCH"
    input:
        json_files = lambda wc: expand(
            f"{CONNECTIVITY_OUT}/unrefined/{{dataset}}/{{other_dir}}/downsampled/{{res}}px/{{image_id}}_metrics.json",
            dataset=wc.dataset,
            other_dir=wc.other_dir,
            res=wc.res,
            image_id=get_image_ids_from_vascx(f"{VASCX_UNREFINED_ROOT}/{wc.dataset}/{wc.other_dir}/downsampled/{wc.res}px")
        ),
        csv_files = lambda wc: expand(
            f"{CONNECTIVITY_OUT}/unrefined/{{dataset}}/{{other_dir}}/downsampled/{{res}}px/{{image_id}}_components.csv",
            dataset=wc.dataset,
            other_dir=wc.other_dir,
            res=wc.res,
            image_id=get_image_ids_from_vascx(f"{VASCX_UNREFINED_ROOT}/{wc.dataset}/{wc.other_dir}/downsampled/{wc.res}px")
        )
    output:
        summary = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/{other_dir}/downsampled/{res}px/summary.csv",
        components = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/{other_dir}/downsampled/{res}px/all_components.csv",
        stats = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/{other_dir}/downsampled/{res}px/statistics.json"
    log:
        "logs/connectivity/aggregate_unrefined_{dataset}_{other_dir}_{res}px.log"
    script:
        "../scripts/aggregate_connectivity_metrics.py"


# ==============================================================================
# COMPARE - Before/after refinement comparison
# ==============================================================================

rule compare_connectivity_simple:
    """Compare unrefined vs refined connectivity for simple datasets."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, SIMPLE_DATASETS)) if SIMPLE_DATASETS else "NO_MATCH"
    input:
        unrefined = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/downsampled/{res}px/summary.csv",
        refined = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/summary.csv"
    output:
        comparison = f"{CONNECTIVITY_OUT}" + "/comparison/{dataset}/k{k}/downsampled/{res}px/refinement_comparison.csv",
        stats = f"{CONNECTIVITY_OUT}" + "/comparison/{dataset}/k{k}/downsampled/{res}px/refinement_stats.json",
        plot = f"{CONNECTIVITY_OUT}" + "/comparison/{dataset}/k{k}/downsampled/{res}px/refinement_plot.pdf"
    params:
        metrics = ['num_components', 'proportion_nodes_connected', 'proportion_length_connected']
    log:
        "logs/connectivity/compare_{dataset}_k{k}_{res}px.log"
    script:
        "../scripts/compare_connectivity_refinement.py"


rule compare_connectivity_otherdir:
    """Compare unrefined vs refined connectivity for datasets with intermediate directory."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, OTHERDIR_DATASETS)) if OTHERDIR_DATASETS else "NO_MATCH"
    input:
        unrefined = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/{other_dir}/downsampled/{res}px/summary.csv",
        refined = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/summary.csv"
    output:
        comparison = f"{CONNECTIVITY_OUT}" + "/comparison/{dataset}/{other_dir}/k{k}/downsampled/{res}px/refinement_comparison.csv",
        stats = f"{CONNECTIVITY_OUT}" + "/comparison/{dataset}/{other_dir}/k{k}/downsampled/{res}px/refinement_stats.json",
        plot = f"{CONNECTIVITY_OUT}" + "/comparison/{dataset}/{other_dir}/k{k}/downsampled/{res}px/refinement_plot.pdf"
    params:
        metrics = ['num_components', 'proportion_nodes_connected', 'proportion_length_connected']
    log:
        "logs/connectivity/compare_{dataset}_{other_dir}_k{k}_{res}px.log"
    script:
        "../scripts/compare_connectivity_refinement.py"


# ==============================================================================
# VISUALIZE - Create plots
# ==============================================================================

rule plot_connectivity_refined_simple:
    """Create visualization plots for refined simple datasets."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, SIMPLE_DATASETS)) if SIMPLE_DATASETS else "NO_MATCH"
    input:
        summary = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/summary.csv",
        components = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/all_components.csv"
    output:
        cc_dist = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/plots/cc_distribution.pdf",
        size_dist = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/plots/size_distribution.pdf",
        od_conn = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/plots/od_connectivity.pdf"
    log:
        "logs/connectivity/plot_refined_{dataset}_k{k}_{res}px.log"
    script:
        "../scripts/plot_connectivity_metrics.py"


rule plot_connectivity_refined_otherdir:
    """Create visualization plots for refined datasets with intermediate directory."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, OTHERDIR_DATASETS)) if OTHERDIR_DATASETS else "NO_MATCH"
    input:
        summary = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/summary.csv",
        components = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/all_components.csv"
    output:
        cc_dist = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/plots/cc_distribution.pdf",
        size_dist = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/plots/size_distribution.pdf",
        od_conn = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/plots/od_connectivity.pdf"
    log:
        "logs/connectivity/plot_refined_{dataset}_{other_dir}_k{k}_{res}px.log"
    script:
        "../scripts/plot_connectivity_metrics.py"
