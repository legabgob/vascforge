# workflow/rules/connectivity_metrics_vascx.smk
"""
Calculate connectivity metrics from VascX datasets using BATCH processing.
Processes entire datasets at once instead of per-image to avoid job explosion.
"""

import re
from pathlib import Path

# Read config values
VASCX_VIEW_ROOT = config.get("vascx", {}).get("dataset_view_root", "results/vascx_datasets")
VASCX_UNREFINED_ROOT = config.get("vascx", {}).get("dataset_unrefined_root", "results/vascx_datasets_unrefined")
CONNECTIVITY_OUT = config.get("vascx", {}).get("connectivity_out", "results/connectivity_metrics")

# Expects: SIMPLE_DATASETS, OTHERDIR_DATASETS, RESOLUTIONS, K_VALUES, OTHERDIRS


# ==============================================================================
# REFINED DATASETS - Batch connectivity calculation
# ==============================================================================

rule connectivity_refined_simple:
    """Calculate connectivity for entire refined dataset (simple, no intermediate dir)."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, SIMPLE_DATASETS)) if SIMPLE_DATASETS else "NO_MATCH"
    input:
        ds_dir = f"{VASCX_VIEW_ROOT}" + "/{dataset}/k{k}/downsampled/{res}px"
    output:
        summary = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/summary.csv",
        stats = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/statistics.json"
    params:
        output_dir = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px"
    log:
        "logs/connectivity/refined/{dataset}_k{k}_{res}px.log"
    shell:
        """
        python workflow/scripts/calculate_connectivity_vascx_batch.py \
            {input.ds_dir} \
            {params.output_dir} \
            --av-subfolder av \
            --fundus-subfolder rgb \
            2>&1 | tee {log}
        """


rule connectivity_refined_otherdir:
    """Calculate connectivity for entire refined dataset (with intermediate directory)."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, OTHERDIR_DATASETS)) if OTHERDIR_DATASETS else "NO_MATCH"
    input:
        ds_dir = f"{VASCX_VIEW_ROOT}" + "/{dataset}/{other_dir}/k{k}/downsampled/{res}px"
    output:
        summary = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/summary.csv",
        stats = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/statistics.json"
    params:
        output_dir = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px"
    log:
        "logs/connectivity/refined/{dataset}_{other_dir}_k{k}_{res}px.log"
    shell:
        """
        python workflow/scripts/calculate_connectivity_vascx_batch.py \
            {input.ds_dir} \
            {params.output_dir} \
            --av-subfolder av \
            --fundus-subfolder rgb \
            2>&1 | tee {log}
        """


# ==============================================================================
# UNREFINED DATASETS - Batch connectivity calculation
# ==============================================================================

rule connectivity_unrefined_simple:
    """Calculate connectivity for entire unrefined dataset (simple, no intermediate dir)."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, SIMPLE_DATASETS)) if SIMPLE_DATASETS else "NO_MATCH"
    input:
        ds_dir = f"{VASCX_UNREFINED_ROOT}" + "/{dataset}/downsampled/{res}px"
    output:
        summary = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/downsampled/{res}px/summary.csv",
        stats = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/downsampled/{res}px/statistics.json"
    params:
        output_dir = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/downsampled/{res}px"
    log:
        "logs/connectivity/unrefined/{dataset}_{res}px.log"
    shell:
        """
        python workflow/scripts/calculate_connectivity_vascx_batch.py \
            {input.ds_dir} \
            {params.output_dir} \
            --av-subfolder av \
            --fundus-subfolder rgb \
            2>&1 | tee {log}
        """


rule connectivity_unrefined_otherdir:
    """Calculate connectivity for entire unrefined dataset (with intermediate directory)."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, OTHERDIR_DATASETS)) if OTHERDIR_DATASETS else "NO_MATCH"
    input:
        ds_dir = f"{VASCX_UNREFINED_ROOT}" + "/{dataset}/{other_dir}/downsampled/{res}px"
    output:
        summary = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/{other_dir}/downsampled/{res}px/summary.csv",
        stats = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/{other_dir}/downsampled/{res}px/statistics.json"
    params:
        output_dir = f"{CONNECTIVITY_OUT}" + "/unrefined/{dataset}/{other_dir}/downsampled/{res}px"
    log:
        "logs/connectivity/unrefined/{dataset}_{other_dir}_{res}px.log"
    shell:
        """
        python workflow/scripts/calculate_connectivity_vascx_batch.py \
            {input.ds_dir} \
            {params.output_dir} \
            --av-subfolder av \
            --fundus-subfolder rgb \
            2>&1 | tee {log}
        """


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
        stats = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/statistics.json"
    output:
        cc_dist = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/plots/cc_distribution.pdf",
        size_dist = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/plots/size_distribution.pdf",
        od_conn = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/k{k}/downsampled/{res}px/plots/od_connectivity.pdf"
    log:
        "logs/connectivity/plot_refined_{dataset}_k{k}_{res}px.log"
    script:
        "../scripts/plot_connectivity_metrics_batch.py"


rule plot_connectivity_refined_otherdir:
    """Create visualization plots for refined datasets with intermediate directory."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, OTHERDIR_DATASETS)) if OTHERDIR_DATASETS else "NO_MATCH"
    input:
        summary = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/summary.csv",
        stats = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/statistics.json"
    output:
        cc_dist = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/plots/cc_distribution.pdf",
        size_dist = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/plots/size_distribution.pdf",
        od_conn = f"{CONNECTIVITY_OUT}" + "/refined/{dataset}/{other_dir}/k{k}/downsampled/{res}px/plots/od_connectivity.pdf"
    log:
        "logs/connectivity/plot_refined_{dataset}_{other_dir}_k{k}_{res}px.log"
    script:
        "../scripts/plot_connectivity_metrics_batch.py"
