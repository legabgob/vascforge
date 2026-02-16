#!/usr/bin/env python3
"""
Compare connectivity metrics before and after Greco refinement.

Performs statistical tests and generates comparison summaries.
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns


def load_and_merge_data(unrefined_path: Path, refined_path: Path) -> pd.DataFrame:
    """
    Load and merge unrefined and refined data.
    
    Returns DataFrame with matched pairs for statistical testing.
    """
    unrefined = pd.read_csv(unrefined_path)
    refined = pd.read_csv(refined_path)
    
    # Merge on image_id
    merged = unrefined.merge(
        refined,
        on='image_id',
        suffixes=('_unrefined', '_refined')
    )
    
    return merged


def calculate_changes(df: pd.DataFrame, metrics: list) -> pd.DataFrame:
    """
    Calculate changes (delta and percent change) for each metric.
    """
    for metric in metrics:
        unrefined_col = f"{metric}_unrefined"
        refined_col = f"{metric}_refined"
        
        if unrefined_col in df.columns and refined_col in df.columns:
            # Absolute change
            df[f"{metric}_delta"] = df[refined_col] - df[unrefined_col]
            
            # Percent change
            df[f"{metric}_pct_change"] = (
                (df[refined_col] - df[unrefined_col]) / 
                df[unrefined_col].replace(0, np.nan) * 100
            )
    
    return df


def perform_statistical_tests(df: pd.DataFrame, metrics: list) -> dict:
    """
    Perform paired statistical tests for each metric.
    
    Tests:
    - Paired t-test (parametric)
    - Wilcoxon signed-rank test (non-parametric)
    - Effect size (Cohen's d)
    """
    results = {}
    
    for metric in metrics:
        unrefined_col = f"{metric}_unrefined"
        refined_col = f"{metric}_refined"
        
        if unrefined_col not in df.columns or refined_col not in df.columns:
            continue
        
        unrefined = df[unrefined_col].dropna()
        refined = df[refined_col].dropna()
        
        # Ensure same length for paired tests
        common_idx = df[[unrefined_col, refined_col]].dropna().index
        unrefined = df.loc[common_idx, unrefined_col]
        refined = df.loc[common_idx, refined_col]
        
        if len(unrefined) < 2:
            continue
        
        # Paired t-test
        t_stat, t_pval = stats.ttest_rel(unrefined, refined)
        
        # Wilcoxon signed-rank test
        w_stat, w_pval = stats.wilcoxon(unrefined, refined)
        
        # Effect size (Cohen's d for paired samples)
        diff = refined - unrefined
        cohens_d = diff.mean() / diff.std()
        
        # Summary statistics
        results[metric] = {
            'n_pairs': len(unrefined),
            'unrefined_mean': float(unrefined.mean()),
            'unrefined_std': float(unrefined.std()),
            'refined_mean': float(refined.mean()),
            'refined_std': float(refined.std()),
            'mean_delta': float(diff.mean()),
            'std_delta': float(diff.std()),
            'median_delta': float(diff.median()),
            't_statistic': float(t_stat),
            't_pvalue': float(t_pval),
            'wilcoxon_statistic': float(w_stat),
            'wilcoxon_pvalue': float(w_pval),
            'cohens_d': float(cohens_d),
            'significant_t': t_pval < 0.05,
            'significant_w': w_pval < 0.05,
        }
        
        # Interpretation of effect size
        abs_d = abs(cohens_d)
        if abs_d < 0.2:
            effect_size_interp = 'negligible'
        elif abs_d < 0.5:
            effect_size_interp = 'small'
        elif abs_d < 0.8:
            effect_size_interp = 'medium'
        else:
            effect_size_interp = 'large'
        
        results[metric]['effect_size_interpretation'] = effect_size_interp
    
    return results


def create_comparison_plot(df: pd.DataFrame, metrics: list, output_path: Path):
    """
    Create visualization comparing unrefined vs refined metrics.
    """
    n_metrics = len(metrics)
    fig, axes = plt.subplots(n_metrics, 2, figsize=(12, 4*n_metrics))
    
    if n_metrics == 1:
        axes = axes.reshape(1, -1)
    
    for i, metric in enumerate(metrics):
        unrefined_col = f"{metric}_unrefined"
        refined_col = f"{metric}_refined"
        delta_col = f"{metric}_delta"
        
        if unrefined_col not in df.columns or refined_col not in df.columns:
            continue
        
        # Left plot: Before vs After scatter
        ax1 = axes[i, 0]
        ax1.scatter(df[unrefined_col], df[refined_col], alpha=0.5)
        
        # Add y=x line
        min_val = min(df[unrefined_col].min(), df[refined_col].min())
        max_val = max(df[unrefined_col].max(), df[refined_col].max())
        ax1.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, label='y=x')
        
        ax1.set_xlabel(f'{metric} (Unrefined)', fontsize=10)
        ax1.set_ylabel(f'{metric} (Refined)', fontsize=10)
        ax1.set_title(f'{metric}: Before vs After Refinement', fontsize=11)
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Right plot: Distribution of changes
        ax2 = axes[i, 1]
        
        if delta_col in df.columns:
            df[delta_col].hist(bins=30, ax=ax2, edgecolor='black')
            ax2.axvline(0, color='r', linestyle='--', linewidth=2, label='No change')
            ax2.axvline(df[delta_col].mean(), color='g', linestyle='-', 
                       linewidth=2, label=f'Mean: {df[delta_col].mean():.3f}')
            
            ax2.set_xlabel(f'Change in {metric}', fontsize=10)
            ax2.set_ylabel('Frequency', fontsize=10)
            ax2.set_title(f'Distribution of Changes in {metric}', fontsize=11)
            ax2.legend()
            ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def main(snakemake):
    """Main function for Snakemake."""
    
    # Load data
    print("Loading data...")
    df = load_and_merge_data(
        Path(snakemake.input.unrefined),
        Path(snakemake.input.refined)
    )
    
    print(f"Loaded {len(df)} matched image pairs")
    
    # Get metrics to compare
    metrics = snakemake.params.get('metrics', [
        'num_components',
        'proportion_nodes_connected',
        'proportion_length_connected'
    ])
    
    # Calculate changes
    print("Calculating changes...")
    df = calculate_changes(df, metrics)
    
    # Save comparison dataframe
    df.to_csv(snakemake.output.comparison, index=False)
    print(f"Saved comparison to {snakemake.output.comparison}")
    
    # Perform statistical tests
    print("Performing statistical tests...")
    stats_results = perform_statistical_tests(df, metrics)
    
    # Save statistics
    with open(snakemake.output.stats, 'w') as f:
        json.dump(stats_results, f, indent=2)
    print(f"Saved statistics to {snakemake.output.stats}")
    
    # Create plots
    print("Creating plots...")
    create_comparison_plot(df, metrics, Path(snakemake.output.plot))
    print(f"Saved plot to {snakemake.output.plot}")
    
    # Print summary
    print("\n" + "="*70)
    print("REFINEMENT COMPARISON SUMMARY")
    print("="*70)
    
    for metric, result in stats_results.items():
        print(f"\n{metric}:")
        print(f"  Unrefined: {result['unrefined_mean']:.3f} ± {result['unrefined_std']:.3f}")
        print(f"  Refined:   {result['refined_mean']:.3f} ± {result['refined_std']:.3f}")
        print(f"  Change:    {result['mean_delta']:.3f} ± {result['std_delta']:.3f}")
        print(f"  Paired t-test: t={result['t_statistic']:.3f}, p={result['t_pvalue']:.4f} {'*' if result['significant_t'] else ''}")
        print(f"  Wilcoxon test: W={result['wilcoxon_statistic']:.1f}, p={result['wilcoxon_pvalue']:.4f} {'*' if result['significant_w'] else ''}")
        print(f"  Effect size (Cohen's d): {result['cohens_d']:.3f} ({result['effect_size_interpretation']})")
    
    print("="*70)
    print("* p < 0.05 (statistically significant)")
    print("="*70)


if __name__ == '__main__':
    main(snakemake)
