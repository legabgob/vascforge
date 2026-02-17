#!/usr/bin/env python3
"""
Compare unrefined vs refined connectivity metrics.
Fixed to handle boolean JSON serialization.
"""

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path


def convert_for_json(obj):
    """Convert numpy/pandas types to JSON-serializable types."""
    if isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.bool_, bool)):  # ← Fix for boolean serialization
        return bool(obj)
    elif isinstance(obj, dict):
        return {k: convert_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_for_json(item) for item in obj]
    else:
        return obj


def calculate_effect_size(before, after):
    """Calculate Cohen's d effect size."""
    pooled_std = np.sqrt((before.std()**2 + after.std()**2) / 2)
    if pooled_std == 0:
        return 0
    return (after.mean() - before.mean()) / pooled_std


def perform_statistical_tests(before, after, metric_name):
    """Perform statistical tests comparing before and after."""
    results = {
        'metric': metric_name,
        'n_samples': len(before),
        'before_mean': float(before.mean()),
        'before_std': float(before.std()),
        'after_mean': float(after.mean()),
        'after_std': float(after.std()),
        'mean_change': float(after.mean() - before.mean()),
        'mean_change_pct': float((after.mean() - before.mean()) / before.mean() * 100) if before.mean() != 0 else 0,
    }
    
    # Paired t-test
    t_stat, t_pvalue = stats.ttest_rel(before, after)
    results['paired_ttest'] = {
        't_statistic': float(t_stat),
        'p_value': float(t_pvalue),
        'significant': bool(t_pvalue < 0.05),  # ← Convert to bool explicitly
    }
    
    # Wilcoxon signed-rank test (non-parametric alternative)
    try:
        w_stat, w_pvalue = stats.wilcoxon(before, after)
        results['wilcoxon'] = {
            'statistic': float(w_stat),
            'p_value': float(w_pvalue),
            'significant': bool(w_pvalue < 0.05),  # ← Convert to bool explicitly
        }
    except Exception as e:
        results['wilcoxon'] = {
            'error': str(e)
        }
    
    # Effect size (Cohen's d)
    effect_size = calculate_effect_size(before, after)
    results['effect_size'] = {
        'cohens_d': float(effect_size),
        'interpretation': (
            'negligible' if abs(effect_size) < 0.2 else
            'small' if abs(effect_size) < 0.5 else
            'medium' if abs(effect_size) < 0.8 else
            'large'
        )
    }
    
    return results


def create_comparison_plot(df, metrics, output_path):
    """Create before/after comparison plots."""
    n_metrics = len(metrics)
    fig, axes = plt.subplots(1, n_metrics, figsize=(6*n_metrics, 5))
    
    if n_metrics == 1:
        axes = [axes]
    
    for ax, metric in zip(axes, metrics):
        before = df[f'{metric}_unrefined']
        after = df[f'{metric}_refined']
        
        # Paired plot
        for i in range(len(df)):
            ax.plot([0, 1], [before.iloc[i], after.iloc[i]], 
                   'o-', alpha=0.2, color='gray')
        
        # Mean with error bars
        ax.errorbar([0, 1], 
                   [before.mean(), after.mean()],
                   yerr=[before.std(), after.std()],
                   fmt='o-', linewidth=3, markersize=10,
                   color='red', label='Mean ± SD')
        
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Unrefined', 'Refined'])
        ax.set_ylabel(metric.replace('_', ' ').title())
        ax.set_title(f'{metric.replace("_", " ").title()}\nChange: {after.mean() - before.mean():.2f}')
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def main(snakemake):
    """Main function called by Snakemake."""
    # Read inputs
    print("Loading data...")
    df_unrefined = pd.read_csv(snakemake.input.unrefined)
    df_refined = pd.read_csv(snakemake.input.refined)
    
    # Match images by ID
    merged = pd.merge(
        df_unrefined,
        df_refined,
        on='image_id',
        suffixes=('_unrefined', '_refined')
    )
    
    print(f"Loaded {len(merged)} matched image pairs")
    
    # Calculate changes
    print("Calculating changes...")
    metrics_to_compare = snakemake.params.metrics
    
    comparison_data = []
    for _, row in merged.iterrows():
        entry = {'image_id': row['image_id']}
        for metric in metrics_to_compare:
            unrefined_val = row[f'{metric}_unrefined']
            refined_val = row[f'{metric}_refined']
            entry[f'{metric}_unrefined'] = unrefined_val
            entry[f'{metric}_refined'] = refined_val
            entry[f'{metric}_change'] = refined_val - unrefined_val
            if unrefined_val != 0:
                entry[f'{metric}_change_pct'] = (refined_val - unrefined_val) / unrefined_val * 100
            else:
                entry[f'{metric}_change_pct'] = 0
        comparison_data.append(entry)
    
    comparison_df = pd.DataFrame(comparison_data)
    
    # Save comparison CSV
    comparison_df.to_csv(snakemake.output.comparison, index=False)
    print(f"Saved comparison to {snakemake.output.comparison}")
    
    # Perform statistical tests
    print("Performing statistical tests...")
    stats_results = {}
    
    for metric in metrics_to_compare:
        before = merged[f'{metric}_unrefined']
        after = merged[f'{metric}_refined']
        
        test_results = perform_statistical_tests(before, after, metric)
        stats_results[metric] = test_results
    
    # Convert all values for JSON serialization
    stats_results = convert_for_json(stats_results)
    
    # Save statistics JSON
    with open(snakemake.output.stats, 'w') as f:
        json.dump(stats_results, f, indent=2)
    print(f"Saved statistics to {snakemake.output.stats}")
    
    # Create comparison plot
    print("Creating plots...")
    create_comparison_plot(comparison_df, metrics_to_compare, snakemake.output.plot)
    print(f"Saved plot to {snakemake.output.plot}")
    
    # Print summary
    print("\n" + "="*70)
    print("COMPARISON SUMMARY")
    print("="*70)
    for metric in metrics_to_compare:
        result = stats_results[metric]
        print(f"\n{metric}:")
        print(f"  Before: {result['before_mean']:.2f} ± {result['before_std']:.2f}")
        print(f"  After:  {result['after_mean']:.2f} ± {result['after_std']:.2f}")
        print(f"  Change: {result['mean_change']:.2f} ({result['mean_change_pct']:.1f}%)")
        print(f"  Paired t-test: p={result['paired_ttest']['p_value']:.4f} {'*' if result['paired_ttest']['significant'] else 'ns'}")
        if 'wilcoxon' in result and 'p_value' in result['wilcoxon']:
            print(f"  Wilcoxon test: p={result['wilcoxon']['p_value']:.4f} {'*' if result['wilcoxon']['significant'] else 'ns'}")
        print(f"  Effect size (Cohen's d): {result['effect_size']['cohens_d']:.3f} ({result['effect_size']['interpretation']})")
    print("="*70)


if __name__ == '__main__':
    # For standalone testing
    import sys
    if len(sys.argv) > 1:
        class FakeSnakemake:
            pass
        snakemake = FakeSnakemake()
        # You can set up test inputs here
        main(snakemake)
