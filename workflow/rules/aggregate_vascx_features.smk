# workflow/rules/aggregate_vascx_features.smk
import pandas as pd
from pathlib import Path

VASCX_AGGREGATED_OUT = config.get("vascx", {}).get("aggregated_out", "results/vascx_aggregated")

def get_vascx_inputs_simple(wc):
    """Get all vascx feature files for a simple dataset at a given resolution"""
    inputs = []
    
    # Add unrefined (k=0) if enabled
    if config.get("vascx", {}).get("unrefined_features_enabled", False):
        unrefined_out = config.get("vascx", {}).get("features_unrefined_out", "results/vascx_features_unrefined")
        inputs.append(f"{unrefined_out}/{wc.dataset}/downsampled/{wc.res}px/vascx_features.tsv")
    
    # Add refined (k=3,4,...,8) if enabled
    if config.get("vascx", {}).get("refined_features_enabled", True):
        refined_out = config.get("vascx", {}).get("features_refined_out", "results/vascx_features_refined")
        for k in K_VALUES:
            inputs.append(f"{refined_out}/{wc.dataset}/k{k}/downsampled/{wc.res}px/vascx_features.tsv")
    
    return inputs

def get_vascx_inputs_otherdir(wc):
    """Get all vascx feature files for an otherdir dataset at a given resolution"""
    inputs = []
    
    # Add unrefined (k=0) if enabled
    if config.get("vascx", {}).get("unrefined_features_enabled", False):
        unrefined_out = config.get("vascx", {}).get("features_unrefined_out", "results/vascx_features_unrefined")
        inputs.append(f"{unrefined_out}/{wc.dataset}/{wc.other_dir}/downsampled/{wc.res}px/vascx_features.tsv")
    
    # Add refined (k=3,4,...,8) if enabled
    if config.get("vascx", {}).get("refined_features_enabled", True):
        refined_out = config.get("vascx", {}).get("features_refined_out", "results/vascx_features_refined")
        for k in K_VALUES:
            inputs.append(f"{refined_out}/{wc.dataset}/{wc.other_dir}/k{k}/downsampled/{wc.res}px/vascx_features.tsv")
    
    return inputs

rule aggregate_vascx_simple:
    """Aggregate vascx features across all k values (including unrefined k=0) for simple datasets"""
    wildcard_constraints:
        dataset="|".join(map(re.escape, SIMPLE_DATASETS)) if SIMPLE_DATASETS else "NO_MATCH"
    input:
        get_vascx_inputs_simple
    output:
        aggregated = f"{VASCX_AGGREGATED_OUT}" + "/{dataset}/downsampled/{res}px/vascx_features_aggregated.tsv"
    run:
        import pandas as pd
        from pathlib import Path
        
        all_dfs = []
        
        # Determine k values based on input files
        unrefined_enabled = config.get("vascx", {}).get("unrefined_features_enabled", False)
        refined_enabled = config.get("vascx", {}).get("refined_features_enabled", True)
        
        # Process unrefined (k=0) first if it exists
        if unrefined_enabled:
            unrefined_path = Path(input[0])
            if unrefined_path.exists():
                df_unrefined = pd.read_csv(unrefined_path, sep='\t')
                df_unrefined['k'] = 0
                all_dfs.append(df_unrefined)
        
        # Process refined (k=3,4,...,8)
        if refined_enabled:
            start_idx = 1 if unrefined_enabled else 0
            for i, k in enumerate(K_VALUES):
                refined_path = Path(input[start_idx + i])
                if refined_path.exists():
                    df_refined = pd.read_csv(refined_path, sep='\t')
                    df_refined['k'] = k
                    all_dfs.append(df_refined)
        
        if not all_dfs:
            raise ValueError(f"No vascx feature files found for {wildcards.dataset} at {wildcards.res}px")
        
        # Concatenate all dataframes
        df_combined = pd.concat(all_dfs, ignore_index=True)
        
        # Add non_nan_count column (count non-NaN values per row, excluding 'k' column)
        feature_cols = [col for col in df_combined.columns if col != 'k']
        df_combined['non_nan_count'] = df_combined[feature_cols].notna().sum(axis=1)
        
        # Reorder columns: put k and non_nan_count first
        cols = ['k', 'non_nan_count'] + [col for col in df_combined.columns if col not in ['k', 'non_nan_count']]
        df_combined = df_combined[cols]
        
        # Save aggregated features
        output_path = Path(output.aggregated)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df_combined.to_csv(output_path, sep='\t', index=False)
        
        print(f"[aggregate_vascx_simple] Aggregated {len(all_dfs)} feature sets for {wildcards.dataset} at {wildcards.res}px ({len(df_combined)} total rows)")

rule aggregate_vascx_otherdir:
    """Aggregate vascx features across all k values (including unrefined k=0) for otherdir datasets"""
    wildcard_constraints:
        dataset="|".join(map(re.escape, OTHERDIR_DATASETS)) if OTHERDIR_DATASETS else "NO_MATCH"
    input:
        get_vascx_inputs_otherdir
    output:
        aggregated = f"{VASCX_AGGREGATED_OUT}" + "/{dataset}/{other_dir}/downsampled/{res}px/vascx_features_aggregated.tsv"
    run:
        import pandas as pd
        from pathlib import Path
        
        all_dfs = []
        
        # Determine k values based on input files
        unrefined_enabled = config.get("vascx", {}).get("unrefined_features_enabled", False)
        refined_enabled = config.get("vascx", {}).get("refined_features_enabled", True)
        
        # Process unrefined (k=0) first if it exists
        if unrefined_enabled:
            unrefined_path = Path(input[0])
            if unrefined_path.exists():
                df_unrefined = pd.read_csv(unrefined_path, sep='\t')
                df_unrefined['k'] = 0
                all_dfs.append(df_unrefined)
        
        # Process refined (k=3,4,...,8)
        if refined_enabled:
            start_idx = 1 if unrefined_enabled else 0
            for i, k in enumerate(K_VALUES):
                refined_path = Path(input[start_idx + i])
                if refined_path.exists():
                    df_refined = pd.read_csv(refined_path, sep='\t')
                    df_refined['k'] = k
                    all_dfs.append(df_refined)
        
        if not all_dfs:
            raise ValueError(f"No vascx feature files found for {wildcards.dataset}/{wildcards.other_dir} at {wildcards.res}px")
        
        # Concatenate all dataframes
        df_combined = pd.concat(all_dfs, ignore_index=True)
        
        # Add non_nan_count column (count non-NaN values per row, excluding 'k' column)
        feature_cols = [col for col in df_combined.columns if col != 'k']
        df_combined['non_nan_count'] = df_combined[feature_cols].notna().sum(axis=1)
        
        # Reorder columns: put k and non_nan_count first
        cols = ['k', 'non_nan_count'] + [col for col in df_combined.columns if col not in ['k', 'non_nan_count']]
        df_combined = df_combined[cols]
        
        # Save aggregated features
        output_path = Path(output.aggregated)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df_combined.to_csv(output_path, sep='\t', index=False)
        
        print(f"[aggregate_vascx_otherdir] Aggregated {len(all_dfs)} feature sets for {wildcards.dataset}/{wildcards.other_dir} at {wildcards.res}px ({len(df_combined)} total rows)")
