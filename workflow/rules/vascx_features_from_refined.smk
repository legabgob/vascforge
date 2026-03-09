# workflow/rules/vascx_features_from_refined.smk
import re

VASCX_VIEW_ROOT = config.get("vascx", {}).get("dataset_view_root", "results/vascx_datasets")
FEATURES_OUT = config.get("vascx", {}).get("features_refined_out", "results/vascx_features_refined")
N_JOBS = int(config.get("vascx", {}).get("n_jobs", 64))

# Expects these to come from vascx_dataset_from_refined.smk:
# SIMPLE_DATASETS, OTHERDIR_DATASETS, OTHERDIRS, RESOLUTIONS, K_VALUES

rule vascx_features_refined_simple:
    wildcard_constraints:
        dataset="|".join(map(re.escape, SIMPLE_DATASETS)) if SIMPLE_DATASETS else "NO_MATCH"
    input:
        ds_dir = f"{VASCX_VIEW_ROOT}" + "/{dataset}/k{k}/downsampled/{res}px"
    output:
        features = f"{FEATURES_OUT}" + "/{dataset}/k{k}/downsampled/{res}px/vascx_features.tsv"
    benchmark:
        "benchmarks/vascx_features_refined_simple/{dataset}_k{k}_{res}px.tsv"
    params:
        name = "vascx_features_refined_simple",
        time = "08:00:00",
        mem = 16000,
        threads = 8,
    log:
        "logs/vascx_features/{dataset}_k{k}_{res}px.log"
    shell:
        r"""
        python ./workflow/scripts/run_full_pipeline.py {input.ds_dir} \
            --skip-preprocessing \
            --skip-segmentation \
            --n-jobs {N_JOBS} \
            2>&1 | tee {log}
        
        # Move output to expected location (run_full_pipeline.py creates dated filename)
        mkdir -p $(dirname {output.features})
        mv {input.ds_dir}/../extracted_features/*_vascx_features.csv {output.features}
        """

rule vascx_features_refined_otherdir:
    wildcard_constraints:
        dataset="|".join(map(re.escape, OTHERDIR_DATASETS)) if OTHERDIR_DATASETS else "NO_MATCH"
    input:
        ds_dir = f"{VASCX_VIEW_ROOT}" + "/{dataset}/{other_dir}/k{k}/downsampled/{res}px"
    output:
        features = f"{FEATURES_OUT}" + "/{dataset}/{other_dir}/k{k}/downsampled/{res}px/vascx_features.tsv"
    benchmark:
        "benchmarks/vascx_features_refined_otherdir/{dataset}_{other_dir}_k{k}_{res}px.tsv"
    params:
        name = "vascx_features_refined_otherdir",
        time = "08:00:00",
        mem = 16000,
        threads = 8,
    log:
        "logs/vascx_features/{dataset}_{other_dir}_k{k}_{res}px.log"
    shell:
        r"""
        python ./workflow/scripts/run_full_pipeline.py {input.ds_dir} \
            --skip-preprocessing \
            --skip-segmentation \
            --n-jobs {N_JOBS} \
            2>&1 | tee {log}
        
        # Move output to expected location (run_full_pipeline.py creates dated filename)
        mkdir -p $(dirname {output.features})
        mv {input.ds_dir}/../extracted_features/*_vascx_features.csv {output.features}
        """
