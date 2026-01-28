# workflow/rules/vascx_features_from_refined.smk
import re

VASCX_VIEW_ROOT = config.get("vascx", {}).get("dataset_view_root", "results/vascx_datasets")
FEATURES_OUT = config.get("vascx", {}).get("features_refined_out", "results/vascx_features_refined")

VASCX_CFG = config.get("vascx", {})
N_JOBS = int(VASCX_CFG.get("n_jobs", 64))
FEATURE_SET = str(VASCX_CFG.get("feature_set", "bergmann"))
AV_SUBFOLDER = str(VASCX_CFG.get("av_subfolder", "av"))

# Expects these to come from vascx_dataset_from_refined.smk:
# SIMPLE_DATASETS, OTHERDIR_DATASETS, OTHERDIRS, RESOLUTIONS, K_VALUES

rule vascx_features_refined_simple:
    wildcard_constraints:
        dataset="|".join(map(re.escape, SIMPLE_DATASETS)) if SIMPLE_DATASETS else "NO_MATCH"
    input:
        # directory() is outputs-only; keep as plain string for input
        ds_dir = f"{VASCX_VIEW_ROOT}" + "/{dataset}/k{k}/downsampled/{res}px"
    output:
        features = f"{FEATURES_OUT}" + "/{dataset}/k{k}/downsampled/{res}px/vascx_features.tsv"
    params:
        n_jobs = N_JOBS,
        feature_set = FEATURE_SET,
        av_subfolder = AV_SUBFOLDER,
        sep = "\t",
        na_rep = "NaN",
    log:
        err = "../logs/vascx_features_refined_simple_{dataset}_k{k}_{res}px.err",
        out = "../logs/vascx_features_refined_simple_{dataset}_k{k}_{res}px.out"
    script:
        "../scripts/vascx_feature_extract_smk.py"

rule vascx_features_refined_otherdir:
    wildcard_constraints:
        dataset="|".join(map(re.escape, OTHERDIR_DATASETS)) if OTHERDIR_DATASETS else "NO_MATCH"
    input:
        ds_dir = f"{VASCX_VIEW_ROOT}" + "/{dataset}/{other_dir}/k{k}/downsampled/{res}px"
    output:
        features = f"{FEATURES_OUT}" + "/{dataset}/{other_dir}/k{k}/downsampled/{res}px/vascx_features.tsv"
    params:
        n_jobs = N_JOBS,
        feature_set = FEATURE_SET,
        av_subfolder = AV_SUBFOLDER,
        sep = "\t",
        na_rep = "NaN",
    log:
        err = "../logs/vascx_features_refined_otherdir_{dataset}_{other_dir}_k{k}_{res}px.err",
        out = "../logs/vascx_features_refined_otherdir_{dataset}_{other_dir}_k{k}_{res}px.out"
    script:
        "../scripts/vascx_feature_extract_smk.py"
