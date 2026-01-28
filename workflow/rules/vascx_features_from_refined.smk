# workflow/rules/vascx_features_standalone.smk

VASCX_VIEW_ROOT = config.get("vascx", {}).get("dataset_view_root", "results/vascx_datasets")
FEATURES_OUT = config.get("vascx", {}).get("features_refined_out", "results/vascx_features_refined")
N_JOBS = int(config.get("vascx", {}).get("n_jobs", 64))
AV_SUBFOLDER = str(config.get("vascx", {}).get("av_subfolder", "av"))

rule vascx_features_refined_otherdir:
    input:
        ds_dir = f"{VASCX_VIEW_ROOT}" + "/{dataset}/{other_dir}/k{k}/downsampled/{res}px"
    output:
        features = f"{FEATURES_OUT}" + "/{dataset}/{other_dir}/k{k}/downsampled/{res}px/vascx_features.tsv"
    params:
        n_jobs = N_JOBS,
        av_subfolder = AV_SUBFOLDER
    log:
        "logs/vascx_features/{dataset}_{other_dir}_k{k}_{res}px.log"
    shell:
        """
        python workflow/scripts/vascx_extract_standalone.py \
            {input.ds_dir} \
            {output.features} \
            --n-jobs {params.n_jobs} \
            --av-subfolder {params.av_subfolder} \
            2>&1 | tee {log}
        """

