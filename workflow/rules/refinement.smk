# workflow/rules/refinement.smk
from snakemake.io import directory

rule refinement:
    input:
        weights = lambda wc: (
            "data/weights/rrwnet_HRF_0.pth" if wc.res == "1024"
            else "data/weights/rrwnet_RITE_refinement.pth"
        ),
        # Use filtered directories that contain only square images
        segmentations = "data/{dataset}/downsampled/{res}px/segs_converted_square",
        masks = "data/{dataset}/downsampled/{res}px/roi_masks_binarized_square",
    output:
        refined = directory("results/refined/{dataset}/k{k}/downsampled/{res}px")
    resources:
        gpu_jobs = 1  # Limit GPU concurrent jobs to prevent OOM
    shell:
        r"""
        {config[rrwnet][python]} {config[rrwnet][script]} \
            --weights {input.weights} \
            --images-path {input.segmentations} \
            --masks-path {input.masks} \
            --save-path {output.refined} \
            --k {wildcards.k} \
            --refine
        """
