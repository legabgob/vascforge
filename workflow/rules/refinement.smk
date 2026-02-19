# workflow/rules/refinement.smk
from snakemake.io import directory

rule refinement:
    input:
        weights = lambda wc: config["weights"][str(wc.res)],
        # These are produced by the downsample rule (directory outputs)
        segmentations = "data/{dataset}/downsampled/{res}px/segs_converted",
        masks = "data/{dataset}/downsampled/{res}px/roi_masks_binarized",
    output:
        refined = directory("results/refined/{dataset}/k{k}/downsampled/{res}px")
    resources:
        gpu_jobs=1
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

