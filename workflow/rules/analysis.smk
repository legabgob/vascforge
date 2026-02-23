# workflow/rules/analysis.smk
import re

# ── Slug helpers ──────────────────────────────────────────────────────────────

def _slug(name: str) -> str:
    """Replace non-alphanumeric characters with '_' for use as Snakemake input keys."""
    return re.sub(r"[^A-Za-z0-9]", "_", name)


def _collect_dice_inputs(wc):
    """
    Build the dice CSV input dict dynamically from the metrics datasets defined
    in dice.smk (METRICS_AV_SIMPLE, METRICS_AV_OTHERDIR, METRICS_AV_OTHERDIR_VALUES).

    Key convention (must match analysis_smk.R):
      - simple dataset:        Fundus_AVSeg
      - otherdir dataset/split: leuven_haifa__train
        (double underscore separates dataset from other_dir)
    """
    res = str(ANALYSIS_CFG.get("res", "1024"))
    inputs = {}

    for ds in METRICS_AV_SIMPLE:
        inputs[_slug(ds)] = f"{METRICS_DIR}/{ds}/metrics_{res}.csv"

    for ds in METRICS_AV_OTHERDIR:
        for od in METRICS_AV_OTHERDIR_VALUES.get(ds, []):
            key = f"{_slug(ds)}__{_slug(od)}"
            inputs[key] = f"{METRICS_DIR}/{ds}/{od}/metrics_{res}.csv"

    return inputs


# ── Config ────────────────────────────────────────────────────────────────────

ANALYSIS_CFG = config.get("analysis", {})

k_start, k_end = config.get("k_range", [3, 9])
_K_REFINED = list(range(int(k_start), int(k_end)))

# ── Rule ──────────────────────────────────────────────────────────────────────

rule analysis:
    """
    Combined analysis rule: connectivity plots, top-CC extraction, and DICE plots.
    DICE datasets are collected automatically from the same metrics datasets
    defined in dice.smk — no manual listing needed.
    """
    input:
        unpack(_collect_dice_inputs),
        connectivity = "results/connectivity/master_connectivity.csv",
    output:
        out_dir = directory("results/analysis/plots"),
    conda:
        "../envs/r_analysis.yaml"
    params:
        k_refined  = _K_REFINED,
        k_compare  = ANALYSIS_CFG.get("k_compare", 5),
        k_top      = ANALYSIS_CFG.get("k_top",     8),
        top_n      = ANALYSIS_CFG.get("top_n",     3),
    script:
        "../scripts/analysis_smk.R"
