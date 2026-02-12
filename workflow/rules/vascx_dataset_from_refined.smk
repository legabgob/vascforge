# workflow/rules/vascx_dataset_from_refined.smk
import os
import re
import shutil
from pathlib import Path
from snakemake.io import glob_wildcards, directory

LEGACY_ROOT = config.get("legacy_root", ".")
VASCX_VIEW_ROOT = config.get("vascx", {}).get("dataset_view_root", "results/vascx_datasets")

# Where your grayscale label maps live (produced by refined_rgb_to_labels)
REF_LABEL_ROOT = config.get("labels_out_root", "results/refined_labels")

# --- Discover datasets (and other_dir) from seg_legacy/original ---
pat_simple = os.path.join(LEGACY_ROOT, "{dataset}", "seg_legacy", "original")
pat_other  = os.path.join(LEGACY_ROOT, "{dataset}", "{other_dir}", "seg_legacy", "original")

ds1 = []
if os.path.exists(os.path.join(LEGACY_ROOT)):
    ds1 = glob_wildcards(pat_simple)[0]

ds2, od2 = [], []
try:
    ds2, od2 = glob_wildcards(pat_other)
except ValueError:
    ds2, od2 = [], []

OTHERDIR_DATASETS = sorted(set(ds2))
SIMPLE_DATASETS = sorted(set(ds1) - set(ds2))

OTHERDIRS = {}
for d, od in zip(ds2, od2):
    OTHERDIRS.setdefault(d, set()).add(od)
OTHERDIRS = {d: sorted(v) for d, v in OTHERDIRS.items()}

RESOLUTIONS = [str(r) for r in config.get("resolutions", ["576", "1024"])]
k_start, k_end = config.get("k_range", [3, 9])
K_VALUES = list(range(int(k_start), int(k_end)))


def original_dir_simple(wc):
    return os.path.join(LEGACY_ROOT, wc.dataset, "seg_legacy", "original")

def original_dir_other(wc):
    return os.path.join(LEGACY_ROOT, wc.dataset, wc.other_dir, "seg_legacy", "original")


rule make_vascx_view_simple:
    wildcard_constraints:
        dataset="[^/]+",  # No slashes allowed
    input:
        # All required folders from seg_legacy
        original = original_dir_simple,
        av = "results/refined_labels/{dataset}/k{k}/downsampled/{res}px",
        meta = "data/{dataset}/meta/meta_filtered.csv",
        # Additional folders that VascX expects
        ce = lambda wc: os.path.join(LEGACY_ROOT, wc.dataset, "seg_legacy", "ce"),
        rgb = lambda wc: os.path.join(LEGACY_ROOT, wc.dataset, "seg_legacy", "rgb"),
        discs = lambda wc: os.path.join(LEGACY_ROOT, wc.dataset, "seg_legacy", "discs"),
        fovea = lambda wc: os.path.join(LEGACY_ROOT, wc.dataset, "seg_legacy", "fovea.csv"),
    output:
        view = directory(f"{VASCX_VIEW_ROOT}" + "/{dataset}/k{k}/downsampled/{res}px"),
    run:
        view = Path(output.view)
        view.mkdir(parents=True, exist_ok=True)

        # Create symlinks for all required components
        links = {
            "original": input.original,
            "av": input.av,
            "meta.csv": input.meta,
            "ce": input.ce,
            "rgb": input.rgb,
            "discs": input.discs,
            "fovea.csv": input.fovea,
        }
        
        for link_name, target in links.items():
            link_path = view / link_name
            
            # Remove existing link/file if present
            if link_path.exists() or link_path.is_symlink():
                link_path.unlink()
            
            # Create symlink with absolute path
            os.symlink(os.path.abspath(target), link_path)
        
        print(f"[make_vascx_view_simple] Created symlinks for {len(links)} components")


rule make_vascx_view_otherdir:
    wildcard_constraints:
        dataset="[^/]+",
        other_dir="[^/]+",
    input:
        # All required folders from seg_legacy
        original = original_dir_other,
        av = "results/refined_labels/{dataset}/k{k}/downsampled/{res}px",
        meta = "data/{dataset}/meta/meta_filtered.csv",
        # Additional folders that VascX expects
        ce = lambda wc: os.path.join(LEGACY_ROOT, wc.dataset, wc.other_dir, "seg_legacy", "ce"),
        rgb = lambda wc: os.path.join(LEGACY_ROOT, wc.dataset, wc.other_dir, "seg_legacy", "rgb"),
        discs = lambda wc: os.path.join(LEGACY_ROOT, wc.dataset, wc.other_dir, "seg_legacy", "discs"),
        fovea = lambda wc: os.path.join(LEGACY_ROOT, wc.dataset, wc.other_dir, "seg_legacy", "fovea.csv"),
    output:
        view = directory(f"{VASCX_VIEW_ROOT}" + "/{dataset}/{other_dir}/k{k}/downsampled/{res}px"),
    run:
        view = Path(output.view)
        view.mkdir(parents=True, exist_ok=True)

        # Create symlinks for all required components
        links = {
            "original": input.original,
            "av": input.av,
            "meta.csv": input.meta,
            "ce": input.ce,
            "rgb": input.rgb,
            "discs": input.discs,
            "fovea.csv": input.fovea,
        }
        
        for link_name, target in links.items():
            link_path = view / link_name
            
            # Remove existing link/file if present
            if link_path.exists() or link_path.is_symlink():
                link_path.unlink()
            
            # Create symlink with absolute path
            os.symlink(os.path.abspath(target), link_path)
        
        print(f"[make_vascx_view_otherdir] Created symlinks for {len(links)} components")
