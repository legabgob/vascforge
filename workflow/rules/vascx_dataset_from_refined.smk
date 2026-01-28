# workflow/rules/vascx_dataset_from_refined.smk
import os
import re
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
        dataset="|".join(map(re.escape, SIMPLE_DATASETS)) if SIMPLE_DATASETS else "NO_MATCH"
    input:
        original = original_dir_simple,
        # directory() is outputs-only; keep as plain path string for input
        av = f"{REF_LABEL_ROOT}" + "/{dataset}/k{k}/downsampled/{res}px",
        meta = "data/{dataset}/meta/meta_filtered.csv",
    output:
        view = directory(f"{VASCX_VIEW_ROOT}" + "/{dataset}/k{k}/downsampled/{res}px"),
    run:
        view = Path(output.view)
        view.mkdir(parents=True, exist_ok=True)

        # Create subdirs
        original_dst = view / "original"
        av_dst = view / "av"
        meta_dst = view / "meta.csv"
        original_dst.mkdir(parents=True, exist_ok=True)
        av_dst.mkdir(parents=True, exist_ok=True)

        # Hardlink files 
        original_src = Path(input.original)
        av_src = Path(input.av)
        meta_src = Path(input.meta)

        # Hardlink original files
        for src_file in original_src.glob("*.png"):
            dst_file = original_dst / src_file.name
            if dst_file.exists():
                dst_file.unlink()
            os.link(src_file, dst_file)
  
        # Hardlink av files
        for src_file in av_src.glob("*.png"):
            dst_file = av_dst / src_file.name
            if dst_file.exists():
                dst_file.unlink()
            os.link(src_file, dst_file)
        # Hardlink meta files
        if meta_dst.exists():
            meta_dst.unlink()
        os.link(meta_src, meta_dst)

        print(f"[make_vascx_view_simple] Hardlinked {len(list(original_dst.glob('*.png')))} original files.")
        print(f"[make_vascx_view_simple] Hardlinked {len(list(av_dst.glob('*.png')))} av files.")
        print(f"[make_vascx_view_simple] Hardlinked meta file.")

    #orig_link = view / "original"
    #   av_link = view / "av"
    #   meta_link = view / "meta.csv"
    #
    #   for link in (orig_link, av_link, meta_link):
    #            if link.exists() or link.is_symlink():
    #           link.unlink()
    #
    #   os.symlink(os.path.abspath(input.original), orig_link)
    #   os.symlink(os.path.abspath(input.av), av_link)
  #   os.symlink(os.path.abspath(input.meta), meta_link)


rule make_vascx_view_otherdir:
    wildcard_constraints:
        dataset="|".join(map(re.escape, OTHERDIR_DATASETS)) if OTHERDIR_DATASETS else "NO_MATCH"
    input:
        original = original_dir_other,
        av = f"{REF_LABEL_ROOT}" + "/{dataset}/k{k}/downsampled/{res}px",
        meta = "data/{dataset}/meta/meta_filtered.csv",
    output:
        view = directory(f"{VASCX_VIEW_ROOT}" + "/{dataset}/{other_dir}/k{k}/downsampled/{res}px"),
    run:
        view = Path(output.view)
        view.mkdir(parents=True, exist_ok=True)

        original_dst = view / "original"
        av_dst = view / "av"
        meta_dst = view / "meta.csv"

        original_dst.mkdir(parents=True, exist_ok=True)
        av_dst.mkdir(parents=True, exist_ok=True)

        # Hardlink files 
        original_src = Path(input.original)
        av_src = Path(input.av)
        meta_src = Path(input.meta)

        # Hardlink original files
        for src_file in original_src.glob("*.png"):
            dst_file = original_dst / src_file.name
            if dst_file.exists():
                dst_file.unlink()
            os.link(src_file, dst_file)

        # Hardlink av files
        for src_file in av_src.glob("*.png"):
            dst_file = av_dst / src_file.name
            if dst_file.exists():
                dst_file.unlink()
            os.link(src_file, dst_file)

        # Hardlink meta files
        if meta_dst.exists():
            meta_dst.unlink()
        os.link(meta_src, meta_dst)

        print(f"[make_vascx_view_otherdir] Hardlinked {len(list(original_dst.glob('*.png')))} original files.")
        print(f"[make_vascx_view_otherdir] Hardlinked {len(list(av_dst.glob('*.png')))} av files.")
        print(f"[make_vascx_view_otherdir] Hardlinked meta file.")


        

    #for link in (orig_link, av_link, meta_link):
    #       if link.exists() or link.is_symlink():
    #           link.unlink()
    #
    #   os.symlink(os.path.abspath(input.original), orig_link)
    #   os.symlink(os.path.abspath(input.av), av_link)
    #   os.symlink(os.path.abspath(input.meta), meta_link)

