# workflow/rules/vascx_features_from_unrefined.smk
import os
import re
from pathlib import Path
from snakemake.io import glob_wildcards, directory

LEGACY_ROOT = config.get("legacy_root", ".")
VASCX_UNREFINED_ROOT = config.get("vascx", {}).get("dataset_unrefined_root", "results/vascx_datasets_unrefined")
FEATURES_UNREFINED_OUT = config.get("vascx", {}).get("features_unrefined_out", "results/vascx_features_unrefined")
N_JOBS = int(config.get("vascx", {}).get("n_jobs", 64))

# Reuse dataset discovery variables from vascx_dataset_from_refined.smk
# SIMPLE_DATASETS, OTHERDIR_DATASETS should already be defined

RESOLUTIONS = [str(r) for r in config.get("resolutions", ["576", "1024"])]

def seg_legacy_simple(wc):
    return os.path.join(LEGACY_ROOT, wc.dataset, "seg_legacy")

def seg_legacy_other(wc):
    return os.path.join(LEGACY_ROOT, wc.dataset, wc.other_dir, "seg_legacy")


rule make_vascx_view_unrefined_simple:
    wildcard_constraints:
        dataset="[^/]+",
    input:
        seg_legacy = seg_legacy_simple,
        meta = "data/{dataset}/meta/meta_filtered_square.csv",
    output:
        view = directory(f"{VASCX_UNREFINED_ROOT}" + "/{dataset}/downsampled/{res}px"),
    run:
        from PIL import Image
        import shutil
        
        view = Path(output.view)
        view.mkdir(parents=True, exist_ok=True)
        seg_legacy = Path(input.seg_legacy)
        
        # Symlink folders that don't need processing
        for folder in ["original", "ce", "rgb", "discs"]:
            src = seg_legacy / folder
            if not src.exists():
                print(f"Warning: {src} does not exist, skipping")
                continue
            link_path = view / folder
            if link_path.exists() or link_path.is_symlink():
                link_path.unlink()
            os.symlink(os.path.abspath(src), link_path)
        
        # Symlink fovea.csv
        fovea_src = seg_legacy / "fovea.csv"
        if fovea_src.exists():
            fovea_link = view / "fovea.csv"
            if fovea_link.exists() or fovea_link.is_symlink():
                fovea_link.unlink()
            os.symlink(os.path.abspath(fovea_src), fovea_link)
        
        # Symlink meta.csv
        meta_link = view / "meta.csv"
        if meta_link.exists() or meta_link.is_symlink():
            meta_link.unlink()
        os.symlink(os.path.abspath(input.meta), meta_link)
        
        # Downsample av/ folder to target resolution
        av_src = seg_legacy / "av"
        av_dst = view / "av"
        av_dst.mkdir(exist_ok=True)
        
        target_res = int(wildcards.res)
        count = 0
        
        for img_file in av_src.glob("*.png"):
            img = Image.open(img_file)
            
            # Downsample if needed
            if img.size[0] != target_res or img.size[1] != target_res:
                img = img.resize((target_res, target_res), Image.LANCZOS)
            
            img.save(av_dst / img_file.name)
            count += 1
        
        print(f"[make_vascx_view_unrefined_simple] Processed {count} av images at {target_res}px")


rule make_vascx_view_unrefined_otherdir:
    wildcard_constraints:
        dataset="[^/]+",
        other_dir="[^/]+",
    input:
        seg_legacy = seg_legacy_other,
        meta = "data/{dataset}/meta/meta_filtered_square.csv",
    output:
        view = directory(f"{VASCX_UNREFINED_ROOT}" + "/{dataset}/{other_dir}/downsampled/{res}px"),
    run:
        from PIL import Image
        import shutil
        
        view = Path(output.view)
        view.mkdir(parents=True, exist_ok=True)
        seg_legacy = Path(input.seg_legacy)
        
        # Symlink folders that don't need processing
        for folder in ["original", "ce", "rgb", "discs"]:
            src = seg_legacy / folder
            if not src.exists():
                print(f"Warning: {src} does not exist, skipping")
                continue
            link_path = view / folder
            if link_path.exists() or link_path.is_symlink():
                link_path.unlink()
            os.symlink(os.path.abspath(src), link_path)
        
        # Symlink fovea.csv
        fovea_src = seg_legacy / "fovea.csv"
        if fovea_src.exists():
            fovea_link = view / "fovea.csv"
            if fovea_link.exists() or fovea_link.is_symlink():
                fovea_link.unlink()
            os.symlink(os.path.abspath(fovea_src), fovea_link)
        
        # Symlink meta.csv
        meta_link = view / "meta.csv"
        if meta_link.exists() or meta_link.is_symlink():
            meta_link.unlink()
        os.symlink(os.path.abspath(input.meta), meta_link)
        
        # Downsample av/ folder to target resolution
        av_src = seg_legacy / "av"
        av_dst = view / "av"
        av_dst.mkdir(exist_ok=True)
        
        target_res = int(wildcards.res)
        count = 0
        
        for img_file in av_src.glob("*.png"):
            img = Image.open(img_file)
            
            # Downsample if needed
            if img.size[0] != target_res or img.size[1] != target_res:
                img = img.resize((target_res, target_res), Image.LANCZOS)
            
            img.save(av_dst / img_file.name)
            count += 1
        
        print(f"[make_vascx_view_unrefined_otherdir] Processed {count} av images at {target_res}px")


rule vascx_features_unrefined_simple:
    wildcard_constraints:
        dataset="|".join(map(re.escape, SIMPLE_DATASETS)) if SIMPLE_DATASETS else "NO_MATCH"
    input:
        ds_dir = f"{VASCX_UNREFINED_ROOT}" + "/{dataset}/downsampled/{res}px"
    output:
        features = f"{FEATURES_UNREFINED_OUT}" + "/{dataset}/downsampled/{res}px/vascx_features.tsv"
    log:
        "logs/vascx_features_unrefined/{dataset}_{res}px.log"
    shell:
        r"""
        python ./workflow/scripts/run_full_pipeline.py {input.ds_dir} \
            --skip-preprocessing \
            --skip-segmentation \
            --n-jobs {N_JOBS} \
            2>&1 | tee {log}
        
        mkdir -p $(dirname {output.features})
        mv {input.ds_dir}/../extracted_features/*_vascx_features.csv {output.features}
        """


rule vascx_features_unrefined_otherdir:
    wildcard_constraints:
        dataset="|".join(map(re.escape, OTHERDIR_DATASETS)) if OTHERDIR_DATASETS else "NO_MATCH"
    input:
        ds_dir = f"{VASCX_UNREFINED_ROOT}" + "/{dataset}/{other_dir}/downsampled/{res}px"
    output:
        features = f"{FEATURES_UNREFINED_OUT}" + "/{dataset}/{other_dir}/downsampled/{res}px/vascx_features.tsv"
    log:
        "logs/vascx_features_unrefined/{dataset}_{other_dir}_{res}px.log"
    shell:
        r"""
        python ./workflow/scripts/run_full_pipeline.py {input.ds_dir} \
            --skip-preprocessing \
            --skip-segmentation \
            --n-jobs {N_JOBS} \
            2>&1 | tee {log}
        
        mkdir -p $(dirname {output.features})
        mv {input.ds_dir}/../extracted_features/*_vascx_features.csv {output.features}
        """
