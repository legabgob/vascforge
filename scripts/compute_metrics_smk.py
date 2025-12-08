from pathlib import Path
import re
import numpy as np
import pandas as pd
from PIL import Image

# ---------------------- Common utilities ---------------------- #

EXT = ".png"


def dice_bool(a: np.ndarray, b: np.ndarray, eps: float = 1e-8) -> float:
    inter = np.logical_and(a, b).sum(dtype=np.float64)
    return float((2.0 * inter + eps) / (a.sum(dtype=np.float64) + b.sum(dtype=np.float64) + eps))


def dice_bool_masked(a: np.ndarray, b: np.ndarray, roi: np.ndarray, eps: float = 1e-8) -> float:
    """
    Dice between boolean masks a and b, but only within roi==True.
    If ROI has zero True pixels, returns np.nan (no area to evaluate).
    """
    if roi.dtype != bool:
        roi = roi.astype(bool)
    n_roi = roi.sum()
    if n_roi == 0:
        return np.nan
    aR = np.logical_and(a, roi)
    bR = np.logical_and(b, roi)
    inter = np.logical_and(aR, bR).sum(dtype=np.float64)
    return float((2.0 * inter + eps) / (aR.sum(dtype=np.float64) + bR.sum(dtype=np.float64) + eps))


def find_k_dirs_for_res(refined_root: Path, res_subdir: Path) -> dict[int, Path]:
    """
    Finds k* subfolders under refined_root that contain the given res_subdir.
    Returns {K: Path(.../kK/res_subdir)}.
    """
    rex = re.compile(r"^[kK]\s*=?\s*(\d+)$")
    out: dict[int, Path] = {}
    for d in refined_root.iterdir():
        if d.is_dir():
            m = rex.match(d.name)
            if m:
                K = int(m.group(1))
                p = d / res_subdir
                if p.is_dir():
                    out[K] = p
    return dict(sorted(out.items()))

# ---------------------- A/V GT (Fundus-AVSeg style) ---------------------- #

def load_rgb_arr(p: Path) -> np.ndarray:
    return np.array(Image.open(p).convert("RGB"), dtype=np.uint8)


def load_roi_bool(p: Path) -> np.ndarray:
    """
    Load ROI as boolean mask: True = inside ROI (nonzero).
    """
    arr = np.array(Image.open(p).convert("L"))
    return arr > 0


def av_masks(arr: np.ndarray):
    """A from Red>0, V from Green>0. (Blue ignored for GT comparisons.)"""
    A = arr[..., 0] > 0
    V = arr[..., 1] > 0
    return A, V


def rgb_bool_channels(arr: np.ndarray):
    """R/G/B channel > 0 as booleans (for change metrics)."""
    return arr[..., 0] > 0, arr[..., 1] > 0, arr[..., 2] > 0


def av_macro_dice(A1, V1, A2, V2) -> float:
    """Macro-average Dice over A and V."""
    return 0.5 * (dice_bool(A1, A2) + dice_bool(V1, V2))


def av_macro_dice_masked(A1, V1, A2, V2, roi: np.ndarray) -> float:
    return 0.5 * (dice_bool_masked(A1, A2, roi) + dice_bool_masked(V1, V2, roi))


def build_df_av_for_res(
    gt_dir: Path,
    unref_dir: Path,
    k_dirs: dict[int, Path],
    roi_dir: Path | None = None,
) -> pd.DataFrame:
    """
    A/V macro dice metrics (Fundus-AVSeg style), for a given resolution.
    """
    if not k_dirs:
        raise SystemExit("No k* subfolders provided")

    unref_names = {p.name for p in unref_dir.glob(f"*{EXT}") if p.is_file()}
    gt_names = {p.name for p in gt_dir.glob(f"*{EXT}") if p.is_file()}
    names = sorted(unref_names & gt_names)
    if not names:
        raise SystemExit(f"No overlapping {EXT} files between GT {gt_dir} and UNREF {unref_dir}")

    rows = []
    for name in names:
        p_unref = unref_dir / name
        p_gt = gt_dir / name

        arr_unref = load_rgb_arr(p_unref)
        arr_gt = load_rgb_arr(p_gt)

        if arr_gt.shape != arr_unref.shape:
            raise SystemExit(f"Size mismatch for {name}: GT {arr_gt.shape} vs UNREF {arr_unref.shape}")

        # Load ROI if provided; else full-True mask
        if roi_dir is not None:
            p_roi = roi_dir / name
            if not p_roi.exists():
                raise SystemExit(f"ROI mask missing: {p_roi}")
            roi = load_roi_bool(p_roi)
            if roi.shape != arr_unref.shape[:2]:
                raise SystemExit(
                    f"ROI size mismatch for {name}: ROI {roi.shape} vs UNREF {arr_unref.shape[:2]}"
                )
        else:
            roi = np.ones(arr_unref.shape[:2], dtype=bool)

        Au, Vu = av_masks(arr_unref)
        Ag, Vg = av_masks(arr_gt)

        row: dict[str, float | str] = {"image": name}
        # Unref vs GT (macro + per-class) within ROI
        row["dice_unref_vs_gt_macro"] = av_macro_dice_masked(Au, Vu, Ag, Vg, roi)
        row["dice_unref_vs_gt_A"] = dice_bool_masked(Au, Ag, roi)
        row["dice_unref_vs_gt_V"] = dice_bool_masked(Vu, Vg, roi)

        # Precompute UNREF per-channel bools (for change metrics)
        r_u, g_u, b_u = rgb_bool_channels(arr_unref)

        for K, dK in k_dirs.items():
            p_ref = dK / name
            if not p_ref.exists():
                row[f"dice_k{K}_vs_gt_macro"] = np.nan
                row[f"dice_k{K}_vs_gt_A"] = np.nan
                row[f"dice_k{K}_vs_gt_V"] = np.nan
                row[f"dice_unref_vs_k{K}_rgbmacro"] = np.nan
                row[f"dice_unref_vs_k{K}_avmacro"] = np.nan
                row[f"dice_unref_vs_k{K}_A"] = np.nan
                row[f"dice_unref_vs_k{K}_V"] = np.nan
                continue

            arr_ref = load_rgb_arr(p_ref)
            if arr_ref.shape != arr_unref.shape:
                raise SystemExit(
                    f"Size mismatch for {name} at k={K}: REF {arr_ref.shape} vs UNREF {arr_unref.shape}"
                )

            Ar, Vr = av_masks(arr_ref)

            # Refined vs GT (macro + per-class), masked
            row[f"dice_k{K}_vs_gt_macro"] = av_macro_dice_masked(Ar, Vr, Ag, Vg, roi)
            row[f"dice_k{K}_vs_gt_A"] = dice_bool_masked(Ar, Ag, roi)
            row[f"dice_k{K}_vs_gt_V"] = dice_bool_masked(Vr, Vg, roi)

            # Change metrics within ROI
            r_r, g_r, b_r = rgb_bool_channels(arr_ref)
            row[f"dice_unref_vs_k{K}_rgbmacro"] = float(
                (
                    dice_bool_masked(r_u, r_r, roi)
                    + dice_bool_masked(g_u, g_r, roi)
                    + dice_bool_masked(b_u, b_r, roi)
                )
                / 3.0
            )
            row[f"dice_unref_vs_k{K}_avmacro"] = 0.5 * (
                dice_bool_masked(r_u, r_r, roi) + dice_bool_masked(g_u, g_r, roi)
            )
            row[f"dice_unref_vs_k{K}_A"] = dice_bool_masked(r_u, r_r, roi)
            row[f"dice_unref_vs_k{K}_V"] = dice_bool_masked(g_u, g_r, roi)

        rows.append(row)

    df = pd.DataFrame(rows).set_index("image")
    k_list = sorted(k_dirs.keys())
    cols: list[str] = ["dice_unref_vs_gt_macro", "dice_unref_vs_gt_A", "dice_unref_vs_gt_V"]
    for K in k_list:
        cols += [f"dice_k{K}_vs_gt_macro", f"dice_k{K}_vs_gt_A", f"dice_k{K}_vs_gt_V"]
    for K in k_list:
        cols += [
            f"dice_unref_vs_k{K}_rgbmacro",
            f"dice_unref_vs_k{K}_avmacro",
            f"dice_unref_vs_k{K}_A",
            f"dice_unref_vs_k{K}_V",
        ]
    cols = [c for c in cols if c in df.columns]
    return df[cols]

# ---------------------- Vessel-only GT (FIVES style) ---------------------- #

def load_gray_bool(p: Path) -> np.ndarray:
    """Load grayscale (L) as boolean mask: True = foreground (nonzero)."""
    a = np.array(Image.open(p).convert("L"))
    return a > 0


def load_rgb_bools(p: Path):
    """Load RGB and return three boolean masks for (R>0), (G>0), (B>0)."""
    a = np.array(Image.open(p).convert("RGB"))
    return a[..., 0] > 0, a[..., 1] > 0, a[..., 2] > 0


def rgb_macro_dice(p_unref: Path, p_ref: Path) -> float:
    """Macro-average DICE over R/G/B channels between two RGB label images."""
    r1, g1, b1 = load_rgb_bools(p_unref)
    r2, g2, b2 = load_rgb_bools(p_ref)
    return float(
        (dice_bool(r1, r2) + dice_bool(g1, g2) + dice_bool(b1, b2)) / 3.0
    )


def compute_seg_vs_gt_map(gt_dir: Path, seg_dir: Path) -> dict[str, float]:
    """
    Compute DICE between GT (grayscale) and segmentation (grayscale) at native res.
    Returns {filename: dice}.
    """
    gt_names = {p.name for p in gt_dir.glob(f"*{EXT}") if p.is_file()}
    seg_names = {p.name for p in seg_dir.glob(f"*{EXT}") if p.is_file()}
    common = sorted(gt_names & seg_names)
    if not common:
        raise SystemExit(f"No common {EXT} between {gt_dir} and {seg_dir}")
    out: dict[str, float] = {}
    for name in common:
        gt = load_gray_bool(gt_dir / name)
        seg = load_gray_bool(seg_dir / name)
        out[name] = dice_bool(gt, seg)
    return out


def build_df_vessel_for_res(
    seg_vs_gt_map: dict[str, float],
    unref_dir: Path,
    k_dirs: dict[int, Path],
) -> pd.DataFrame:
    """
    Build one DataFrame for given resolution (vessel-only):
      - 'dice_seg_vs_gt' comes from seg_vs_gt_map (native resolution).
      - One column per K: dice between unrefined RGB and refined RGB (macro over R/G/B).
    """
    if not k_dirs:
        raise SystemExit("No k* subfolders provided to build_df_for_res")

    unref_names = {p.name for p in unref_dir.glob(f"*{EXT}") if p.is_file()}
    if not unref_names:
        raise SystemExit(f"No {EXT} files in {unref_dir}")

    # union of refined names across all K
    refined_union = set()
    for d in k_dirs.values():
        refined_union |= {p.name for p in d.glob(f"*{EXT}") if p.is_file()}

    names = sorted(unref_names & refined_union)
    if not names:
        raise SystemExit(
            f"No overlapping {EXT} files between {unref_dir} and any provided k* directory"
        )

    rows = []
    for name in names:
        row: dict[str, float | str] = {
            "image": name,
            "dice_seg_vs_gt": seg_vs_gt_map.get(name, np.nan),
        }
        p_unref = unref_dir / name
        for K, dK in k_dirs.items():
            p_ref = dK / name
            row[f"dice_unref_vs_k{K}"] = (
                rgb_macro_dice(p_unref, p_ref) if p_ref.exists() else np.nan
            )
        rows.append(row)

    df = pd.DataFrame(rows).set_index("image")
    k_cols = [f"dice_unref_vs_k{K}" for K in sorted(k_dirs.keys())]
    return df[["dice_seg_vs_gt"] + k_cols]

# ---------------------- Snakemake entrypoint ---------------------- #

# Expected from Snakemake:
#   wildcards.dataset              -> dataset name (e.g. "Fundus-AVSeg", "FIVES")
#   params.has_av_gt               -> bool/int: 1 if A/V GT available, 0 if only vessel GT
#   params.refined_root            -> base refined root path
#   params.unref_1024_dir          -> unrefined segs at 1024px
#   params.unref_576_dir           -> unrefined segs at 576px
#
#   If has_av_gt:
#       params.gt_1024_dir, params.gt_576_dir, params.roi_1024_dir, params.roi_576_dir
#   Else:
#       params.gt_native_dir, params.seg_native_dir
#
#   output.csv_1024, output.csv_576 -> where to write CSVs

wc = snakemake.wildcards
params = snakemake.params

dataset = str(wc.dataset)
has_av_gt = bool(int(params.has_av_gt))

refined_root = Path(str(params.refined_root))

# Resolution subdirs under refined_root/kX/
res_subdir_1024 = Path("downsampled/1024px")
res_subdir_576 = Path("downsampled/576px")

# Find k* dirs for both resolutions
k_dirs_1024 = find_k_dirs_for_res(refined_root, res_subdir_1024)
if not k_dirs_1024:
    raise SystemExit(f"No k* with {res_subdir_1024} under {refined_root}")

k_dirs_576 = find_k_dirs_for_res(refined_root, res_subdir_576)
if not k_dirs_576:
    raise SystemExit(f"No k* with {res_subdir_576} under {refined_root}")

out_csv_1024 = Path(str(snakemake.output.csv_1024))
out_csv_576 = Path(str(snakemake.output.csv_576))

unref_1024_dir = Path(str(params.unref_1024_dir))
unref_576_dir = Path(str(params.unref_576_dir))

if has_av_gt:
    # A/V GT branch
    gt_1024_dir = Path(str(params.gt_1024_dir))
    gt_576_dir = Path(str(params.gt_576_dir))
    roi_1024_dir = Path(str(params.roi_1024_dir)) if params.roi_1024_dir else None
    roi_576_dir = Path(str(params.roi_576_dir)) if params.roi_576_dir else None

    df_1024 = build_df_av_for_res(gt_1024_dir, unref_1024_dir, k_dirs_1024, roi_dir=roi_1024_dir)
    df_576 = build_df_av_for_res(gt_576_dir, unref_576_dir, k_dirs_576, roi_dir=roi_576_dir)

    df_1024.to_csv(out_csv_1024, float_format="%.6f")
    df_576.to_csv(out_csv_576, float_format="%.6f")

    print(f"[{dataset}] A/V GT mode: saved {out_csv_1024} ({df_1024.shape}) and {out_csv_576} ({df_576.shape})")

else:
    # Vessel-only GT branch
    gt_native_dir = Path(str(params.gt_native_dir))
    seg_native_dir = Path(str(params.seg_native_dir))

    seg_vs_gt = compute_seg_vs_gt_map(gt_native_dir, seg_native_dir)

    df_1024 = build_df_vessel_for_res(seg_vs_gt, unref_1024_dir, k_dirs_1024)
    df_576 = build_df_vessel_for_res(seg_vs_gt, unref_576_dir, k_dirs_576)

    df_1024.to_csv(out_csv_1024, float_format="%.6f")
    df_576.to_csv(out_csv_576, float_format="%.6f")

    print(f"[{dataset}] vessel-only GT mode: saved {out_csv_1024} ({df_1024.shape}) and {out_csv_576} ({df_576.shape})")

