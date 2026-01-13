# scripts/refined_rgb_to_labels_dir_smk.py
from pathlib import Path
import numpy as np
from PIL import Image


def box_mean_3x3(arr: np.ndarray) -> np.ndarray:
    """
    Fast 3x3 mean using reflect padding.
    Returns float32 array of same shape.
    """
    arr = arr.astype(np.float32, copy=False)
    p = np.pad(arr, 1, mode="reflect")

    s = (
        p[0:-2, 0:-2] + p[0:-2, 1:-1] + p[0:-2, 2:] +
        p[1:-1, 0:-2] + p[1:-1, 1:-1] + p[1:-1, 2:] +
        p[2:,   0:-2] + p[2:,   1:-1] + p[2:,   2:]
    )
    return s / 9.0


def rgb_to_vascx_labels(
    img_rgb: np.ndarray,
    tb: int = 100,
    th: int = 200,
    delta: int = 15,
    tie_vote: bool = True,
) -> np.ndarray:
    """
    Convert RGB refined segmentation (uint8) -> labels {0,1,2,3}.

    Priority (first match wins):
      1) Crossing (3): B>=th & R>=th & G>=th
      2) Artery   (1): B>=tb & R>=G+delta
      3) Vein     (2): B>=tb & G>=R+delta
      4) Ties/ambiguous with B>=tb and |R-G|<delta:
           - if tie_vote: 3x3 local mean vote between R and G
           - else: deterministic tiebreak
      5) Background (0): everything else
    """
    if img_rgb.ndim != 3 or img_rgb.shape[2] != 3:
        raise ValueError("Expected an RGB image of shape (H, W, 3).")

    R = img_rgb[..., 0].astype(np.uint8)
    G = img_rgb[..., 1].astype(np.uint8)
    B = img_rgb[..., 2].astype(np.uint8)

    H, W = R.shape
    labels = np.zeros((H, W), dtype=np.uint8)

    cross = (B >= th) & (R >= th) & (G >= th)
    labels[cross] = 3

    vessel = (B >= tb) & ~cross

    art = vessel & (R.astype(np.int16) >= (G.astype(np.int16) + delta))
    vein = vessel & (G.astype(np.int16) >= (R.astype(np.int16) + delta))

    labels[art] = 1
    labels[vein] = 2

    ties = vessel & ~(art | vein)

    if tie_vote:
        meanR = box_mean_3x3(R)
        meanG = box_mean_3x3(G)

        promote_cross = ties & (meanR >= th) & (meanG >= th) & (B >= th)
        labels[promote_cross] = 3

        still_ties = ties & ~promote_cross
        labels[still_ties & (meanR > meanG)] = 1
        labels[still_ties & (meanG >= meanR)] = 2
    else:
        labels[ties & (R > G)] = 1
        labels[ties & (G >= R)] = 2

    return labels


def convert_one(src: Path, dst: Path, tb: int, th: int, delta: int, tie_vote: bool) -> None:
    img = Image.open(src).convert("RGB")
    arr = np.array(img, dtype=np.uint8)
    labels = rgb_to_vascx_labels(arr, tb=tb, th=th, delta=delta, tie_vote=tie_vote)
    dst.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(labels, mode="L").save(dst)


# ---------------- Snakemake entrypoint ----------------

in_dir = Path(str(snakemake.input.in_dir))
out_dir = Path(str(snakemake.output.out_dir))

ext = str(getattr(snakemake.params, "ext", ".png"))
tb = int(getattr(snakemake.params, "tb", 100))
th = int(getattr(snakemake.params, "th", 200))
delta = int(getattr(snakemake.params, "delta", 15))
tie_vote = bool(int(getattr(snakemake.params, "tie_vote", 1)))

files = sorted([p for p in in_dir.glob(f"*{ext}") if p.is_file()])
if not files:
    raise SystemExit(f"No '{ext}' files found in {in_dir}")

out_dir.mkdir(parents=True, exist_ok=True)

for src in files:
    dst = out_dir / (src.stem + ".png")
    convert_one(src, dst, tb=tb, th=th, delta=delta, tie_vote=tie_vote)

print(f"[refined_rgb_to_labels] {len(files)} file(s): {in_dir} -> {out_dir}")

