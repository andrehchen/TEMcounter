#!/usr/bin/env python3
"""Measure nanoparticle diameters in TEM micrographs.

Reads a .tif TEM image of (mesoporous silica) nanoparticles, calibrates the
pixel size from the scale bar burnt into the image, segments the particles and
reports the diameter of N of them in nanometres.

Outputs per image: a CSV of per-particle measurements, an annotated overlay
PNG, and a histogram PNG. A summary row per image is printed and, for multiple
inputs, written to summary.csv.

Usage:
    python temcounter.py IMAGE.tif -n 100
    python temcounter.py folder/ -n 50 -o results/
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass, asdict
from datetime import date
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from skimage import filters, measure, morphology, segmentation
from skimage.feature import peak_local_max

Image.MAX_IMAGE_PIXELS = None

# Segmentation runs on the image resampled to roughly this pixel size, so the
# same parameters behave the same way regardless of magnification.
TARGET_NM_PER_PX = 0.75

IMAGE_SUFFIXES = {".tif", ".tiff", ".png", ".jpg", ".jpeg", ".bmp"}

UNIT_TO_NM = {"nm": 1.0, "um": 1000.0, "µm": 1000.0, "μm": 1000.0, "mm": 1e6}


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------


def load_image(path: Path) -> tuple[np.ndarray, dict]:
    """Return the image as float32 in [0, 1] plus any TIFF description tags."""
    im = Image.open(path)
    meta = {}
    tags = getattr(im, "tag_v2", None)
    if tags is not None:
        for tag in (270, 65500):
            if tag in tags:
                meta[tag] = str(tags[tag])
    arr = np.asarray(im)
    if arr.ndim == 3:  # colour -> luminance
        arr = arr[..., :3].mean(axis=2)
    arr = arr.astype(np.float32)
    hi = float(arr.max())
    if hi > 0:
        arr /= hi
    return arr, meta


def find_info_bar(img: np.ndarray) -> int:
    """Return the first row index of the instrument info bar at the bottom.

    The bar is a solid white strip carrying black text; micrograph rows almost
    never contain saturated-white pixels. Returns img.shape[0] if there is no
    bar (i.e. the whole array is micrograph).
    """
    white = (img >= 0.97).mean(axis=1)
    h = img.shape[0]

    # The longest predominantly white band in the lower part of the frame. Not
    # simply "rows from the bottom up", because some instruments print a dark
    # separator strip below the text.
    runs, start = [], None
    for i in range(int(h * 0.55), h):
        if white[i] >= 0.5 and start is None:
            start = i
        elif white[i] < 0.5 and start is not None:
            runs.append((start, i))
            start = None
    if start is not None:
        runs.append((start, h))
    runs = [r for r in runs if r[1] - r[0] >= 20]
    if not runs:
        return h

    top, bottom = max(runs, key=lambda r: r[1] - r[0])
    # It has to reach nearly the bottom of the frame, and be a strip rather
    # than half the image.
    if bottom < h * 0.85 or h - top > 0.4 * h:
        return h
    return top


# --------------------------------------------------------------------------
# scale bar
# --------------------------------------------------------------------------


@dataclass
class Scale:
    nm_per_px: float
    length_px: float
    length_nm: float
    label: str
    source: str

    def describe(self) -> str:
        return (
            f"{self.nm_per_px:.4f} nm/px  "
            f"(bar {self.length_px:.1f} px = {self.length_nm:g} nm, "
            f"label {self.label!r}, via {self.source})"
        )


def _dark_components(bar: np.ndarray):
    """Label dark (text/ink) components in the info bar."""
    dark = bar < 0.35
    lab, _ = ndi.label(dark)
    return measure.regionprops(lab)


@dataclass
class Ruler:
    span_px: float
    x0: float
    x1: float
    y_bottom: float


def _tick_ruler(props, bar_w: int) -> Ruler | None:
    """Find an evenly spaced row of tick marks (the Hitachi ruler style)."""
    ticks = []
    for p in props:
        minr, minc, maxr, maxc = p.bbox
        h, w = maxr - minr, maxc - minc
        # Tick marks are narrow, taller than wide, and solidly filled.
        if 3 <= h <= 60 and 1 <= w <= 12 and h >= 1.5 * w and p.solidity > 0.7:
            ticks.append((p.centroid[0], p.centroid[1], float(maxr)))
    if len(ticks) < 5:
        return None

    ticks.sort(key=lambda t: t[0])
    best = None
    for i, (y0, _, _) in enumerate(ticks):
        row = [t for t in ticks[i:] if abs(t[0] - y0) <= 4]
        if len(row) < 5:
            continue
        row.sort(key=lambda t: t[1])
        xs = np.array([t[1] for t in row])
        gaps = np.diff(xs)
        med = float(np.median(gaps))
        if med <= 0:
            continue
        # Keep the longest run of consistently spaced ticks.
        run, runs = [0], []
        for j, g in enumerate(gaps, 1):
            if abs(g - med) <= max(2.0, 0.2 * med):
                run.append(j)
            else:
                runs.append(run)
                run = [j]
        runs.append(run)
        run = max(runs, key=len)
        if len(run) < 5:
            continue
        span = xs[run[-1]] - xs[run[0]]
        if best is None or span > best.span_px:
            best = Ruler(
                span_px=float(span),
                x0=float(xs[run[0]]),
                x1=float(xs[run[-1]]),
                y_bottom=max(row[j][2] for j in run),
            )
    if best is None or best.span_px < 0.05 * bar_w:
        return None
    return best


def _solid_bars(props, bar_w: int) -> list[Ruler]:
    """The other common style: a plain solid horizontal bar. Widest first.

    Capped at half the frame width, because a full-width rule separating the
    info bar from the image is the same shape as a scale bar but much longer.
    """
    out = []
    for p in props:
        minr, minc, maxr, maxc = p.bbox
        h, w = maxr - minr, maxc - minc
        if (
            max(40, 0.03 * bar_w) <= w <= 0.5 * bar_w
            and w >= 6 * h
            and p.solidity > 0.85
        ):
            out.append(Ruler(float(w), float(minc), float(maxc), float(maxr)))
    return sorted(out, key=lambda r: -r.span_px)


_TESSERACT_HINT = (
    "Tesseract OCR is needed to read the scale-bar label but was not found.\n"
    "  macOS:   brew install tesseract\n"
    "  Windows: install from https://github.com/UB-Mannheim/tesseract/wiki\n"
    "           (then either tick 'Add to PATH' during setup, or pass\n"
    '            --tesseract "C:\\Program Files\\Tesseract-OCR\\tesseract.exe")\n'
    "  Linux:   sudo apt install tesseract-ocr\n"
    "Alternatively pass --nm-per-px to skip the scale bar and set the scale by hand."
)

# Standard install locations, so a Windows user who did not tick "Add to PATH"
# still works without extra arguments.
_TESSERACT_GUESSES = (
    r"C:\Program Files\Tesseract-OCR\tesseract.exe",
    r"C:\Program Files (x86)\Tesseract-OCR\tesseract.exe",
    "/opt/homebrew/bin/tesseract",
    "/usr/local/bin/tesseract",
    "/usr/bin/tesseract",
)


def locate_tesseract(explicit: str | None = None) -> str | None:
    """Return a usable tesseract executable path, or None."""
    import shutil

    if explicit:
        return explicit if Path(explicit).exists() else None
    found = shutil.which("tesseract")
    if found:
        return found
    for guess in _TESSERACT_GUESSES:
        if Path(guess).exists():
            return guess
    return None


def _ocr(crop: np.ndarray, text_height: float | None = None) -> str:
    """OCR a label crop. Tries single-line then block mode."""
    import pytesseract

    pil = Image.fromarray((np.clip(crop, 0, 1) * 255).astype(np.uint8))
    # Tesseract wants glyphs around 40 px tall; scale from the actual text
    # height when we know it, otherwise from the crop.
    ref = text_height if text_height else crop.shape[0]
    scale = int(np.clip(round(40 / max(ref, 1)), 1, 8))
    if scale > 1:
        pil = pil.resize((pil.width * scale, pil.height * scale), Image.LANCZOS)

    whitelist = "-c tessedit_char_whitelist=0123456789.numµμ"
    for psm in (7, 6):
        try:
            text = pytesseract.image_to_string(
                pil, config=f"--psm {psm} {whitelist}"
            ).strip()
        except Exception:
            return ""
        if _parse_label(text):
            return text
    return text


_LABEL_RE = re.compile(r"(\d+(?:\.\d+)?)\s*(nm|um|µm|μm|mm)", re.I)


def _parse_label(text: str) -> tuple[float, str] | None:
    m = _LABEL_RE.search(text.replace(" ", ""))
    if not m:
        return None
    value = float(m.group(1))
    unit = m.group(2).lower()
    return value * UNIT_TO_NM[unit], m.group(0)


def read_scale(img: np.ndarray, bar_top: int, verbose: bool = False) -> Scale | None:
    """Calibrate nm/px from the burnt-in scale bar."""
    if bar_top >= img.shape[0]:
        return None
    bar = img[bar_top:]
    props = _dark_components(bar)

    candidates = []
    ticks = _tick_ruler(props, bar.shape[1])
    if ticks is not None:
        candidates.append(("tick ruler", ticks))
    candidates += [("solid bar", r) for r in _solid_bars(props, bar.shape[1])[:3]]

    # Try each candidate bar until one has a label next to it that parses as a
    # length. A candidate with no readable label is not a scale bar.
    for kind, ruler in candidates:
        for crop, box, text_h in _label_crops(bar, props, ruler):
            text = _ocr(crop, text_h)
            if verbose:
                print(
                    f"    {kind} span={ruler.span_px:.0f}px, OCR {box} -> {text!r}",
                    file=sys.stderr,
                )
            parsed = _parse_label(text)
            if parsed:
                length_nm, label = parsed
                return Scale(
                    length_nm / ruler.span_px,
                    ruler.span_px,
                    length_nm,
                    label,
                    f"scale bar ({kind}) + OCR",
                )
    return None


def _label_crops(bar: np.ndarray, props, ruler: Ruler):
    """Yield candidate crops of the scale-bar label, most likely first.

    The label is a line of text near the bar -- usually underneath it, but some
    instruments put it above. So group the nearby dark components into text
    lines and offer them up nearest-first, rather than assuming a position.
    """
    h, w = bar.shape
    pad = 0.4 * ruler.span_px
    glyphs = [
        p
        for p in props
        if ruler.x0 - pad <= p.centroid[1] <= ruler.x1 + pad
        and 4 <= p.bbox[2] - p.bbox[0] <= 120
        and not (p.bbox[0] < ruler.y_bottom <= p.bbox[2])  # not the bar itself
    ]

    if glyphs:
        med_h = float(np.median([p.bbox[2] - p.bbox[0] for p in glyphs]))
        lines: list[list] = []
        for p in sorted(glyphs, key=lambda p: p.bbox[0]):
            for line in lines:
                if abs(p.bbox[0] - line[0].bbox[0]) <= 0.7 * med_h:
                    line.append(p)
                    break
            else:
                lines.append([p])

        # Nearest to the bar first, preferring the line below it.
        def distance(line):
            top = min(p.bbox[0] for p in line)
            below = top >= ruler.y_bottom
            return (0 if below else 1, abs(top - ruler.y_bottom))

        for line in sorted(lines, key=distance)[:4]:
            r0 = min(p.bbox[0] for p in line) - 4
            r1 = max(p.bbox[2] for p in line) + 4
            c0 = min(p.bbox[1] for p in line) - 8
            c1 = max(p.bbox[3] for p in line) + 8
            out = _clip(bar, r0, r1, c0, c1)
            if out is not None:
                yield out[0], out[1], med_h

    yb = int(ruler.y_bottom)
    for box in (
        (yb + 2, yb + 110, int(ruler.x0 - pad), int(ruler.x1 + pad)),
        (yb + 2, h, int(ruler.x0 - 2 * pad), w),
        (0, h, int(ruler.x0 - 2 * pad), w),
    ):
        out = _clip(bar, *box)
        if out is not None:
            yield out[0], out[1], None


def _clip(bar: np.ndarray, r0, r1, c0, c1):
    h, w = bar.shape
    r0, c0 = max(0, int(r0)), max(0, int(c0))
    r1, c1 = min(h, int(r1)), min(w, int(c1))
    if r1 - r0 < 5 or c1 - c0 < 5:
        return None
    crop = bar[r0:r1, c0:c1]
    if (crop < 0.35).sum() < 20:
        return None
    return crop, (r0, r1, c0, c1)


def scale_from_metadata(meta: dict, img_w: int) -> Scale | None:
    """Last-resort fallback for Hitachi HT7700 TIFFs (see README caveat)."""
    text = " ".join(meta.values())
    m = re.search(r"Magnification=(\d+)", text)
    if not m:
        return None
    mag = float(m.group(1))
    if mag <= 0:
        return None
    # Empirical constant from the scale bars of this instrument/camera pair.
    nm_per_px = 10949.0 / mag
    return Scale(nm_per_px, float("nan"), float("nan"), f"x{mag:g}", "TIFF metadata")


# --------------------------------------------------------------------------
# segmentation
# --------------------------------------------------------------------------


@dataclass
class Particle:
    id: int
    x_px: float
    y_px: float
    diameter_nm: float
    major_nm: float
    minor_nm: float
    area_nm2: float
    circularity: float
    roundness: float
    solidity: float
    aspect_ratio: float
    cluster_size: int


def _drop_small(bw: np.ndarray, min_area: int) -> np.ndarray:
    """Remove connected components smaller than min_area pixels."""
    lab, n = ndi.label(bw)
    if n == 0:
        return bw
    keep = np.bincount(lab.ravel()) >= min_area
    keep[0] = False
    return keep[lab]


def segment(
    img: np.ndarray,
    nm_per_px: float,
    *,
    dark_particles: bool = True,
    min_nm: float | None = None,
    max_nm: float | None = None,
    min_circularity: float = 0.80,
    min_roundness: float = 0.85,
    min_solidity: float = 0.90,
    max_aspect: float = 1.35,
    max_cluster: int = 1,
    threshold_scale: float = 1.0,
    keep_edge: bool = False,
    verbose: bool = False,
):
    """Segment particles; return (particles, label image, downsample factor)."""
    ds = max(1, int(round(TARGET_NM_PER_PX / nm_per_px)))
    work = img[::ds, ::ds]
    px = nm_per_px * ds  # nm per working pixel

    if not dark_particles:
        work = 1.0 - work

    # Flatten the illumination/thickness gradient, then blur away the pore
    # texture so a particle reads as one solid blob.
    bg = filters.gaussian(work, sigma=max(8.0, 40.0 / ds), preserve_range=True)
    flat = work / np.maximum(bg, 1e-6)
    smooth = filters.gaussian(flat, sigma=max(1.5, 2.5 / px))

    thresh = filters.threshold_otsu(smooth) * threshold_scale
    bw = smooth < thresh
    bw = ndi.binary_fill_holes(bw)

    min_area_px = 30.0
    if min_nm:
        min_area_px = math.pi * (min_nm / 2.0 / px) ** 2 * 0.5
    bw = _drop_small(bw, int(max(20, min_area_px)))
    bw = morphology.opening(bw, morphology.disk(max(1, int(round(2.0 / px)))))
    if not bw.any():
        return [], np.zeros_like(bw, dtype=np.int32), ds, 0

    dist = ndi.distance_transform_edt(bw)

    # Estimate a typical particle radius from the distance transform so the
    # watershed seed spacing adapts to the particle size in this image.
    p90 = float(np.percentile(dist[bw], 90))
    radius = max(2.0, p90 / 0.68)
    if min_nm:
        radius = max(radius, min_nm / 2.0 / px)
    min_distance = int(max(3, round(0.7 * radius)))

    # Without an explicit floor, reject anything under a third of the typical
    # particle size in this image: stain specks and debris are small but
    # perfectly round, so no shape filter catches them.
    if min_nm is None:
        min_nm = 0.35 * 2 * radius * px
        if verbose:
            print(f"    auto minimum diameter {min_nm:.1f} nm", file=sys.stderr)

    if verbose:
        print(
            f"    ds={ds} px={px:.3f} nm  est. radius {radius * px:.1f} nm  "
            f"seed spacing {min_distance} px",
            file=sys.stderr,
        )

    coords = peak_local_max(
        dist, min_distance=min_distance, labels=bw, exclude_border=False
    )
    markers = np.zeros(dist.shape, dtype=np.int32)
    if len(coords):
        markers[tuple(coords.T)] = np.arange(1, len(coords) + 1)
    labels = segmentation.watershed(-dist, markers, mask=bw)

    # How many particles share each blob of touching material. A blob holding
    # more than one is an aggregate (or at least particles in contact), and its
    # members are measured across a contact line rather than a free edge.
    cluster_size = _cluster_sizes(labels, bw)

    if not keep_edge:
        labels = segmentation.clear_border(labels)

    particles = []
    n_aggregated = 0
    kept = np.zeros_like(labels)
    for prop in measure.regionprops(labels):
        area_px = float(prop.area)
        d_px = 2.0 * math.sqrt(area_px / math.pi)
        d_nm = d_px * px
        if min_nm and d_nm < min_nm:
            continue
        if max_nm and d_nm > max_nm:
            continue
        perim = float(prop.perimeter_crofton) or 1.0
        circ = min(1.0, 4.0 * math.pi * area_px / (perim * perim))
        round_ = _roundness(prop.image, d_px)
        major = float(prop.axis_major_length) * px
        minor = float(prop.axis_minor_length) * px
        ar = (major / minor) if minor > 0 else float("inf")
        if (
            circ < min_circularity
            or round_ < min_roundness
            or prop.solidity < min_solidity
            or ar > max_aspect
        ):
            continue
        csize = int(cluster_size[prop.label])
        if max_cluster and csize > max_cluster:
            n_aggregated += 1
            continue
        particles.append(
            Particle(
                id=0,
                x_px=float(prop.centroid[1]) * ds,
                y_px=float(prop.centroid[0]) * ds,
                diameter_nm=d_nm,
                major_nm=major,
                minor_nm=minor,
                area_nm2=area_px * px * px,
                circularity=circ,
                roundness=round_,
                solidity=float(prop.solidity),
                aspect_ratio=ar,
                cluster_size=csize,
            )
        )
        kept[labels == prop.label] = len(particles)

    return particles, kept, ds, n_aggregated


def _cluster_sizes(labels: np.ndarray, bw: np.ndarray) -> np.ndarray:
    """Particles per blob of touching material, indexed by watershed label."""
    n = int(labels.max())
    out = np.ones(n + 1, dtype=int)
    if n == 0:
        return out
    blobs, n_blobs = ndi.label(bw)
    idx = np.arange(1, n + 1)
    # Every pixel of a watershed region lies in exactly one blob, so the max
    # blob id over the region is that blob's id.
    parent = np.nan_to_num(
        ndi.maximum(blobs, labels=labels, index=idx), nan=0
    ).astype(int)
    per_blob = np.bincount(parent, minlength=n_blobs + 1)
    out[idx] = per_blob[parent]
    return out


def _roundness(mask: np.ndarray, equiv_d_px: float) -> float:
    """Diameter of the largest inscribed circle over the equivalent diameter.

    1.0 for a perfect disk and ~0.7 for the two shapes that most often survive
    the other filters: a merged pair of touching particles, and the crescent
    left when the watershed cuts one particle in the wrong place.
    """
    if equiv_d_px <= 0:
        return 0.0
    padded = np.pad(mask, 1)
    inscribed = float(ndi.distance_transform_edt(padded).max())
    return min(1.0, 2.0 * inscribed / equiv_d_px)


# --------------------------------------------------------------------------
# selection and reporting
# --------------------------------------------------------------------------


def select(particles, n: int | None, mode: str, seed: int):
    if n is None or n >= len(particles):
        chosen = list(particles)
    elif mode == "random":
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(len(particles), size=n, replace=False))
        chosen = [particles[i] for i in idx]
    elif mode == "best":
        chosen = sorted(particles, key=lambda p: -p.circularity)[:n]
    elif mode == "first":
        chosen = sorted(particles, key=lambda p: (p.y_px, p.x_px))[:n]
    else:
        raise ValueError(f"unknown selection mode {mode!r}")
    chosen.sort(key=lambda p: (p.y_px, p.x_px))
    for i, p in enumerate(chosen, 1):
        p.id = i
    return chosen


def stats(particles) -> dict:
    d = np.array([p.diameter_nm for p in particles], dtype=float)
    if d.size == 0:
        return {"n": 0}
    out = {
        "n": int(d.size),
        "mean_nm": float(d.mean()),
        "sd_nm": float(d.std(ddof=1)) if d.size > 1 else 0.0,
        "cv_pct": float(d.std(ddof=1) / d.mean() * 100) if d.size > 1 else 0.0,
        "median_nm": float(np.median(d)),
        "min_nm": float(d.min()),
        "max_nm": float(d.max()),
        "normal_p": None,
    }
    # Is a normal a fair description of this distribution? Shapiro-Wilk needs a
    # handful of points to say anything useful.
    if d.size >= 8:
        from scipy import stats as sps

        out["normal_p"] = float(sps.shapiro(d).pvalue)
    return out


def print_histogram(particles, bins: int = 12, width: int = 40, indent: str = "  ") -> None:
    """Draw the size distribution in the terminal, with the normal fit marked."""
    d = np.array([p.diameter_nm for p in particles], dtype=float)
    if d.size < 2 or d.max() == d.min():
        return
    counts, edges = np.histogram(d, bins=bins)
    peak = int(counts.max())
    if peak == 0:
        return
    mu, sigma = float(d.mean()), float(d.std(ddof=1))
    bin_w = edges[1] - edges[0]

    print(f"{indent}{'diameter (nm)':>13}  {'':{width}}  count")
    for i, n in enumerate(counts):
        bar = ["#"] * round(n / peak * width)
        bar += [" "] * (width - len(bar))
        if sigma > 0:
            centre = (edges[i] + edges[i + 1]) / 2
            expected = (
                math.exp(-0.5 * ((centre - mu) / sigma) ** 2)
                / (sigma * math.sqrt(2 * math.pi))
                * d.size
                * bin_w
            )
            # Clamp so the marker still shows when the fit overshoots the
            # tallest bar, rather than silently falling off the end.
            col = min(width - 1, max(0, round(expected / peak * width)))
            bar[col] = "|"
        label = f"{edges[i]:5.1f}-{edges[i + 1]:5.1f}"
        print(f"{indent}{label:>13}  {''.join(bar)}  {n}")
    print(f"{indent}{'':>13}  # measured   | normal fit (mu={mu:.1f}, sigma={sigma:.1f})")


def write_particles_csv(path: Path, rows: list[dict]) -> None:
    """Every particle from the whole run, one file, with an image column."""
    if not rows:
        return
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for row in rows:
            w.writerow(
                {k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}
            )


def save_overlay(path: Path, img: np.ndarray, labels, ds: int, particles, scale) -> None:
    from PIL import ImageDraw

    view = img[::ds, ::ds]
    lo, hi = np.percentile(view, [0.5, 99.5])
    view = np.clip((view - lo) / max(hi - lo, 1e-6), 0, 1)
    rgb = np.repeat((view * 255).astype(np.uint8)[..., None], 3, axis=2)

    keep_ids = {p.id for p in particles}
    mask = np.isin(labels, list(keep_ids)) if keep_ids else np.zeros_like(labels, bool)
    edges = mask ^ morphology.erosion(mask, morphology.disk(2))
    rgb[edges] = (255, 40, 40)

    other = (labels != 0) & ~mask
    other_edges = other ^ morphology.erosion(other, morphology.disk(1))
    rgb[other_edges] = (70, 150, 255)

    pil = Image.fromarray(rgb)
    draw = ImageDraw.Draw(pil)
    for p in particles:
        draw.text((p.x_px / ds + 4, p.y_px / ds - 6), str(p.id), fill=(255, 220, 0))

    # Burn in a fresh scale bar as a visual check on the calibration.
    if scale is not None:
        bar_nm = 100.0
        while bar_nm / scale.nm_per_px / ds < 0.08 * pil.width:
            bar_nm *= 2
        length = bar_nm / scale.nm_per_px / ds
        x1, y1 = pil.width - 30, pil.height - 30
        draw.rectangle([x1 - length, y1 - 8, x1, y1], fill=(255, 255, 0))
        draw.text((x1 - length, y1 - 26), f"{bar_nm:g} nm", fill=(255, 255, 0))

    pil.save(path)


def save_histogram(path: Path, particles, title: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = np.array([p.diameter_nm for p in particles], dtype=float)
    s = stats(particles)
    fig, ax = plt.subplots(figsize=(6.5, 4.2), dpi=150)
    bins = max(8, min(30, int(np.sqrt(d.size) * 1.5)))
    counts, edges, _ = ax.hist(
        d, bins=bins, color="#4C72B0", edgecolor="white", label="measured"
    )

    mu, sigma = s["mean_nm"], s["sd_nm"]
    if sigma > 0 and d.size > 1:
        # Scale the PDF to counts so it sits on the same axis as the bars.
        width = edges[1] - edges[0]
        x = np.linspace(edges[0] - width, edges[-1] + width, 400)
        pdf = np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * math.sqrt(2 * math.pi))
        ax.plot(
            x,
            pdf * d.size * width,
            color="#C44E52",
            lw=2,
            label=f"normal fit\nµ={mu:.1f}, σ={sigma:.1f} nm",
        )
    ax.axvline(mu, color="#C44E52", lw=1, ls="--")
    ax.set_xlabel("Diameter (nm)")
    ax.set_ylabel("Count")
    subtitle = (
        f"n={s['n']}  mean={mu:.1f} ± {sigma:.1f} nm (CV {s['cv_pct']:.1f}%)"
    )
    if s.get("normal_p") is not None:
        verdict = "consistent with normal" if s["normal_p"] >= 0.05 else "not normal"
        subtitle += f"\nShapiro-Wilk p={s['normal_p']:.3f} — {verdict}"
    ax.set_title(f"{title}\n{subtitle}", fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# --------------------------------------------------------------------------


@dataclass
class RunDirs:
    """Where one run's output goes."""

    root: Path
    name: str

    @property
    def summary(self) -> Path:
        return self.root / f"{self.name}.csv"

    @property
    def particles(self) -> Path:
        return self.root / f"{self.name}_particles.csv"

    @property
    def overlays(self) -> Path:
        return self.root / "overlays"

    @property
    def histograms(self) -> Path:
        return self.root / "histograms"

    def create(self) -> None:
        for d in (self.root, self.overlays, self.histograms):
            d.mkdir(parents=True, exist_ok=True)


def default_run_name(outdir: Path) -> str:
    """YYMMDD_SummaryN, counting up past any runs already in outdir today."""
    stamp = date.today().strftime("%y%m%d")
    n = 1
    while (outdir / f"{stamp}_Summary{n}").exists():
        n += 1
    return f"{stamp}_Summary{n}"


def clean_run_name(raw: str, fallback: str) -> str:
    name = raw.strip().strip("\"'")
    if name.lower().endswith(".csv"):
        name = name[:-4]
    # Characters Windows will not accept in a file name.
    name = re.sub(r'[<>:"/\\|?*]', "_", name).strip(" .")
    return name or fallback


def ask_run_name(outdir: Path, preset: str | None) -> str:
    fallback = default_run_name(outdir)
    if preset is not None:
        name = clean_run_name(preset, fallback)
    elif sys.stdin.isatty():
        try:
            name = clean_run_name(input(f"Name this run [{fallback}]: "), fallback)
        except EOFError:
            name = fallback
    else:
        name = fallback

    # Never silently overwrite a previous run the user named themselves.
    if (outdir / name).exists() and name != fallback:
        base, n = name, 2
        while (outdir / f"{base}_{n}").exists():
            n += 1
        name = f"{base}_{n}"
        print(f"that name is taken, using {name} instead")
    return name


def gather_inputs(paths) -> list[Path]:
    out = []
    for raw in paths:
        p = Path(raw)
        if p.is_dir():
            out += sorted(
                q for q in p.iterdir() if q.suffix.lower() in IMAGE_SUFFIXES
            )
        else:
            out.append(p)
    return out


def process(path: Path, args, dirs: RunDirs) -> dict | None:
    print(f"\n{path.name}")
    img, meta = load_image(path)
    bar_top = find_info_bar(img)
    micrograph = img[:bar_top]
    print(f"  image {img.shape[1]}x{img.shape[0]} px, micrograph height {bar_top} px")

    if args.nm_per_px:
        scale = Scale(args.nm_per_px, float("nan"), float("nan"), "manual", "--nm-per-px")
    else:
        scale = read_scale(img, bar_top, verbose=args.verbose)
        if scale is None:
            scale = scale_from_metadata(meta, img.shape[1])
            if scale is not None:
                print(
                    "  WARNING: could not read the scale bar. Falling back to the\n"
                    "  TIFF magnification field and a constant fitted to one "
                    "particular\n"
                    "  HT7700 + XR81-B camera. If that is not your instrument these\n"
                    "  numbers are wrong -- check the bar in the overlay PNG, or "
                    "pass\n"
                    "  --nm-per-px explicitly.",
                    file=sys.stderr,
                )
    if scale is None:
        print("  ERROR: could not determine the scale; pass --nm-per-px", file=sys.stderr)
        return None
    print(f"  scale: {scale.describe()}")

    particles, labels, ds, n_aggregated = segment(
        micrograph,
        scale.nm_per_px,
        dark_particles=not args.bright_particles,
        min_nm=args.min_diameter,
        max_nm=args.max_diameter,
        min_circularity=args.min_circularity,
        min_roundness=args.min_roundness,
        min_solidity=args.min_solidity,
        max_aspect=args.max_aspect,
        max_cluster=args.max_cluster,
        threshold_scale=args.threshold_scale,
        keep_edge=args.keep_edge,
        verbose=args.verbose,
    )
    note = ""
    if args.max_cluster:
        limit = (
            "touching another particle"
            if args.max_cluster == 1
            else f"in a group of more than {args.max_cluster}"
        )
        note = f" ({n_aggregated} more rejected as {limit})"
    print(f"  {len(particles)} particles passed the filters{note}")

    chosen = select(particles, args.n, args.select, args.seed)
    if args.n and len(chosen) < args.n:
        print(
            f"  WARNING: only {len(chosen)} usable particles, fewer than the "
            f"{args.n} requested"
        )
    if not chosen:
        print("  ERROR: nothing to measure", file=sys.stderr)
        return None

    # Re-index the label image to the chosen particles for the overlay.
    keep = np.zeros_like(labels)
    lookup = {}
    for old, p in enumerate(particles, 1):
        lookup[old] = p.id if p in chosen else 0
    remap = np.zeros(labels.max() + 1, dtype=labels.dtype)
    for old, new in lookup.items():
        remap[old] = new
    keep = remap[labels]
    rejected = (labels > 0) & (keep == 0)
    keep[rejected] = -1

    stem = path.stem.replace(" ", "_")
    save_overlay(
        dirs.overlays / f"{stem}_overlay.png", micrograph, keep, ds, chosen, scale
    )
    save_histogram(dirs.histograms / f"{stem}_hist.png", chosen, path.name)

    s = stats(chosen)
    print(
        f"  diameter: {s['mean_nm']:.1f} +/- {s['sd_nm']:.1f} nm  "
        f"(median {s['median_nm']:.1f}, range {s['min_nm']:.1f}-{s['max_nm']:.1f}, "
        f"CV {s['cv_pct']:.1f}%, n={s['n']})"
    )
    if s.get("normal_p") is not None:
        verdict = "consistent with normal" if s["normal_p"] >= 0.05 else "not normal"
        print(
            f"  normal fit: mu={s['mean_nm']:.1f} nm, sigma={s['sd_nm']:.1f} nm  "
            f"(Shapiro-Wilk p={s['normal_p']:.3f}, {verdict})"
        )
    if args.hist_bins:
        print_histogram(chosen, bins=args.hist_bins)

    summary = {
        "image": path.name,
        "nm_per_px": round(scale.nm_per_px, 5),
        "scale_source": scale.source,
        "scale_label": scale.label,
        "detected": len(particles),
        "aggregated_rejected": n_aggregated,
        **{k: (round(v, 3) if isinstance(v, float) else v) for k, v in s.items()},
    }
    rows = [{"image": path.name, **asdict(p)} for p in chosen]
    return summary, rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Measure nanoparticle diameters in TEM images using the "
        "burnt-in scale bar for calibration."
    )
    ap.add_argument("images", nargs="+", help=".tif file(s) or a folder of them")
    ap.add_argument("-n", type=int, default=None, help="number of particles to measure")
    ap.add_argument("-o", "--outdir", type=Path, default=Path("results"))
    ap.add_argument(
        "--name",
        default=None,
        help="name for this run's folder and summary CSV. Omitted, you are "
        "prompted for one, and pressing Enter accepts YYMMDD_SummaryN",
    )
    ap.add_argument(
        "--select",
        choices=["random", "best", "first"],
        default="best",
        help="which N to keep when more are found: best = the most circular "
        "(default), random = an unbiased seeded sample, first = top-left first",
    )
    ap.add_argument("--seed", type=int, default=0, help="seed for --select random")
    ap.add_argument("--nm-per-px", type=float, default=None, help="override the scale")
    ap.add_argument(
        "--tesseract",
        default=None,
        help="path to tesseract.exe, if it is not on PATH (Windows)",
    )
    ap.add_argument(
        "--min-diameter",
        type=float,
        default=None,
        help="reject particles below this, in nm. Left unset, a floor of a "
        "third of the typical particle size in each image is applied "
        "automatically to drop specks and debris; pass 0 to disable",
    )
    ap.add_argument("--max-diameter", type=float, default=None, help="reject above, nm")
    ap.add_argument("--min-circularity", type=float, default=0.80)
    ap.add_argument(
        "--min-roundness",
        type=float,
        default=0.85,
        help="inscribed-circle diameter / equivalent diameter; rejects merged "
        "pairs and mis-split fragments (default: 0.85)",
    )
    ap.add_argument("--min-solidity", type=float, default=0.90)
    ap.add_argument("--max-aspect", type=float, default=1.35)
    ap.add_argument(
        "--max-cluster",
        type=int,
        default=1,
        help="largest group of touching particles still counted. 1 (default) "
        "measures only free-standing particles and so excludes aggregates; "
        "0 disables the check",
    )
    ap.add_argument(
        "--threshold-scale",
        type=float,
        default=1.0,
        help="nudge the Otsu threshold to tune where the particle edge falls; "
        "+1%% grows the mean diameter by roughly 1%% (default: 1.0)",
    )
    ap.add_argument(
        "--keep-edge", action="store_true", help="keep particles touching the border"
    )
    ap.add_argument(
        "--bright-particles",
        action="store_true",
        help="particles are brighter than the background",
    )
    ap.add_argument(
        "--hist-bins",
        type=int,
        default=12,
        help="bins in the terminal histogram; 0 to skip it (default: 12)",
    )
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    paths = gather_inputs(args.images)
    if not paths:
        print("no input images found", file=sys.stderr)
        return 1

    # Fail up front rather than silently degrading to the metadata fallback,
    # which is calibrated to one specific microscope.
    if not args.nm_per_px:
        exe = locate_tesseract(args.tesseract)
        if exe is None:
            print(_TESSERACT_HINT, file=sys.stderr)
            return 1
        import pytesseract

        pytesseract.pytesseract.tesseract_cmd = exe
        if args.verbose:
            print(f"using tesseract at {exe}", file=sys.stderr)

    dirs = RunDirs(args.outdir / ask_run_name(args.outdir, args.name), "")
    dirs.name = dirs.root.name
    dirs.create()

    rows, all_particles = [], []
    for p in paths:
        if not p.exists():
            print(f"missing: {p}", file=sys.stderr)
            continue
        result = process(p, args, dirs)
        if result:
            summary, particle_rows = result
            rows.append(summary)
            all_particles += particle_rows

    if not rows:
        return 1

    # Pool every particle in the run into one extra summary row, so a folder of
    # repeats of the same sample gives its combined statistics directly.
    out_rows = list(rows)
    if len(rows) > 1:
        d = np.array([r["diameter_nm"] for r in all_particles], dtype=float)
        pooled = {
            "n": int(d.size),
            "mean_nm": float(d.mean()),
            "sd_nm": float(d.std(ddof=1)),
            "cv_pct": float(d.std(ddof=1) / d.mean() * 100),
            "median_nm": float(np.median(d)),
            "min_nm": float(d.min()),
            "max_nm": float(d.max()),
        }
        out_rows.append(
            {
                "image": f"POOLED ({len(rows)} images)",
                **{k: round(v, 3) if isinstance(v, float) else v
                   for k, v in pooled.items()},
            }
        )

    fields = list(rows[0].keys())
    with dirs.summary.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, restval="")
        w.writeheader()
        w.writerows(out_rows)
    write_particles_csv(dirs.particles, all_particles)

    print(f"\n{len(rows)} image(s), {len(all_particles)} particles measured")
    if len(rows) > 1:
        print(
            f"  pooled: {pooled['mean_nm']:.1f} +/- {pooled['sd_nm']:.1f} nm "
            f"(CV {pooled['cv_pct']:.1f}%, n={pooled['n']})"
        )
    w = max(len(dirs.summary.name), len(dirs.particles.name), 12) + 2
    print(f"  {dirs.root}/")
    print(f"    {dirs.summary.name:<{w}}one row per image, plus a pooled row")
    print(f"    {dirs.particles.name:<{w}}every particle measured")
    print(f"    {'overlays/':<{w}}check these: outlines drawn on the micrograph")
    print(f"    {'histograms/':<{w}}size distribution with the normal fit")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
