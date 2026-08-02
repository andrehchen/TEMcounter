# TEMcounter

Measures nanoparticle diameters in TEM micrographs. Point it at a `.tif`, say how
many particles you want, and it calibrates the pixel size from the scale bar burnt
into the image, segments the particles, and reports diameters in nanometres.

## Setup

```sh
python3 -m venv .venv
.venv/bin/pip install numpy scipy scikit-image pillow matplotlib pytesseract
brew install tesseract        # needed to read the scale-bar label
```

## Usage

```sh
.venv/bin/python temcounter.py "OneDrive_1_2026-07-23/B1 AE_1.tif" -n 100
.venv/bin/python temcounter.py OneDrive_1_2026-07-23 -n 100 -o results
```

For each image it prints the scale it read, the summary statistics, and the
distribution as a terminal histogram, with `|` marking where the normal fit
predicts each bar should land:

```
  diameter: 30.1 +/- 4.3 nm  (median 29.9, range 18.5-39.4, CV 14.3%, n=100)
  normal fit: mu=30.1 nm, sigma=4.3 nm  (Shapiro-Wilk p=0.472, consistent with normal)
  diameter (nm)                                            count
     25.5- 27.2  ###########################|              11
     27.2- 29.0  ####################################|###  16
     29.0- 30.7  ######################################|   15
                 # measured   | normal fit (mu=30.1, sigma=4.3)
```

`--hist-bins N` changes its resolution and `--hist-bins 0` turns it off.

## Output

Each run gets its own folder, so runs never overwrite each other. Before it starts
you are asked to name the run:

```
Name this run [260728_Summary1]:
```

Press Enter to accept the default — `YYMMDD_SummaryN`, where N counts up past any
runs already made today (`260728_Summary1`, `260728_Summary2`, …). Or type a name;
characters Windows will not accept in a filename are replaced, and an existing name
gets a numeric suffix rather than being overwritten. `--name` sets it directly and
skips the prompt, which is also what happens when output is piped.

The CSVs are per run, not per image, so processing a folder gives two of them
however many images went in:

```
results/260728_Summary1/
    260728_Summary1.csv            one row per image, plus a pooled row
    260728_Summary1_particles.csv  every particle measured, with an image column
    overlays/                      one PNG per image
    histograms/                    one PNG per image
```

- **`<run>.csv`** — per image: the calibration it read, how many particles were
  detected and how many were rejected as aggregated, then n, mean, SD, CV, median,
  range and the Shapiro-Wilk p. When there is more than one image, a final
  `POOLED` row gives the statistics over every particle in the run combined.
- **`<run>_particles.csv`** — one row per particle across the whole run: diameter,
  major/minor axis, area, shape metrics, centroid, cluster size. Filter on the
  `image` column to get back to a single micrograph.
- **`overlays/`** — the micrograph with measured particles outlined in red and
  numbered, other detections in blue, and a scale bar redrawn from the derived
  calibration as a visual check.
- **`histograms/`** — diameter distribution with the fitted normal curve.

`diameter_nm` is the equivalent circle diameter, `2·√(A/π)` — the standard measure
for roughly spherical particles. `major_nm`/`minor_nm` are the fitted ellipse axes
if you want the shape as well.

The `POOLED` row combines every particle in the run, so it is only meaningful when
the images are repeats of the *same* sample. Pointing the tool at a folder holding
two different samples pools both — the per-image rows are the ones to read there.

Useful options:

- `--select best|random|first` — which N to keep when more particles are found.
  Default `best`, the most circular ones, which are the cleanest segmentations.
  `random` draws a seeded, reproducible, unbiased sample instead; on the
  reference images the two agree on the mean to within 0.3 nm but `best` narrows
  the SD by 10-20%, so it understates the true spread (see *Accuracy*).
- `--max-cluster` — largest group of touching particles still counted. Default
  `1`, i.e. only free-standing particles are measured (see *Aggregates* below).
  `--max-cluster 0` turns the check off.
- `--min-diameter` / `--max-diameter` in nm, to exclude debris or aggregates.
- `--nm-per-px` to override the calibration entirely.
- `--keep-edge` to include particles clipped by the image border (off by default,
  since their diameter is not measurable).
- `--bright-particles` if the particles are brighter than the background.
- `--threshold-scale` to tune where the particle edge is drawn (see below).
- `-v` to print the OCR attempts and the segmentation parameters chosen.

## How the scale is read

The instrument info bar is found as the longest predominantly white band in the
lower part of the frame — not simply the rows up from the bottom, because some
instruments print a dark separator strip below the text.

Within it, two scale-bar styles are recognised: an evenly spaced row of short
vertical tick marks (Hitachi), and a plain solid bar (AMT). Solid candidates are
capped at half the frame width, since a full-width rule separating the info bar
from the image has exactly the same shape as a scale bar, only longer.

The label is a line of text near the bar. Rather than assuming it sits underneath,
the nearby dark components are grouped into text lines and offered up
nearest-first, so a label above or beside the bar is still found. Each is read
with Tesseract — single-line mode first, then block mode — restricted to the
characters that can appear in a length (`0123456789.numµ`). A candidate bar with
no readable length next to it is rejected, so the bar and its label validate each
other.

For the reference images this gives:

| magnification | bar | label | nm/px |
|---|---|---|---|
| ×30 000 | 548 px | 200 nm | 0.3650 |
| ×20 000 | 918 px | 500 nm | 0.5447 |
| ×12 000 | 548 px | 500 nm | 0.9124 |
| ×10 000 | 918 px | 1.0 µm | 1.0893 |

These four independent calibrations agree to within 0.5% once scaled by
magnification, which confirms that the label is the full tick-to-tick span rather
than a single division.

On the AMT/XR280 images the derived scale can be checked against the `Cal:` value
the microscope prints in its own info bar: 1.0417 vs 1.0370 nm/px at 15000x and
0.5208 vs 0.5183 nm/px at 30000x, both within 0.5%.

If no scale bar can be found, the tool falls back to the `Magnification` field in
the Hitachi TIFF header using an empirical constant fitted to these images
(`nm/px = 10949 / magnification`) and prints a warning. That constant is specific
to this microscope and camera — treat those numbers as provisional.

## How the particles are measured

1. The image is resampled to ≈0.75 nm/px so the same parameters behave identically
   at any magnification.
2. The illumination/thickness gradient is divided out with a large-sigma Gaussian
   background, then the mesopore texture is blurred away so each particle reads as
   one solid blob.
3. Otsu threshold, hole filling, and a small opening give the particle mask.
4. A distance transform sets an adaptive watershed seed spacing (derived from the
   typical particle radius in that image), and the watershed splits touching
   particles.
5. Anything under a third of the typical particle size in that image is dropped.
   Stain specks and debris are small but perfectly round, so no shape filter
   catches them; `--min-diameter` overrides the automatic floor, `0` disables it.
6. Particles touching the border are dropped, and the rest are filtered on shape:
   circularity, solidity, aspect ratio, and **roundness** — the diameter of the
   largest inscribed circle over the equivalent diameter. Roundness is the useful
   one here: it is 1.0 for a disk but falls to ≈0.7 both for a merged pair of
   touching particles and for the crescent left when the watershed cuts a particle
   in the wrong place, which are the two artefacts that survive the other filters.

Always look at the overlay PNG. It is the fastest way to confirm the outlines sit
where you would have drawn them.

## The normal fit

The histogram carries a fitted normal curve, scaled to counts so it overlays the
bars directly. µ and σ are the sample mean and SD — for a normal these *are* the
maximum-likelihood fit, so there is no separate curve-fitting step and the numbers
match the ones printed to the console and stored in `summary.csv`.

Alongside it is a Shapiro-Wilk test of whether a normal is a fair description of
the data at all, reported as `normal_p` (computed when n ≥ 8). p ≥ 0.05 means the
data is consistent with a normal; below that, the curve is still drawn but you
should not lean on it. On the reference images 7 of 10 pass; the three that fail
(p = 0.030–0.043) are marginal, and are all N42, whose distributions are slightly
right-skewed. Nanoparticle sizes are often closer to log-normal than normal, so if
you see low p-values consistently, that is the likely reason.

## Aggregates

Aggregates are excluded by default, in two stages.

Fused or overlapping clumps that segment as one lumpy region are thrown out by the
shape filters, roundness in particular. What that alone does *not* catch is a
particle in a cluster that the watershed splits cleanly: it looks round, so it
passes, but part of its outline is a contact line with its neighbour rather than a
free edge.

So each particle also carries a `cluster_size`: the number of particles sharing its
blob of touching material. `--max-cluster 1` (the default) keeps only particles
whose blob contains nothing else, so anything in contact with a neighbour is
dropped regardless of how well it segmented. The `cluster_size` column is in the
CSV, and the console prints how many particles this rejected.

On the reference images this excludes 13–63% of otherwise-valid detections and
moves the mean by +0.2 to +0.3 nm (~1%) — clustered particles measure slightly
*smaller*, because the watershed cut trims them. So the effect on your averages is
modest, but it is a real bias and it is now gone. Enough particles survive either
way: the sparsest image still yields 176 free-standing ones.

**Known limitation.** Two particles stacked in projection — one directly behind the
other — look like a single round particle and are not detectable by shape or contact.
They should show up as unusually dark for their size, but on these images optical
density per unit diameter has no bimodality and the darkest tenth is not consistently
larger, so a filter on it would discard ~10% of the data for no measurable gain. It
isn't implemented. If your samples aggregate heavily in the beam direction, the
overlay is the check.

## Checked against manual measurements

The 250908 S24P69 AE set has 105 diameters measured by hand in a spreadsheet, so
it is the one real check of the whole pipeline against a human.

| | n | mean | SD | CV |
|---|---|---|---|---|
| manual | 105 | 46.24 nm | 4.00 nm | 8.65% |
| temcounter | 474 | 44.29 nm | 4.02 nm | 9.08% |

**The width of the distribution matches to 0.02 nm.** The mean is 1.95 nm lower
(−4.2%, 95% CI −2.80 to −1.10). That is a clean uniform offset, not a disagreement
about the shape of the distribution: the tool places the particle boundary about
1 nm further in on the radius than the operator did, which is the soft-edge
ambiguity described under *Accuracy* below and nothing else.

`--threshold-scale 1.04` moves the pooled mean to 46.3 nm and reproduces the
manual figures almost exactly. It is not the default, because that would be
calibrating the tool to one operator on one sample. Use it when you need to match
historical numbers; leave it alone when you want internally consistent ones.

## Accuracy

Reference images, 100 particles each:

| sample | per-image means (nm) | overall |
|---|---|---|
| B1 | 29.2, 30.2, 30.1, 31.0, 30.1 | **30.1 ± 0.6 nm** |
| N42 | 53.5, 52.4, 53.1, 55.4, 54.6 | **53.8 ± 1.1 nm** |

Each set spans two magnifications with two different scale-bar labels, so the
±0.5–0.9 nm spread between images is a genuine end-to-end check of calibration and
segmentation together, not just repeatability.

Two things dominate the systematic error:

- **Where the edge is drawn.** Silica has a soft edge in TEM, so the diameter
  depends on the threshold. A ±1% shift in the Otsu threshold moves the mean by
  about ±1% (≈0.4 nm on a 30 nm particle). Tightening the shape filters from the
  defaults down to a third as many particles moved the mean by 3%, so the defaults
  are not being propped up by bad segmentations. If you have a manual measurement
  to calibrate against, `--threshold-scale` is the knob.
- **Which N get picked.** The default `--select best` takes the most circular
  particles. Their mean tracks an unbiased random sample closely (≤0.3 nm apart on
  the reference images), but the reported SD is 10-20% smaller, because rounder
  particles are also the more uniform ones:

  | image | `best` | `random` |
  |---|---|---|
  | B1 AE_1 | 30.1 ± 4.3 | 29.2 ± 4.6 |
  | B1 AE_4 | 31.0 ± 3.7 | 31.0 ± 4.0 |
  | N42 AE_1 | 53.7 ± 5.7 | 53.5 ± 6.1 |
  | N42 AE_4 | 55.2 ± 5.3 | 55.4 ± 6.7 |

  Use `best` for a clean mean, `--select random` if you are quoting a size
  distribution, a polydispersity, or a CV.

- **Selection.** Only free-standing particles are measured, so the result describes
  the dispersed population. If a sample aggregates so completely that few particles
  stand alone, the surviving ones may not represent it — the console tells you how
  many were rejected, so you can see when that is happening.
