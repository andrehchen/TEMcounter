# TEMcounter

Measures nanoparticle diameters in nm from TEM `.tif` micrographs, calibrating the
pixel size from the scale bar burnt into each image. Single script, `temcounter.py`.

## Running it

Use `.venv/bin/python`, never system python — the Homebrew install is externally
managed and rejects `pip install`. Recreate the venv with
`python3 -m venv .venv && .venv/bin/pip install -r requirements.txt`.

Needs Tesseract on the system (`brew install tesseract`) to read the scale-bar
label. The tool refuses to run without it rather than falling back silently,
because a wrong scale produces believable but incorrect diameters.

```sh
.venv/bin/python temcounter.py OneDrive_1_2026-07-23 -n 100
```

## Data

Raw images are gitignored and live only on this machine / OneDrive.

- `OneDrive_1_2026-07-23/` — Hitachi HT7700, 16-bit, tick-mark scale bar, full
  TIFF metadata. Samples B1 (~30 nm) and N42 (~54 nm), 5 images each.
- `Some previous pictures to test out/` — AMT XR280, 8-bit, **solid** scale bar,
  no Hitachi tags, dark separator strip below the info bar. A different format
  that must keep working; it broke three separate assumptions when first added.
- `250908 S24P69 AE size.xlsx` — 105 hand measurements, ground truth for the
  `250908 S24 AE_*` images. **Check against this after touching segmentation.**
  Expected: SD matches to ~0.02 nm, mean runs ~4% low.

## Working habits

- Iterate on one image. A full 10-image run takes ~2 minutes and most changes
  don't need it. Use the whole set only to confirm a finished change.
- **Look at the overlay PNGs.** Every segmentation bug found so far — watershed
  slivers, merged pairs, 5 nm stain specks — was caught by looking, not by a
  metric. `--hist-bins 0` silences the terminal histogram when not needed.
- Never commit `.tif` or `.xlsx`. Raw data stays out of git.

## Decisions that look arbitrary but aren't

- **The scale-bar label is the full tick-to-tick span**, not one division.
  Confirmed by four magnifications with three different labels agreeing to
  <0.5% once scaled by magnification. Don't "fix" this.
- **Roundness** (inscribed-circle diameter ÷ equivalent diameter) exists because
  it is ~0.7 for *both* a merged pair and a mis-split fragment, while circularity
  and solidity pass both. It is the filter that does the real work.
- **`--max-cluster` defaults to 1**, so only free-standing particles are measured.
  Clustered particles read 0–3.6% smaller; excluding them removes that bias at the
  cost of 13–63% of detections. Plenty survive.
- **`--select` defaults to `best`** at the user's request. It tracks an unbiased
  sample on the mean but narrows the SD by 10–20%. If a distribution width, CV or
  polydispersity is being reported, use `--select random`.
- **Automatic minimum diameter** at a third of the typical particle size. Stain
  specks are tiny but perfectly round, so no shape filter catches them.
- **The metadata fallback** (`nm/px = 10949 / magnification`) is fitted to this
  one HT7700 + XR81-B pair. It is not a general calibration and warns loudly.
- **The ~4% offset vs hand measurement** is where the soft silica edge gets
  thresholded, not a bug. `--threshold-scale 1.04` reproduces the manual numbers;
  it is deliberately not the default.

## Git

Public repo `andrehchen/TEMcounter`, branch `main`. Only the four source files are
tracked. `TEMcounter-for-colleague.zip` is a rebuilt artifact for Windows users
(code + `SETUP-WINDOWS.txt`, no venv, no images).
