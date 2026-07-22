# PltC-intoxication-image-analysis-script

# PltC–CD31 colocalization analysis — reproducible scripts

Deterministic image-analysis pipeline for the 2D dataset (Alexa Fluor 488-labelled
PltC intoxication, CLARITY imaging). Given the same input files, every script returns
bit-identical results.

## Files
- `pltc_analysis_common.py` — shared routines (image loading, segmentation, metrics). Imported by the two scripts below.
- `analysis_villus_restricted.py` — **main analysis**. Readout restricted to the villus interior (excludes the epithelial surface).
- `analysis_whole_tissue.py` — **robustness control**. Same pipeline over the whole DAPI+ tissue (villus restriction removed).

## Requirements
Python 3.10 with: `numpy`, `scipy`, `scikit-image`, `h5py`
```
pip install numpy scipy scikit-image h5py
```

## Input
The `2D Images` folder with four numbered condition subfolders, each containing the
raw Imaris `.ims` files (single optical planes):
```
2D Images/
  1. B6_NoToxin/            (50 .ims)
  2. B6_Intoxication/       (50 .ims)
  3. Tek-Cre_NoToxin/       (50 .ims)
  4. Tek-Cre_Intoxication/  (50 .ims)
```
Channel order in each .ims: ch0 = DAPI (405), ch1 = CD31 (640/AF647), ch2 = PltC (488).

## Run
```
python analysis_villus_restricted.py "2D Images"
python analysis_whole_tissue.py      "2D Images"
```
Each writes a per-image CSV (`villus_metrics.csv` / `whole_tissue_metrics.csv`) with
columns: condition, file, vessel_area_um2, n_puncta, puncta_density, enrichment, M1, M2,
coverage — and prints a per-condition summary plus two-sided Mann–Whitney comparisons
(B6_Intoxication vs each other condition).

## Metrics
- **puncta_density** — discrete 488-specific PltC puncta on the CD31 endothelium, per 1000 µm² vessel area (primary readout).
- **M2** — Manders' coefficient: fraction of CD31 signal overlapping PltC.
- **M1** — fraction of total PltC signal on CD31 (dominated by diffuse background here; reported for completeness).
- **enrichment** — mean on-vessel / off-vessel 488 signal.
- **coverage** — % of CD31 area covered by PltC.

Broadband autofluorescence (bright in both 488 and 640) is removed, and each image is
thresholded against its own off-vessel 488 floor (99.9th percentile).

## Note on the two scripts
The conclusion (B6_Intoxication > all other conditions across puncta density, M2 and
enrichment) holds in **both** the villus-restricted and whole-tissue analyses; the
villus restriction increases the effect size and specificity but does not change the
result.
