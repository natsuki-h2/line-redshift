# line-redshift

Find **frequency-consistent redshift hypotheses** for astronomical spectral-line pairs, draw `lines.freq` redshift ladders, and announce the results in the terminal or with an optional **local desktop notification**.

This is a reusable successor to the `all_spectra_plot/plot_spectra_ALMAcy7.py` and `plot_natsuki.py` workflow. Instead of editing a hard-coded anchor frequency and visually reading off every possible identification, supply a pair catalogue and receive all compatible assignments. Single-anchor plotting is also retained.

**The program proposes hypotheses, not secure redshifts.** It does not detect lines, infer their reliability, fit centroids, rank by astrophysical probability, or independently establish that paired signals belong to the same object. No astronomical observations or unpublished candidate catalogue are included; the example is synthetic.

## Install and run

Python 3.9 or newer:

```bash
git clone https://github.com/natsuki-h2/line-redshift.git
cd line-redshift
python -m venv .venv
source .venv/bin/activate
python -m pip install .
line-redshift --lines lines.freq --pairs examples/pairs.csv --out results/example
```

The legacy-style entry point is equivalent:

```bash
python all_spectra_plot.py --lines lines.freq --pairs examples/pairs.csv \
  --tolerance-kms 100 --z-min 0 --z-max 10 --out results/pairs --notify
```

`--notify` requests a local macOS/Linux desktop notification after the run. Availability depends on the OS/session; failure falls back to the terminal. There is no email, messaging-service integration, background scheduler, or network access in the analysis code.

## Inputs

### Rest frequencies

Legacy two-column whitespace format, frequency then a whitespace-free transition name:

```text
# Rest frequencies in MHz
230538.0000 CO21
345795.9899 CO32
```

`lines.freq` uses **MHz**, while observed frequencies use **GHz**. Other rest-frequency units require `--line-unit Hz`, `kHz`, or `GHz`. Comments start with `#`. Duplicate names and non-positive/non-finite frequencies are rejected.

The bundled 30-entry table is retained from the original personal plotting workflow. It includes CO, CI, H2O and fine-structure transitions, but **is not a complete or newly validated spectroscopic catalogue**. Several values are rounded, and original per-entry bibliographic references were not present. Verify rest frequencies and species/transition labels against a suitable laboratory catalogue before publication. Supply your own expanded/corrected table as needed; output JSON records its SHA-256 hash. Near-identical rest-frequency aliases are not counted as distinct transitions.

### Observed pairs

```csv
pair_id,frequency_a_ghz,frequency_b_ghz,frame_a,frame_b
example,76.846,115.26532996666667,LSRK,LSRK
```

The three ID/frequency columns are required. IDs must be unique. Use measured centroids where possible; a blind-search window centre is only a preliminary estimate. Pair selection, positional association, duplicate detection removal, and observational independence are upstream responsibilities. A pair with a negative or absorption signal may be investigated too; this tool does not assign significance from its sign or S/N.

Frame columns should state the actual frequency reference frame. A mixed pair such as LSRK/BARYCENT is **rejected by default**. `--allow-mixed-frames` is an explicit exploratory override; it flags results as unverified and does **not** convert frames. Missing frames are also flagged. Proper transformations require observing metadata and are outside this tool.

### Optional spectra and coverage

`--spectra-manifest spectra.csv` associates each pair with one or more spectrum segments:

```csv
pair_id,path,flux_unit
example,spectra/band3.dat,mJy/beam
example,spectra/band4.dat,mJy/beam
```

Each spectrum is a two-column, whitespace-delimited numeric file: observed frequency **GHz**, flux density. A `#` header/comment is allowed. Paths are relative to the manifest. Spectra for a pair must already share a flux unit; the tool neither rescales fluxes nor matches beams. Segments are plotted separately, avoiding lines connecting separate files. Split any within-file spectral gaps into separate files. In anchor mode use `pair_id=anchor`.

Specify actual usable coverage, including gaps, with repeated options:

```bash
line-redshift --lines lines.freq --pairs examples/pairs.csv \
  --coverage 75:78 --coverage 113:117 --out results/coverage
```

Coverage flags in `hypotheses.json` only indicate whether a predicted frequency lies within a supplied interval. **Coverage does not imply sufficient sensitivity; an undetected predicted line is not automatically an exclusion.** No coverage is inferred from the extent of plotted spectra.

## Matching rule

For each assignment of **two distinct transitions** to the two observed lines:

```text
z_a = rest_a / observed_a - 1
z_b = rest_b / observed_b - 1
z   = (z_a + z_b) / 2
velocity_mismatch = c * abs(z_a - z_b) / (1 + z)
```

Keep the assignment if `z_min <= z <= z_max` and `velocity_mismatch <= tolerance_kms`.

The tolerance is the **total pair-wise velocity mismatch**, not a +/- tolerance for each line. The default 100 km/s is an illustrative starting value, not a calibrated error or a line-width-based probability. Choose it from centroid errors and plausible relative kinematics; document sensitivity to that choice. The reported `z` is the arithmetic mean of the two implied redshifts, not a fitted redshift with an uncertainty.

All passing assignments are retained, sorted by mismatch for convenience. That ordering is not a posterior ranking. One line is assigned to at most one transition in an assignment, and the same transition cannot explain both inputs. Close but non-overlapping observed features may still be components of a single transition; the caller must check that interpretation. This version evaluates pairs independently; it does not jointly optimize an entire multi-line source.

## Outputs

- `hypotheses.csv`: every passing assignment, implied redshifts, mismatch, predicted frequencies, frame status.
- `hypotheses.json`: the same results, no-match inputs, parameters, input/catalogue hashes and coverage predictions.
- `report.md`: readable per-pair results and links to figures.
- `pair_001.png` / `.pdf`, etc.: spectra/frequency markers plus redshift-ladder hypotheses. The report maps numeric figure names to input IDs.
- Optional `pair_001_spectra.png` / `.pdf`: zoomed individual spectrum segments, with at most three hypothesis overlays.
- Terminal summary for every input, including zero matches. Optional local completion notification.

`--plot-limit 12` limits only plotted rows; CSV/JSON/report retain every hypothesis. `--no-plots` produces reports without figures. Use a new output directory for each run; matching filenames are overwritten, and stale files from an older longer run are not automatically deleted.

No match means no solution within the supplied catalogue, redshift range and tolerance. It does **not** prove a source is spurious. A match does not prove a redshift either: line flux ratios, other covered transitions and limits, spectral morphology, calibration, confusion and independent observations require scientific follow-up.

## Original single-anchor workflow

```bash
line-redshift --lines lines.freq --anchor-ghz 105.01 \
  --z-min 0 --z-max 10 --coverage 95:110 --coverage 137:160 \
  --out results/anchor
```

Every allowed transition is tried as the anchor identification and the complete line list is shifted accordingly. No hand-coded highlighted redshift interval or arbitrary slice of transitions is used.

## Tests and provenance

```bash
python -m unittest discover -s tests -v
```

Tests cover synthetic redshift recovery, reversed pair ordering, exact velocity-tolerance boundaries, ambiguous solutions, unit parsing, invalid input, same-transition alias rejection, mixed/unknown frames, and CLI outputs. GitHub Actions runs tests on Python 3.9 and 3.12 and generates the synthetic example.

The original plotting approach and `lines.freq` came from Natsuki Hayatsu's existing `all_spectra_plot` workflow. The implementation here is a portable rewrite: personal paths, observing products, hard-coded target frequencies/redshift highlights and automatic shell `open` calls have been removed. The table is intentionally preserved rather than silently replaced.

Code is MIT licensed. Rest-frequency values are scientific data; their inclusion is not a claim of new laboratory measurements or authoritative catalogue validation.

## 日本語メモ

`lines.freq` を使う従来の図示に、輝線ペアからの候補探索・一覧出力・実行時通知を加えたものです。公開版の例は合成データです。今回の実天体カタログやスペクトルは含めていません。

`--tolerance-kms` は二本の輝線が示す赤方偏移のずれを速度換算した許容値です。「100 km/s以内だから赤方偏移確定」という判定はしません。周波数基準系の混在、線幅と中心周波数の不確かさ、同一遷移の速度成分、別の遷移の検出限界は別途検証してください。
