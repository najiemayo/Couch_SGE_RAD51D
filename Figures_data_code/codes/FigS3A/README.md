# RAD51D Supplementary Figure S3, panel A -- R source

R implementation of the **refined six-row panel**, with smaller bold text and no
leader lines. This is not Figure 4B or a molecular-structure rendering.

**Provenance:** The earlier refined TIFF was produced with Python/matplotlib.
This R translation preserves its workbook processing, category membership,
marker coordinates, colors, font sizes and layout algorithm. Fonts and rendering
engines may produce small visual differences; use the TIFF produced by this R
script when depositing it as that image's source code.

**Testing status:** R is not installed in the authoring environment and an
installation attempt could not reach the package repositories. This R script
and its R tests have **not been executed** here. The workbook-derived data were
rechecked against the previous figure's five source-data tables in Python.
See `validation/VALIDATION.md` for the exact scope. Local R execution and visual
inspection are required before journal submission.

## Quick start in RStudio

Unzip the repository, set the working directory to `RAD51D_SuppFigS3_PanelA_R`,
and place the original workbook in `data/`:

```
data/Supplementary_Data_2_updatedVarCall092826.xlsx
```

Install the packages once, then run:

```r
install.packages(c("readxl", "ragg", "systemfonts", "digest", "svglite"))

source("scripts/panel_a.R")

result <- make_panel_a(
  workbook = "data/Supplementary_Data_2_updatedVarCall092826.xlsx",
  output_dir = "output",
  reference_dir = "validation/reference"
)
```

The source workbook is not bundled or uploaded to an external repository by this
package. The script never installs dependencies automatically or fetches data.

`reference_dir` activates exact source-data checks against the preceding figure;
these CSVs are **tests only**, not input to the plot. A future workbook update
will normally fail those checks. To make a new figure from a future workbook,
run with `reference_dir = NULL`, audit the new outputs, and version its new
reference files separately. Do not replace references silently.

## Command-line execution

From the repository root:

```bash
Rscript scripts/panel_a.R \
  --workbook "data/Supplementary_Data_2_updatedVarCall092826.xlsx" \
  --outdir output \
  --reference-dir validation/reference
```

Available options: `--sheet NAME`, `--font FAMILY`, `--no-svg`, `--overwrite`.
For example, `--font Arial` requests an installed Arial family. Existing outputs
are protected unless `overwrite = TRUE` or `--overwrite` is specified.

## Data definition

The source is `Supplementary_Data_2_updatedVarCall092826.xlsx`, worksheet
`Supplementary_Data_2_updatedVar`. Header names are validated, so a reorganized
column order does not redirect a selection to the wrong field.

- Include rows whose `variant class` is `missense` (after whitespace trimming
  and case normalization). Distinct genomic/cDNA SNVs remain distinct even if
  they have the same amino-acid substitution.
- Parse `AA`, cross-check its position against `c.`, and require a consistent
  reference amino acid at each position.
- Union nonempty `Greenhouph Nature 2023` and `Rawal Nature 2023` annotations
  across all included SNVs at each position. `Greenhouph` is the workbook's
  actual column spelling and is intentionally preserved in the importer.
- Color a marker red (`#FF0000`) when **at least one** included SNV at that
  position has explicit `VarCall functional category` **P strong**. Otherwise
  color it black (`#000000`). Black does **not** mean benign.
- The category is taken from column AR's named header, not recomputed from
  posterior probabilities. No SpliceAI filter or VarCall refit is applied.
- Start/stop-loss records remain as classified by the source workbook; they do
  not acquire structural annotations merely because they are included in the
  missense-labeled input. There is no synthetic point for an unannotated residue.

The upper four rows use the Greenhough column; the lower two use the Rawal column:

| Row | Exact workbook rule | Marker |
| --- | --- | --- |
| Interaction with ATP | Any Greenhough annotation containing `ATP` | Down triangle |
| Interaction with XRCC2 | Any Greenhough annotation starting `XRCC2 (` | Up triangle |
| Interaction with RAD51C | Any Greenhough annotation containing `RAD51C (` | Square |
| Monomer stability | Any Greenhough annotation containing `Predicted monomer instability` | Circle |
| ATP binding/hydrolysis | Rawal annotation is `ATP binding`, `Mg ion coordination`, or `ATPase site` | Diamond |
| ssDNA binding | Rawal annotation is `ssDNA binding` | Down triangle |

These are the rules in the executed refined Python workflow. They intentionally
include Q115 in both upper ATP/RAD51C rows and include I311 in the upper ATP row.
They do not treat `ATP (XRCC2, ...)` alone as a direct protein-contact annotation.
Other structural contexts are exported to an omitted-context audit, not silently
forced into a displayed row. No new literature or atomic-contact analysis is
performed. The inherited `ATP binding/hydrolysis` row title is a category label,
not a new claim that RAD51D catalyzes ATP hydrolysis.

## Typography and layout

Default canvas: **12 x 6.6 inches**. TIFF: **7200 x 3960 pixels at 600 dpi**,
8-bit channels, white background, lossless LZW compression.

| Text element | Default size at the export dimensions |
| --- | ---: |
| Residue labels | 9.5 pt, bold |
| Structural-context labels | 11.5 pt, bold |
| Axis ticks | 11.5 pt, bold |
| Axis titles | 13 pt, bold |
| Legend/source headings | 10 pt, bold |
| Panel letter | 25 pt, bold |

Labels are centered on the exact residue x coordinate. The script measures text
on the active graphics device and uses compact alternating vertical tiers to
avoid overlap. No leader lines are drawn. The prior +/-0.030-inch categorical
marker separation for adjacent residues is preserved; marker y is not a score.

Font preference: Liberation Sans, then Arial, then DejaVu Sans. Fonts are not
bundled. A different font can alter collision-free label tier choices.
The script stops on detected text overlap or clipping rather than silently
exporting an invalid final layout. Reducing the final placed image width also
reduces its apparent font size; 600 dpi does not protect text from print scaling.

Adjust sizes without changing the data:

```r
style <- panel_a_style(
  residue_size_pt = 9.0,
  row_size_pt = 11.0,
  tick_size_pt = 11.0,
  title_size_pt = 12.5,
  legend_size_pt = 9.5,
  dpi = 600
)

result <- make_panel_a(
  workbook = "data/Supplementary_Data_2_updatedVarCall092826.xlsx",
  output_dir = "output_smaller_text",
  reference_dir = "validation/reference",
  style = style
)
```

Changing font sizes may require inspecting the new label layout. The main layout
has a fixed 12 x 6.6-inch coordinate system, matching the refined reference.

## Outputs

- `RAD51D_SuppFigS3_PanelA_refined_R_600dpi.tiff`
- `RAD51D_SuppFigS3_PanelA_refined_R_preview.png` (150 dpi)
- `RAD51D_SuppFigS3_PanelA_refined_R.svg` (editable text/vector objects)
- The five reconstructed source-data CSVs, including plotted and omitted rows
- `PanelA_label_layout.csv` and `PanelA_text_bounds.csv`
- `PanelA_provenance.txt`: workbook SHA-256, selection rules, counts, font sizes,
  export settings, data-check results and timestamp
- `sessionInfo.txt`: R and loaded package versions

Output files are staged temporarily and copied to the destination only after
all requested exports and checks succeed. The script never overwrites the input.
An editable SVG's appearance may depend on the font available in the editor.

## Validation tests

```bash
Rscript tests/check_panel_a.R \
  "data/Supplementary_Data_2_updatedVarCall092826.xlsx"

# Also exercise rendering and compare marker positions/shapes/colors:
Rscript tests/check_panel_a.R \
  "data/Supplementary_Data_2_updatedVarCall092826.xlsx" --render
```

Tests compare all columns of five frozen data tables. They also check current
counts, important category changes, annotation transfer across SNVs at a residue,
no posterior reclassification, and rejection of malformed/duplicate input.
The optional render test checks the 68 marker identities/coordinates/shapes/colors
against the preceding plot. It does not require identical label tiers or pixels
between R and Python. Read `validation/VALIDATION.md` before interpreting the tests.

For journal/GitHub deposition, keep the code, reviewed input-data access details,
source-data outputs and tested session information together. Select a repository
license under your institution's policy; none is assigned here on your behalf.

## Software documentation

- readxl: https://readxl.tidyverse.org/reference/read_excel.html
- ragg TIFF device: https://ragg.r-lib.org/reference/agg_tiff.html
- R grid text metrics: https://stat.ethz.ch/R-manual/R-devel/library/grid/html/grobWidth.html
- svglite: https://svglite.r-lib.org/
