# GM mixture-model analysis

This folder contains the GM analysis code and its prepared RAD51D input table.
The preprocessing script is included for reference; the prepared TSV file
means preprocessing is not required to run the analysis.

## Requirements

Use R or RStudio. The analysis sources `funs.R`, which loads these R packages:
`dplyr`, `tidyr`, `stringr`, `ggplot2`, `mixtools`, `rmarkdown`, `readr`, and
`rlang`. Install any packages that are not already available:

```r
install.packages(c(
  "dplyr", "tidyr", "stringr", "ggplot2", "mixtools",
  "rmarkdown", "readr", "rlang"
))
```

## Run

Set the R working directory to this `GM` folder. For example, in RStudio,
open this folder as the working directory or run:

```r
setwd("path/to/For_reviewer/Analysis/GM")
```

Then run the following in the R console:

```r
if (!dir.exists("result")) dir.create("result")
source("runGM_RAD51D_paper.R")
```

The runner sources `funs.R` and `config_RAD51D_GM_paper.R`. The configuration
uses the included `RAD51D_input_noE10_github_gm.tsv` input, analyzes the
`T2_T1` ratio, and enables the HTML report.

## Expected outputs

Results are written under `Analysis/GM/result/`:

- `MM_result.csv` — variant-level mixture-model results.
- `MM.RData` — saved R objects (`norm_ratios`, `FS_data`, and `mymm`) for
  further inspection.
- `reports/MM.html` — analysis report and visualizations.


The output directory is relative to the GM folder. Running the code again
overwrites files with the same names.
