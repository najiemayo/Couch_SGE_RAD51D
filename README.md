
# RAD51D SGE 

This repository provides the VarCall model tool used in the analysis of high throughput CRISPR based saturated genome editing (SGE) data, one type of multiplexed assays of variant effects (MAVEs). 

This folder contains the analysis and figure generation code used for the RAD51D
saturated-genome-editing (SGE) analyses and figure generation. The analysis to run GM and VarCall are included in their respective folders.

## Folder structure

```text
For_reviewer/
├── README.md
├── Analysis/
│   ├── GM/
│   │   ├── README.md
│   │   ├── config_RAD51D_GM_paper.R
│   │   ├── runGM_RAD51D_paper.R
│   │   ├── funs.R
│   │   └── RAD51D_input_noE10_github_gm.tsv
│   └── VarCall/
│       ├── README.md
│       ├── RAD51Dmave25.noPosNorm.ldaErrVarModel.ClinVarExclude.Rtex
│       ├── aggregated_RAD51D_E1_10_new3_SNV_091626.tsv
│       ├── clinVarExclusions.csv
│       └── preamble.tex
├── Figures_code/
│   ├── README.md
│   ├── Figure2.R
│   ├── Figure4A.R
│   ├── Figure4B.R
│   └── FigureS3A.R
└── data/
    └── Supplementary Data.xlsx
```

Each analysis folder also contains its helper scripts and any other required
inputs. Generated results are not included; they are created when the code is
run.

## Running the analyses

- **GM model analysis:** Follow [Analysis/GM/README.md](Analysis/GM/README.md).
  Run it from the `Analysis/GM` directory.
- **VarCall:** Follow [Analysis/VarCall/README.md](Analysis/VarCall/README.md).
  Run it from the `Analysis/VarCall` directory. This Bayesian analysis requires
  JAGS and may take about 12 hours on the long setting which is used in the publication. The current setting is for a short run and finishes in around 20 minutes.
- **Figures:** Follow [Figures_code/README.md](Figures_code/README.md). It
  explains the R packages, workbook selection, commands, and generated figure
  and data files. The workbook is provided in `data/`.

The supplied GM and VarCall input data are generated from Supplementary Data.xlsx, from running the preprocessing_GM.R or preprocessing_VarCall.R . Outputs are written to the locations described in each analysis
README.

