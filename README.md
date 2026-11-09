# Kyrgyzstan AMR Burden Estimates

This repository contains code used to generate estimates of the burden of antimicrobial resistance (AMR) in Kyrgyzstan.

## Repository structure

```text
├── code/
│   └── ...
└── README.md
```

The `code/` folder contains the scripts and analytical code used to produce the AMR burden estimates.

## About

The analyses in this repository support estimation of the health burden associated with antimicrobial resistance in Kyrgyzstan. Code and documentation may be updated as the analyses are developed and refined.

**Reproducing the analysis.** All analyses, tables and figures are reproduced by the R pipeline in `code/` (R ≥ 4.5; tested with 4.5.1). The survey data are not public, so place the three input files (`Database_with_confirmed_ID_AMR_survey_KYR_261026.xlsx`, `amr-survey_kgz_2026-08-25.xlsx` and `NewSamplingFrameDec2022.csv`) in `data/`, or set `KGZ_DATA_DIR` to the folder that holds them. Install the required packages once with `Rscript install_packages.R`, then run `Rscript run_all.R`. The pipeline runs 15 steps in order, each in a clean R process: species naming and fungal exclusion, harmonised susceptibility results (intermediate counted as resistant), MDR/XDR/PDR, design weights and calibration, prevalence and incidence, risk factors, mortality and length-of-stay models, diagnostic stewardship, and all figures as 500-dpi TIFF. A full run takes about 30–40 minutes. Tables are written to `outputs/` and `legacy/out/`, figures to `outputs/figures_TIFF_500dpi/`, and logs and `sessionInfo()` to `logs/`. The last step checks that the headline numbers are reproduced and stops with an error if they are not. Random steps (multiple imputation, bootstrap) use fixed seeds.

```bash
Rscript install_packages.R                       # once
Rscript run_all.R                                # full run
KGZ_DATA_DIR=/path/to/inputs Rscript run_all.R   # inputs stored elsewhere
Rscript run_all.R --from=09                      # restart from a step
```

## Citation

Citation information will be added when the associated analysis or publication is available.

## License
MIT License 
