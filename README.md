# Guanaco Tierra del Fuego — SPLATCHE3 + ABC pipeline

Scripts to reproduce the SPLATCHE3 spatial simulations and the ABC analyses
reported in the guanaco Tierra del Fuego manuscript
(Peña-Monroy et al., *in prep.*).

This package contains **only the SPLATCHE + ABC portion** of the pipeline.
Upstream population-genetic analyses (DAPC, PopCluster, FEEMS/FEEMSmix, LD
decay, sampling map) are archived separately.

## Repository layout

```
.
├── 01_splatche/ SPLATCHE3 landscape parameterisation and runs
│ ├── par/
│ │ ├── splatche_thin1Mb_MAFvar.par Final SPLATCHE3 marker/recombination definition
│ │ └── splatche_thin1Mb_MAFvar_sitios.txt Sampling-site definitions
│ ├── createScenario1_C.sh Build carrying-capacity time series for Scenario 1 (Refugium)
│ ├── createScenario2_C.sh Build carrying-capacity time series for Scenario 2 (Human bridge)
│ ├── runSplat_C.sh Wrapper to run SPLATCHE3 for one scenario
│ ├── run_parallel_C.sh Parallel driver over calibration replicates
│ ├── runCalib.sh Calibration wrapper (immigrant counts per m)
│ ├── runSplat_calib.sh SPLATCHE3 driver used during calibration
│ ├── K_scenario1_eruption.tif, K_scenario2_corridor.tif Rasterised K surfaces
│ ├── make_bemmels_maps.R Figure 3 — K surfaces per phase
│ ├── make_map_templates.sh Helper to prepare Figure 3 templates
│ └── fix_bemmels_v3.R Post-processing of Figure 3 panels
│
├── 02_observed/ Observed summary statistics from the empirical VCF
│ ├── popmap_dataset_C.txt Individual → sampling area map
│ ├── calc_observed_C.R Compute observed statistics (H, F_ST, F_IS…)
│ ├── calc_stats_C.R Same statistics, applied per simulation chunk
│ ├── calc_extra_stats.R Extra statistics for the 953-SNP ascertained panel
│ ├── stats_C.py Fast Python computation of raw statistics
│ └── sfs_all.py Folded SFS bootstrap for Figure 6
│
├── 03_simulated_stats/ Convert SPLATCHE .arp output into statistics
│ ├── ARP2pegas.R Import SPLATCHE .arp into pegas (Claudio S. Quilodrán, 2025)
│ └── run_calc_parallel_C.sh Parallel driver over simulation chunks
│
├── 04_abc/ ABC inference and posterior figures/tables
│ ├── final_H_isl_R1.R Main analysis (H_ISL + R1 design) — Figs 7 and 8; Tables S6, S7, S9, S13
│ ├── figS7_S8_alt_pair.R Alternative pair (H_ISL + F_ST(I–M)) — Figs S7 and S8; Tables S10–S12
│ ├── figS5.R Fig S5 — summary-statistic distributions
│ ├── figS6.R Fig S6 — sensitivity to tolerance
│ └── abc_ratios.R Table S8 — cumulative immigrants across m
│
└── data/ Empirical inputs (populated by data curator)
```

## Order of execution

1. **`01_splatche/`** — simulate under Scenarios 1 and 2 (10,000 replicates each).
   1. Run `createScenario1_C.sh` and `createScenario2_C.sh` to write the phased K files.
   2. Launch `run_parallel_C.sh` for each scenario (calls `runSplat_C.sh`).
   3. Optionally run the calibration (`runCalib.sh` + `runSplat_calib.sh`) to
      match the immigrant regime across scenarios.
2. **`02_observed/`** — compute observed summary statistics from the LD-pruned VCF.
3. **`03_simulated_stats/`** — convert each simulation's `.arp` output into a
   row of summary statistics via `ARP2pegas.R` (parallelised by
   `run_calc_parallel_C.sh`).
4. **`04_abc/`** — run the ABC inference (`final_H_isl_R1.R`) and produce
   supplementary figures and tables.

## Script → manuscript mapping

| Manuscript item                                           | Script                                    |
|-----------------------------------------------------------|-------------------------------------------|
| Figure 3 — SPLATCHE3 K surfaces                           | `01_splatche/make_bemmels_maps.R` (+ `fix_bemmels_v3.R`, `make_map_templates.sh`) |
| Figure 6 — Ne trajectory from folded SFS                  | `02_observed/sfs_all.py`                  |
| Figure 7 — Model adequacy (H_ISL + R1)                    | `04_abc/final_H_isl_R1.R`                 |
| Figure 8 — Prior/posterior (Refugium, H_ISL + R1)         | `04_abc/final_H_isl_R1.R`                 |
| Figure S5 — Summary-statistic distributions               | `04_abc/figS5.R`                          |
| Figure S6 — Sensitivity to tolerance                      | `04_abc/figS6.R`                          |
| Figures S7, S8 — Alternative pair H_ISL + F_ST(I–M)       | `04_abc/figS7_S8_alt_pair.R`              |
| Table S4 — SPLATCHE3 K per environmental class            | (documented in `01_splatche/par/`)        |
| Table S5 — Genetic diversity indices                      | `02_observed/calc_observed_C.R` (+ `calc_extra_stats.R`) |
| Tables S6, S7, S9, S13 — Posteriors and model choice (main) | `04_abc/final_H_isl_R1.R`               |
| Table S8 — Cumulative immigrants across m                 | `04_abc/abc_ratios.R`                     |
| Tables S10–S12 — Alternative pair posteriors and model choice | `04_abc/figS7_S8_alt_pair.R`          |

## Software requirements

- **SPLATCHE3** (v3.0+)
- **R ≥ 4.3** with packages: `abc`, `pegas`, `adegenet`, `hierfstat`,
  `terra`, `ggplot2`, `dplyr`
- **Python ≥ 3.9** with `numpy`, `pandas`, `scipy`, `cyvcf2`
- **PLINK v1.9** (for LD pruning of the observed panel)
- Standard Unix utilities (`bash`, `awk`, GNU parallel is optional)

## Credits

- `ARP2pegas.R` was written by **Claudio S. Quilodrán** (Univ. of Geneva,
  24.06.2025) and is redistributed here with permission. See the header
  of that file for licence and contact information.
- All other scripts by the manuscript authors.

## Known issues

See `KNOWN_ISSUES.md` for hardcoded paths that need adjustment before the
pipeline can run in a new environment.

## Licence

To be finalised before publication. A permissive licence such as MIT is
suggested for the code; empirical data will be released under CC-BY 4.0
via a separate archive.
