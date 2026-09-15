# Known issues to address before publication

## Hardcoded absolute paths

The following scripts contain absolute paths pointing to the original
working directory (`/home/apena/Splatche/…` or `~/Splatche/datasetC`).
These need to be replaced with a portable path convention (for example a
`REPO_ROOT` environment variable, or `here::here()` in R) before the
pipeline can be run in a fresh environment.

| File                                          | Line(s)   | Pattern to replace                                              |
|-----------------------------------------------|-----------|-----------------------------------------------------------------|
| `02_observed/calc_extra_stats.R`              | 8         | `source("/home/apena/Splatche/ARP2pegas.R")`                    |
| `02_observed/calc_extra_stats.R`              | 19        | `/home/apena/Splatche/datasetC/outscenario…_C_PANEL953`         |
| `02_observed/calc_extra_stats.R`              | 25        | `/home/apena/Splatche/datasetC/param.txt`                       |
| `02_observed/calc_extra_stats.R`              | 111       | `/home/apena/Splatche/datasetC/extra_stats_scen…`               |
| `02_observed/calc_observed_C.R`               | 10        | `/home/apena/Splatche/datasetC/datasetC_bi_plink.ped`           |
| `02_observed/calc_observed_C.R`               | 54        | `/home/apena/Splatche/datasetC/observed_C.obs`                  |
| `02_observed/calc_observed_C.R`               | 61        | `/home/apena/Splatche/datasetC/basic_stats_observed_C.rds`      |
| `02_observed/calc_stats_C.R`                  | 8         | `source("/home/apena/Splatche/ARP2pegas.R")`                    |
| `02_observed/calc_stats_C.R`                  | 15        | `/home/apena/Splatche/datasetC/outscenario…_C`                  |
| `02_observed/calc_stats_C.R`                  | 16        | `/home/apena/Splatche/datasetC/param.txt`                       |
| `02_observed/calc_stats_C.R`                  | 43        | `/home/apena/Splatche/datasetC/chunks_scen…`                    |
| `03_simulated_stats/run_calc_parallel_C.sh`   | 6         | `cd /home/apena/Splatche/datasetC`                              |
| `04_abc/abc_ratios.R`                         | 2         | `setwd(path.expand("~/Splatche/datasetC"))`                     |
| `04_abc/figS7_S8_alt_pair.R`                  | 8         | `setwd(path.expand("~/Splatche/datasetC"))`                     |
| `04_abc/final_H_isl_R1.R`                     | 15        | `setwd(path.expand("~/Splatche/datasetC"))`                     |
| `01_splatche/runSplat_C.sh`                   | 5         | `outdir="/home/apena/Splatche/datasetC/outscenario${scen}_C"`   |
| `01_splatche/make_bemmels_maps.R`             | 9         | `base_dir <- "/home/apena/Splatche/datasetC"`                   |
| `01_splatche/fix_bemmels_v3.R`                | 9         | `base_dir <- "/home/apena/Splatche/datasetC"`                   |
| `01_splatche/runCalib.sh`                     | 3         | `out="/home/apena/Splatche/datasetC/calib_${scen}.txt"`         |
| `01_splatche/runSplat_calib.sh`               | 3         | `outdir="/home/apena/Splatche/datasetC/calib_scen${scen}"`      |
| `01_splatche/make_map_templates.sh`           | 2         | `cd /home/apena/Splatche/datasetC`                              |

Suggested pattern for R scripts:

```r
REPO_ROOT <- Sys.getenv("REPO_ROOT", unset = normalizePath("."))
```

Suggested pattern for shell scripts:

```bash
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
```

## Empirical inputs to add under `data/`

The scripts expect the following files, which are not shipped in this
package and must be provided by the data curator:

- `datasetC_biallelic.recode.vcf` — LD-pruned biallelic SNP panel
- `datasetC_bi_plink.ped`, `datasetC_bi_plink.map` — PLINK format of the
  same panel
- `param.txt` — SPLATCHE3 parameter table (m, K_REFUGE, K_CORRIDOR values
  per replicate)

## Intermediate output directories

The pipeline writes to directories named `outscenario1_C/`,
`outscenario2_C/`, `chunks_scen1/`, `chunks_scen2/`, `calib_scen1/`,
`calib_scen2/`. These should be created under a `results/` folder and
kept out of version control (see `.gitignore`).


