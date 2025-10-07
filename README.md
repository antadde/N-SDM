[![DOI](https://img.shields.io/badge/DOI-10.1111%2Fecog.06540-blue)](https://doi.org/10.1111/ecog.06540)

<img src="https://github.com/N-SDM/N-SDM/blob/main/images/n-sdm_bandeau_v4.png" alt="N-SDM bandeau" align="right" width="40%"/>

<h3>📢 Upcoming Release – September/October 2025</h3>

A new version of N-SDM is coming soon with updated features and improved workflows. Details will be announced here.

# About N-SDM

N-SDM is a high-performance computing pipeline for Nested-Species Distribution Modeling, addressing spatial niche truncation by combining global, coarse-grain models with regional, fine-grain models. This approach preserves the full niche while keeping local resolution and precision, improving biodiversity projections under global change. The pipeline integrates automated covariate selection, multiple modeling algorithms, ensemble methods, and scalable parallel processing for efficient reproducible assessments across species and scenarios.

# Prerequisites

Prerequisites for running N-SDM include:

1. **Linux HPC cluster** equipped with the **Slurm Workload Manager**.
   
2. A **clone of the N-SDM GitHub repository** in the working directory:
   ```bash
   git clone https://github.com/antadde/N-SDM.git nsdm
   ```

3. An installation of the **nsdm2 R package**:
   ```r
   # --- Stable release --- NOT OUT YET !!!
   remotes::install_github("antadde/N-SDM/Rpkg", ref = "v2.0.0")

   # --- Development version ---
   remotes::install_github("antadde/N-SDM/Rpkg")
   ```
   
4. [Optional for running the N-SDM example] Download and unzip the example dataset available at https://zenodo.org/records/17177174 in the `./data` directory. Follow the instructions provided on Zenodo.

## Tested cluster environments

N-SDM v2.0.0 (and nsdm2 R package) has been successfully installed and executed on multiple high-performance computing (HPC) systems. The following configurations describe the module environments and R setups that were verified to work.

#### ETH Zürich – Euler Cluster
***Cluster environment:***
- Operating system: Ubuntu 22.04.5 LTS
- Workload manager: SLURM
- Modules:
  - `stack/2024-06`
  - `gcc/12.2.0`
  - `gdal`
  - `udunits/2.2.28`
  - `r/4.3.2`
- Required modules are specified in the `settings.psv` file through the parameters:
  - `module_r = r/4.3.2`
  - `module_others = stack/2024-06,gcc/12.2.0,gdal,udunits/2.2.28`

***R environment:***
- Platform: R version 4.3.2 (2023-10-31)
- Verified with the following [sessionInfo()](./documentation/R_session_info/ETHZ_Euler.txt)

#### University of Lausanne – Curnagl Cluster
***Cluster environment:***
- Operating system: Linux (Rocky Linux 8)
- Workload manager: SLURM
- Modules:
  - `r-light` (which uses a container)
- Same settings parameters used:
  - `module_r = r-light`
  - `module_others = ""` (no additional modules required)

***R environment:***
- Platform: R 4.3.2
- Verified with the following [sessionInfo()](./documentation/R_session_info/UNIL_Curnagl.txt)

# Example N-SDM run

We will run an applied example to demonstrate the main operations and outputs of N-SDM by modeling the habitat suitability of three species (Larix decidua, Capra ibex, and Cantharellus cibarius) at 100-meter resolution across Switzerland, for the current period (1980–2021) and a future period (2070–2100).

## Study area

Following the spatially nested framework (see Figure 1 and the section 'Highlighted features' of the N-SDM software note for details), we distinguished between 'regional' and 'global' study areas. The regional area encompassed all of Switzerland, covering approximately 40,000 km². For the global area, we used a bounding box spanning the European continent, from 32.60°N to 71.70°N and from 28.56°W to 40.21°E, covering an area of approximately 10 million km².

## Data

Further details on covariate and species data preparation are available in the [DATA_PREPARATION.odt](./documentation/DATA_PREPARATION_20251002.odt) document located in the `documentation` directory.

### Species data

Global-level species occurrence records were obtained from GBIF, and regional-level records from the Swiss Species Information Center. To reduce spatial clustering, records were disaggregated so that no two points were closer than 1 km at the global level and 200 m at the regional level. For each species and level, 10,000 background absence points were randomly generated across the target areas to contrast with the occurrence records.

### Covariate data

We used suite of 374 candidate covariates categorized into six main categories: bioclimatic, land use and cover, edaphic, topographic, population density, transportation, and vegetation. Only bioclimatic covariates were used to fit the global-level model, while all other categories were used for the regional model. To capture environmental conditions beyond the immediate occurrence points, covariates from the 'land use and cover' category were extracted using 13 moving radii, ranging from 25 m to 5 km. All covariates were standardized to have a mean of zero and a unit variance.

## N-SDM settings

Further details on N-SDM settings and hyperparameter tuning options are available in the [SETTINGS_DETAILS.odt](./documentation/SETTINGS_DETAILS_20251002.odt) and [ALGORITHMS_PARAMETERS.odt](./documentation/ALGORITHMS_PARAMETERS_20251003.odt) documents located in the `documentation` directory.

N-SDM settings must be adapted to your computing environment (e.g., paths and module versions) by editing the `settings.psv` file in `./scripts/settings` with an editor that preserves the `"|"` delimiter. You can also customize data and modeling options such as covariate selection, modeling algorithms, and ensembling strategies.  

In the same `./scripts/settings` directory:  
- `param_grid.psv` defines the grid for hyperparameter tuning  
- `expert_table.psv` (pre-filled) allows expert-based prefiltering of candidate covariates  

## Running N-SDM

Position yourself in the `./scripts` directory, where the main N-SDM bash file (nsdm.sh) is stored.

It is recommended to run N-SDM in the background using no hangup mode, so the process continues even if you log out or close the shell. Use the following command:

```bash
nohup bash nsdm.sh > nsdm.out &
```

You can monitor the execution of N-SDM by checking the `nsdm.out` file.

# Contributing

If you have suggestions to improve N-SDM, please fork the repository and submit a pull request.

# Citation

To cite N-SDM or acknowledge its use, cite the Software note as follows, substituting the version of the application that you used for ‘version 1.0':

Adde, A. et al. 2023. N-SDM: a high-performance computing pipeline for Nested Species Distribution Modelling. – Ecography 2023: e06272 (ver. 1.0). https://onlinelibrary.wiley.com/doi/10.1111/ecog.06540

# Contact

antoine.adde@eawag.ch

# Acknowledgments

N-SDM development was conducted within the [Guisan](https://www.unil.ch/ecospat/en/home.html) and [Altermatt](https://www.altermattlab.ch/) labs. This work was supported by the Action Plan of the Swiss Biodiversity Strategy, with funding from the Swiss Federal Office for the Environment to the [ValPar.CH project](https://www.valpar.ch/index_en.php?page=home_en). It was later integrated into the [SPEED2ZERO project](https://speed2zero.ethz.ch/en/), supported by the ETH-Board through the Joint Initiatives scheme.
