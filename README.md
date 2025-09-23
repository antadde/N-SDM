[![DOI](https://img.shields.io/badge/DOI-10.1111%2Fecog.06540-blue)](https://doi.org/10.1111/ecog.06540)

<img src="https://github.com/N-SDM/N-SDM/blob/main/images/n-sdm_bandeau_v4.png" alt="N-SDM bandeau" align="right" width="50%"/>

## 📢 Upcoming Release

A new version of N-SDM is coming soon with updated features and improved workflows.  
Details will be announced here.

# About N-SDM

N-SDM is a high-performance computing pipeline for Nested-Species Distribution Modeling, addressing spatial niche truncation by combining global, coarse-grain models with regional, fine-grain models. This approach preserves the full niche while keeping local resolution and precision, improving biodiversity projections under global change. The pipeline integrates automated covariate selection, multiple modeling algorithms, ensemble methods, and scalable parallel processing for efficient reproducible assessments across species and scenarios.

# Prerequisites

Prerequisites for running N-SDM include:

1. **Linux HPC cluster** equipped with the **Slurm Workload Manager**.
   
2. Availability of the following **modules** (below are the versions used during N-SDM development):
   - `gcc (12.2.0)`
   - `r (4.6.3)`
   - `gdal (3.4.4)`
   - `udunits (2.2.28)`

3. A **clone of the N-SDM GitHub repository** in the working directory:
   ```bash
   git clone https://github.com/antadde/N-SDM.git .
   ```

4. An installation of the **nsdm2 R package**:
   ```r
   remotes::install_github("antadde/N-SDM/Rpkg", upgrade=FALSE)
   ```

5. [Optional for running the N-SDM example] Download and unzip the example dataset available at https://zenodo.org/records/17177174 in the `./data` directory

# Example N-SDM run

We will run an applied example to illustrate the main operations and performance of N-SDM by modeling the habitat suitability of three species—Larix decidua, Capra ibex, and Cantharellus cibarius—at 100-meter resolution across Switzerland, for both the current period (1980–2021) and the future period (2070–2100).

## Study area

Following the spatially nested framework (see Figure 1 and the section 'Highlighted features' of the N-SDM software note for details), we distinguished between 'regional' and 'global' study areas. The regional area encompassed all of Switzerland, covering approximately 40,000 km². For the global area, we used a bounding box spanning the European continent, from 32.60°N to 71.70°N and from 28.56°W to 40.21°E, covering an area of approximately 10 million km².

## Data

### Species data

Global-level species occurrence records were obtained from GBIF (https://www.gbif.org/), while regional-level records were sourced from the Swiss Species Information Center, InfoSpecies (www.infospecies.ch). To minimize spatial clustering effects, occurrence records will be disaggregated, ensuring that no two points are closer than 1 km at the global level and 200 m at the regional level. For each species and level, 10,000 background absence points will be randomly generated across the target areas to contrast with the occurrence records.

### Covariate data

We will use a suite of 453 candidate covariates derived from 42 individual parameters, categorized into six main groups: bioclimatic, land use and cover, edaphic, topographic, population density, transportation, and vegetation. Once the .zip data file is unzipped, detailed information about the covariates can be found in ./data/nsdm-project/covariates/covariates.xlsx. Only bioclimatic covariates were used to fit the global-level model, while all other categories were used for the regional model (see Figure 1 in the companion Software Note for details). To capture environmental conditions beyond the immediate occurrence points, covariates from the 'land use and cover' category were extracted using 13 moving radii, ranging from 25 m to 5 km. All covariates were standardized to have a mean of zero and a unit variance.

## N-SDM settings

N-SDM settings must be adapted to your computing environment (e.g., paths, HPC account, partition, etc.) by editing the settings.psv file located in ./scripts/nsdm-project/main/settings with a text editor (to preserve the "|" delimiter). You can also customize the data and modeling settings (e.g., covariate selection, modeling algorithm, ensembling strategy, etc.). In the same ./scripts/nsdm-project/main/settings directory, you will find the param-grid.xlsx file, which allows you to specify the grid for hyperparameter tuning. A pre-filled example expert-table.xlsx is also included, allowing for expert-based prefiltering of candidate covariates. Additional details on N-SDM settings, hyperparameter tuning, and expert pre-filtering can be found in the N-SDM Software Note.

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
