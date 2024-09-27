[![DOI](https://zenodo.org/badge/488654433.svg)](https://zenodo.org/badge/latestdoi/488654433)

<img src="https://github.com/N-SDM/N-SDM/blob/main/images/n-sdm_bandeau_v4.png" alt="N-SDM bandeau" align="right" width="50%"/>

# About N-SDM

Uniting species distribution modeling (SDM) techniques into a high-performance computing (HPC) pipeline, we developed N-SDM, a platform designed to deliver reproducible outputs for standard biodiversity assessments. N-SDM is built around a spatially nested framework, aimed at facilitating the combined use of species occurrence data from multiple sources and across various spatial scales. N-SDM enables the combination of two models, fitted with species and covariate data retrieved from global to regional scales, which helps address the issue of spatial niche truncation. The set of state-of-the-art SDM features in N-SDM includes a newly devised covariate selection procedure, five modeling algorithms, an algorithm-specific hyperparameter grid search, and the ensemble of small-models approach. N-SDM is designed for HPC environments, enabling the parallel processing of thousands of species simultaneously.

# Prerequisites

Prerequisites for running N-SDM include:

1. **Linux cluster computer** equipped with the **Slurm workload manager**.
   
2. Availability of the following **modules** (with the versions used during N-SDM development):
   - `gcc (9.3.0)`
   - `r (4.0.5)`
   - `proj (5.2.0)`
   - `perl (5.32.1)`
   - `curl (7.76.1)`
   - `geos (3.8.1)`
   - `gdal (2.4.4)`

3. A **clone of the N-SDM GitHub repository** in the working directory:
   ```bash
   git clone https://github.com/N-SDM/N-SDM.git .
   ```

4. An installation of the **nsdm R package**:
   ```r
   remotes::install_github("N-SDM/N-SDM/scripts/nsdm-project/functions", upgrade=FALSE)
   ```

5. [Optional for running the N-SDM example] Download the 6GB zip file containing example species and covariate data into the `./data` folder:
   ```bash
   wget -O ./data/nsdm-project.zip https://unils-my.sharepoint.com/:u:/g/personal/antoine_adde_unil_ch/EQ-B2q08HQ5MuVrav33MnMQBp61DzUF9Eoi3nP_qe1FrOQ?download=1
   ```

6. [Optional for running the N-SDM example] Unzip `nsdm-project.zip` in the `data` folder:
   ```bash
   unzip ./data/nsdm-project.zip -d ./data/nsdm-project/
   ```

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

Position yourself at `cd ./scripts/nsdm-project/main`, where the main N-SDM bash file (`nsdm.sh`) is stored. 

We encourage running N-SDM in a background "no hangup" mode to prevent the command from being aborted if you log out or exit the shell. You can do this with the following command:

```bash
nohup bash nsdm.sh > nsdm.out &
```

You can monitor the execution of N-SDM by checking the `nsdm.out` file.

# Contributing

If you have a suggestion that would make N-SDM better, please fork the repository and create a pull request. You can also simply open an issue with the tag "enhancement".
Thanks!

# Citation

To cite N-SDM or acknowledge its use, cite the Software note as follows, substituting the version of the application that you used for ‘version 1.0':

Adde, A. et al. 2023. N-SDM: a high-performance computing pipeline for Nested Species Distribution Modelling. – Ecography 2023: e06272 (ver. 1.0). https://onlinelibrary.wiley.com/doi/10.1111/ecog.06540

# Contact

antoine.adde@eawag.ch

Project Link: [https://github.com/N-SDM/N-SDM](https://https://github.com/N-SDM/N-SDM)

# Acknowledgments

N-SDM development was conducted within the [Ecospat lab](https://www.unil.ch/ecospat/en/home.html).

We gratefully acknowledge financial support from the Federal Office for the Environment (FOEN) through the Action Plan of the Swiss Biodiversity Strategy, which funded the Valpar.ch and SwissCatchment projects.

The Swiss Species Information Center, [InfoSpecies](https://www.infospecies.ch), provided Swiss-level species occurrence data and expertise on species ecology, and we appreciate their support regarding the database.

This research was also made possible, in part, by the support provided by the [Scientific Computing and Research Unit](https://www.unil.ch/ci/dcsr) of the University of Lausanne.

