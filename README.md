[![DOI](https://zenodo.org/badge/488654433.svg)](https://zenodo.org/badge/latestdoi/488654433)

<img src="https://github.com/N-SDM/N-SDM/blob/main/images/n-sdm_bandeau_v4.png" alt="N-SDM bandeau" align="right" width="50%"/>

# About N-SDM

Uniting species distribution modelling (SDM) techniques into one high-performance computing (HPC) pipeline, we developed N-SDM, an SDM platform aimed at delivering reproducible outputs for standard biodiversity assessments. N-SDM was built around a spatially-nested framework, intended at facilitating the combined use of species occurrence data retrieved from multiple sources and at various spatial scales. N-SDM allows combining two models fitted with species and covariate data retrieved from global to regional scales, which is useful for addressing the issue of spatial niche truncation. The set of state-of-the-art SDM features embodied in N-SDM includes a newly devised covariate selection procedure, five modelling algorithms, an algorithm-specific hyperparameter grid search and the ensemble of small-models approach. N-SDM is designed to be run on HPC environments, allowing the parallel processing of thousands of species at the same time.

# Prerequisites

Prerequisites for running N-SDM include:

1.	Linux cluster computer equipped with the Slurm workload manager

2.	availability of the modules (with versions used for N-SDM development) gcc (9.3.0), r (4.0.5), proj (5.2.0), perl (5.32.1), curl (7.76.1), geos (3.8.1) and gdal (2.4.4)

3.	a clone of the N-SDM/N-SDM GitHub repository in the working directory `git clone https://github.com/N-SDM/N-SDM.git .`

4.	an installation of the nsdm R package `remotes::install_github("N-SDM/N-SDM/scripts/nsdm-project/functions", upgrade=FALSE)`

5. [optional for running the N-SDM example] download the 6GB zip file containing example species and covariate data in the ./data folder `curl -o ./data/nsdm-project.zip  "https://unils-my.sharepoint.com/:u:/g/personal/antoine_adde_unil_ch/EQ-B2q08HQ5MuVrav33MnMQBp61DzUF9Eoi3nP_qe1FrOQ?download=1"`

6. [optional for running the N-SDM example] unzip nsdm-project.zip in the data folder `unzip ./data/nsdm-project.zip -d ./data/nsdm-project/`

# Example N-SDM run

We will run an applied example aimed at illustrating the main operations and performances of N-SDM by modelling the habitat suitability of three species (Larix decidua, Capra Ibex and Cantharellus cibarius) at 100-m resolution across Switzerland for both current (1980–2021) and future (2070–2100) periods.

## Study area

Following the spatially-nested framework (see Figure 1 and section “Highlighted features” of the N-SDM software note for details), we distinguished between “regional-” and “global-”level study areas. The regional-level area included all of Switzerland, with a total area of ≈ 40,000 km². For the global-level area, we used a bounding box covering the European continent, ranging from 32,60 °N to 71,70 °N, and from 28,56 °W to 40,21 °E, for an area of ≈ 10 Mio. km².

## Data

### Species data

Global-level species occurrence records were obtained from GBIF (https://www.gbif.org/). Regional-level records were obtained from the Swiss Species Information Center InfoSpecies (www.infospecies.ch). To limit spatial clustering effects, occurrence records will be disaggregated so that two points cannot be closer than 1 km at the global level and 200 m at the regional level. For each species and level 10,000 background absences aimed at contrasting with occurrence records will be randomly generated across the target areas.

### Covariate data

We will use a suite of 453 candidate covariates derived from 42 individual parameters belonging to 6 main categories (bioclimatic, land use and cover, edaphic, topographic, population density, transportation, and vegetation). Once the .zip data file has been unzipped, detailed information on the covariates can be found in ./data/nsdm-project/covariates/covariates.xlsx.  Only bioclimatic covariates were used for fitting the global-level model and all the others were used for the regional model (see Figure 1 in the companion Software Note for details). To account for environmental conditions within a wider area than the only coordinates of the occurrence records, covariates from the “land use and cover” category were extracted using 13 moving radii ranging from 25 m to 5 km. All covariates were standardized to zero mean and unit variance.

## N-SDM settings

N-SDM settings must be adapted to your own computing environment (e.g. paths, HPC account, partition etc.) by editing the file "settings.csv" located in ./scripts/nsdm-project/main/settings. You can also try to customize the data and/or modelling settings (e.g. covariate selection, modelling algorithm, ensembling strategy etc.). Be careful when saving "settings.csv" to use “;” as delimiters. In the same ./scripts/nsdm-project/main/settings directory you will also find the "param-grid.xslx" file that allows specifying the grid for hyperparameter tunning. An example pre-filled "expert-table.xslx" table allowing for an expert-based prefiltering of the candidate covariates is also included. Additional details on N-SDM settings, hyperparameter tunning or expert pre-filtering are provided in the N-SDM Sofwtare Note.

## Running N-SDM

Position yourself at `cd ./scripts/nsdm-project/main`, where the main N-SDM bash file (nsdm.sh) is stored. We encourage you running N-SDM in a background no hangup mode to prevent the command from being aborted automatically if logging out or exiting the shell, such as: `nohup bash nsdm.sh > nsdm.out &`. You can follow the execution of N-SDM by checking nsdm.out.

# Contributing

If you have a suggestion that would make N-SDM better, please fork the repository and create a pull request. You can also simply open an issue with the tag "enhancement".
Thanks!

# Citation

To cite N-SDM or acknowledge its use, cite the Software note as follows, substituting the version of the application that you used for ‘version 1.0':

Adde, A. et al. 2023. N-SDM: a high-performance computing pipeline for Nested Species Distribution Modelling. – Ecography 2023: e06272 (ver. 1.0). https://onlinelibrary.wiley.com/doi/10.1111/ecog.06540

# Contact

antoine.adde@unil.ch

Project Link: [https://github.com/N-SDM/N-SDM](https://https://github.com/N-SDM/N-SDM)

# Acknowledgments

N-SDM development has been conducted within the Ecospat lab https://www.unil.ch/ecospat/en/home.html.

We gratefully acknowledge financial support through the Action Plan of the Swiss Biodiversity Strategy by the Federal Office for the Environment (FOEN) for financing the Valpar.ch and SwissCatchment projects.

The Swiss Species Information Center InfoSpecies (www.infospecies.ch) supplied Swiss-level species occurrence data and expertise on species’ ecology, and we acknowledge their support regarding the database.

This research was enabled in part by the support provided by the Scientific Computing and Research Unit of Lausanne University (https://www.unil.ch/ci/dcsr).

