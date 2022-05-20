<img src="https://github.com/AnAdde/N-SDM/blob/main/images/n-sdm_bandeau_v4.png" alt="N-SDM bandeau" width="100%"/>

# About N-SDM

Bringing together leading-edge species distribution modelling (SDM) techniques into one high-performance computing (HPC) pipeline, we developed N-SDM, an SDM platform mainly written in R language aimed at delivering reproducible outputs for standard biodiversity assessments. A key characteristic of N-SDM is that it has been built around a spatially-nested framework for facilitating the combined use of species occurrence data retrieved from multiple sources and at different spatial scales. In practice, this is done in N-SDM by combining two models fitted with scale-specific species and covariate data (e.g.: global and local). Among other benefits, this strategy allows addressing the potential issue of niche truncation. A second key characteristic of N-SDM is that it has been explicitly designed for HPC environments, filling the gap in the availability of SDM platforms directly connectable to computing clusters. Among several others, the set of state-of-the-art SDM features embodied in N-SDM includes a newly devised covariate selection procedure, five modelling algorithms, an algorithm-specific hyperparameter grid search and the ensemble of small models approach. Finally, N-SDM has been designed to be easily customizable and allows to directly benefit from the inputs of both the computational modelling community and species experts.

# Prerequisites

Prerequisites for running N-SDM include:

1.	you are working on a Linux HPC cluster equipped with the Slurm workload manager

2.	the list of available modules in your system include gcc, r, proj, perl, curl, geos and gdal

3.	you have cloned the AnAdde/N-SDM repository in your working directory `git clone https://github.com/AnAdde/N-SDM.git .`

4.	you have installed the nsdm R package `remotes::install_github("AnAdde/N-SDM/scripts/nsdm-project/functions", upgrade=FALSE)`

5. [optional for running the N-SDM example] you have downloaded the 6GB zip file containing example species and covariate data in the ./data folder `curl -o ./data/nsdm-project.zip  https://drive.switch.ch/index.php/s/BlN3GV7x8M8CCI7/download`

6. [optional for running the N-SDM example] you have unzipped nsdm-project.zip in the data folder `unzip ./data/nsdm-project.zip -d ./data/nsdm-project/`

# Example N-SDM run

<img align="right" alt="N-SDM logo" src="https://github.com/AnAdde/N-SDM/blob/main/images/n-sdm_logo_v2.png" width="30%"/>

We will run an applied example aimed at illustrating the main operations and performances of N-SDM by modelling the habitat suitability of three species (Larix decidua, Capra Ibex and Cantharellus cibarius) at 100-m resolution across Switzerland.

## Study area

Following the spatially-nested framework described earlier (see Figure 1 and section “Highlighted features” in the N-SDM Software Note), we distinguished between “local” and “global” -level study areas. The local-level one included all of Switzerland, with a total area of ≈ 40,000 km². For the global-level one, we used a bounding box covering the European continent, with coordinates 32,60N, 71,70N, 28,56W and 40,21E.

## Data

### Species data

Global-level species occurrence records were obtained from GBIF (https://www.gbif.org/). Local-level records were obtained from the Swiss Species Information Center InfoSpecies (www.infospecies.ch). To limit spatial clustering effects, occurrence records will be disaggregated so that two points cannot be closer than 1 km at the global level and 200 m at the local level. For each species and level 10,000 background absences aimed at contrasting with occurrence records will be randomly generated across the target areas.

### Covariate data

We will use a suite of 453 candidate covariates derived from 42 individual parameters and belonging to 6 main categories (bioclimatic, land use and cover, edaphic, topographic, population density, transportation, and vegetation). Once the .zip data file has been unzipped, detailed information on the covariates can be found in ./data/nsdm-project/covariates/covariates.xlsx.  Only bioclimatic covariates were used for fitting the global-level model and all the others were used for the local model (see Figure 1 in the companion Software Note for details). To account for environmental conditions within a wider area than the only coordinates of the occurrence records, covariates from the “land use and cover” category were extracted using 13 moving radii ranging from 25 m to 5 km. All covariates were standardized to zero mean and unit variance.

## N-SDM settings

N-SDM settings must be adapted to your own computing environment (e.g. paths, HPC account, partition etc.) by editing the file "settings.csv" located in ./scripts/nsdm-project/main/settings. You can also try to customize the data and/or modelling settings (e.g. covariate selection, modelling algorithm, ensembling strategy etc.). Be careful when saving "settings.csv" to use “;” as delimiters. In the same ./scripts/nsdm-project/main/settings directory you will also find an the "param-grid.xslx" file that allows specifying the grid for hyperparameter tunning. An example pre-filled "expert-table.xslx" table allowing for expert-based prefiltering of taxon-specific candidate covariates is also included. Additional details on N-SDM settings, hyperparameter tunning or expert pre-filtering are provided in the companion N-SDM Sofwtare Note.

## Running N-SDM

Position yourself at `cd ./scripts/nsdm-project/main`, where the main N-SDM bash file (nsdm.sh) is stored. We encourage you running N-SDM in a background no hangup mode to prevent the command from being aborted automatically if logging out or exiting the shell, such as: `nohup bash nsdm.sh > nsdm.out &`. You can follow the execution of N-SDM by checking nsdm.out.

# Contributing

If you have a suggestion that would make N-SDM better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Thanks!

# Citation

To cite N-SDM  or acknowledge its use, cite the companion Software note as follows, substituting the version of the application that you used for ‘version 1.0':

Adde, A. et al. 2022. N-SDM: a high-performance computing  pipeline for fitting Nested Species Distribution Models. – Ecography 2022: XXX (ver. 1.0).

# Contact

Antoine Adde - aadde@unil.ch

Project Link: [https://https://github.com/AnAdde/N-SDM](https://https://github.com/AnAdde/N-SDM)

# Acknowledgments

N-SDM development has been conducted within the ECOSPAT lab https://www.unil.ch/ecospat/en/home.html.

Financial support through the Action Plan of the Swiss Biodiversity Strategy via the Federal Office for the Environment and the www.Valpar.ch project is gratefully acknowledged.

The Swiss Species Information Center InfoSpecies (www.infospecies.ch) supplied Swiss-level species occurrence data and we acknowledge their support regarding the database.

This research was enabled in part by the support provided by the Scientific Computing and Research Unit of Lausanne University (https://www.unil.ch/ci/dcsr).

