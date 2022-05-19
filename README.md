<img src="https://github.com/AnAdde/N-SDM/blob/main/images/n-sdm_bandeau_v4.png" alt="N-SDM bandeau" width="100%"/>

# About N-SDM

Bringing together leading-edge species distribution modelling (SDM) techniques into one high-performance computing (HPC) pipeline, we developed N-SDM, an innovative SDM platform mainly written in R language aimed at delivering scalable and reproducible outputs for standard biodiversity assessments. A key characteristic of N-SDM is that it has been built around a spatially-nested framework aimed at facilitating the combined use of species occurrence data retrieved from multiple sources and at different spatial scales. In practice, this was done in N-SDM by combining two models fitted with scale-specific (e.g.: global and local) species and covariate data. Among other benefits, this strategy allowed addressing the potential issue of niche truncation. A second characteristic of N-SDM is that it has been explictly designed for HPC environments, filling the gap in the availability of SDM platforms directly connectable to computing clusters that are becoming increasingly necessary to handle big-data contexts, the new normal in biodiversity assessments. Among several others, the set of features included in N-SDM encompasses a newly devised covariate selection procedure, an algorithm-specific hyperparameter grid search, five optional modelling algorithms and the ensemble of small models approach. Finally, N-SDM has been designed to be easily customizable and allows directly benefiting of the inputs from both the computational modelling community and species experts.


# Prerequisites

Before running N-SDM, you should make sure you have an appropriate set up with the following requirements:

1.	you are working on a Linux HPC cluster equipped with the Slurm workload manager

2.	the list of available modules in your system include gcc, r, proj, perl, curl, geos and gdal

3.	you have cloned the AnAdde/N-SDM repository in your working directory `git clone https://github.com/AnAdde/N-SDM.git .`

4.	you have installed the nsdm R package `remotes::install_github("AnAdde/N-SDM/scripts/nsdm-project/functions", upgrade=FALSE)`

5. [optional for running the N-SDM example] you have downloaded the 6GB zip file containing example species and covariate data in the ./data folder `curl -o ./data/nsdm-project.zip  https://drive.switch.ch/index.php/s/u7DTEE84oDH4f57/download`

6. [optional for running the N-SDM example] you have unzipped nsdm-project.zip in the data folder `unzip ./data/nsdm-project.zip -d ./data/nsdm-project/`




# Example N-SDM run

<img align="right" alt="N-SDM logo" src="https://github.com/AnAdde/N-SDM/blob/main/images/n-sdm_logo_v2.png" width="300px"/>

In this example N-SDM run we will model the current and future distributions of three species (Larix decidua, Capra Ibex and Cantharellus cibarius) at 100-m resolution across Switzerland.

## Study area

Following the spatially-nested framework of N-SDM, two levels of analysis (local and global) were considered in this example. The local-level area includes all of Switzerland. For the global-level area, we use a bounding box covering the European continent with coordinates 32,60N, 71,70N, 28,56W and 40,21E.

## Data

### Species data

Global-level species occurrence records were obtained from GBIF (https://www.gbif.org/). Local-level records were obtained from the Swiss Species Information Center InfoSpecies (www.infospecies.ch).

### Covariate data

We will use a suite of 453 candidate covariates (Adde et al. in prep) derived from 42 individual parameters and belonging to 6 main categories (bioclimatic, land use and cover, edaphic, topographic, population density, transportation and vegetation). Note that some of these covariates were calculated using focal windows (e.g., land use and cover) and others are temporally dynamic. Once the .zip data file has been unzipped, detailed information on the covariates can be found in `./data/nsdm-project/covariates/covariates.xlsx`. Only bioclimatic covariates (n=19) were used for fitting the global-level model and all the others were used for the local model (Figure 1).

## N-SDM settings

N-SDM settings can be edited by modifying the settings.csv file located in ./scripts/nsdm-project/main/settings. You must edit this file to make N-SDM settings compatible with your computing environment (e.g. paths, HPC account, partition etc.). Be careful when saving settings.csv to use “;” as delimiter. In this same directory, the param-grid.xslx file allows specifying the grid for hyperparameter tunning. The pre-filled expert-table.xslx file allows for expert-based prefiltering of taxon-specific candidate covariates. For more details on N-SDM settings, please have a look at the N-SDM manuscript.

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

N-SDM development has been conducted within the ECOSPAT lab https://www.unil.ch/ecospat/en/home.html

Financial support through the Action Plan of the Swiss Biodiversity Strategy via the Federal Office for the Environment and the www.Valpar.ch project is gratefully acknowledged

