<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/AnAdde/N-SDM">
    <img src="FIG.png" alt="Logo" width="300" height="300">
  </a>

  <h3 align="center">N-SDM</h3>

  <p align="center">
    Nested-Species Distribution Modelling pipeline
    <br />
  </p>
</div>

## About N-SDM

N-SDM is a species distribution modelling computer code, mainly written in R language, specifically developed for high performance computing (HPC) environments. A key characteristic of N-SDM is that it is built around a scale-nesting framework aimed at avoiding climatic niche truncation when doing future projection (Chevalier et al. 2021, REF, REF). Among other innovative SDM features, N-SDM is also equipped with a newly devised covariate selection procedure for selecting best predictors in high dimensional spaces of candidates (Adde et al. in prep), an algorithm-specific hyperparameter grid search for identifying best model parameter values (Vignali et al. 2020, REF, REF) and the ensemble of small models approach for modelling rare species (Breiner et al. 2015, Breiner et al. 2018, Habibzadeh & Ludwig 2019).



## Prerequisites

Before running N-SDM, you should make sure you have an appropriate set up with the following requirements:

•	you are working on a Linux HPC cluster equipped with the Slurm workload manager...

•	... with available modules: gcc; r; proj; perl; curl; geos; gdal

•	the following R packages are installed: c('data.table', 'stringi', 'stringr', 'plyr', 'readxl','writexl', 'parallel', 'sp', 'raster', 'rgdal', 'zoo','fst','tools', 'glmnet', 'gam', 'mgcv', 'randomForest', 'RRF', 'lightgbm', 'ranger', 'maxnet','caret', 'ROCR','ecospat', 'chron','ggpubr')

•	you have installed the nsdm R package to access the bank of custom functions `remotes::install_github("AnAdde/N-SDM/scripts/nsdm-project/functions")`

•	you have downloaded the AnAdde/N-SDM repository as a zip file (e.g. `curl -L http://github.com/AnAdde/N-SDM/master.zip`) and unzipped it in your working directory



## Getting started with an example N-SDM run

In this example N-SDM run we will model the current and future distributions of three species (Larix decidua, Capra Ibex and Cantharellus cibarius) at 100-m resolution across Switzerland by using a suite of more than 100 candidate covariates.

### Study area

Following the spatially-nested framework of N-SDM, two levels of analysis(local and global) were considered in this exmaple. The local-level area includes all of Switzerland. For the global-level area, we use a bounding box covering the European continent with coordinates 32,60N, 71,70N, 28,56W and 40,21E.

### Data

Data required for this example N-SDM run can be downloaded from https://drive.switch.ch/index.php/s/u7DTEE84oDH4f57 and unziped in N-SDM/data>/nsdm-project/.

#### Species data

Global-level species occurrence records were obtained from GBIF (https://www.gbif.org/). Local-level records aggregated at a 100-m spatial resolution were obtained from the Swiss Species Information Center InfoSpecies (www.infospecies.ch).

#### Covariate data

We used a suite of 472 candidate covariates (Adde et al. in prep) derived from 42 individual parameters and belonging to 6 main categories (bioclimatic, land use and cover, edaphic, topographic, population density, transportation and vegetation). Note that some of these covariates are calculated using focal windows (e.g., land use and cover) and others (e.g., bioclimatic) are temporally dynamics. Detailed information on these covariates can be found in the Appendix of the reference manuscript. Only bioclimatic covariates (n=19) were used for fitting the global-level model and all the others were used for the local model (Figure 1).

### N-SDM settings

The edit of N-SDM settings should be done by modifying the settings.csv file available in scripts/nsdm-project/main/settings. You must edit this file to make the settings compatible with your computing environment (e.g. paths, session owner, HPC account, partition etc.). Be careful when saving settings.csv to use “;” as delimiter. In this same directory, the param-grid.xslx file allows specifying the grid for hyperparameter tunning. The pre-filled expert-table.xslx file allows for expert-based prefiltering of taxon-specific candidate covariates.

### Running N-SDM

Position yourself at /scripts/nsdm-project/main, where the main N-SDM bash file (nsdm.sh) is stored. We encourage you running N-SDM in a background no hangup mode to prevent the command from being aborted automatically if logging out or exiting the shell, such as: `nohup bash nsdm.sh &`.


## Contributing

If you have a suggestion that would make N-SDM better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Thanks!



## Contact

Antoine Adde - aadde@unil.ch

Project Link: [https://https://github.com/AnAdde/N-SDM](https://https://github.com/AnAdde/N-SDM)



## Acknowledgments


