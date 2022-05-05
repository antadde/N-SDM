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

We used a new suite of 1,508 swiss-wide candidate covariates (Adde et al. in prep) derived from 157 individual parameters and belonging to 8 main categories (Appendix [put the detailed list is an annex]). Only bioclimatic covariates (n=19) were used for fitting the global-level model and only habitat covariates (n=1,489) were used for the local model (Figure 1). 



<!-- CONTRIBUTING -->
## Contributing

If you have a suggestion that would make N-SDM better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<!-- LICENSE -->
## License

Distributed under the XXX License. See `LICENSE.txt` for more information.



<!-- CONTACT -->
## Contact

Antoine Adde - aadde@unil.ch

Project Link: [https://https://github.com/AnAdde/N-SDM](https://https://github.com/AnAdde/N-SDM)



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments


