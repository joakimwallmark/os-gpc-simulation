To ensure reproducible code, the project uses [renv](https://github.com/rstudio/renv). When the project folder is opened in R-studio on your local machine, the required packages with the correct versions should be automatically installed in a local environment, without affecting your existing global R environment.

# Instructions for replicating simulation study
- Run 0-scenario-setup.R to setup the simulated scenarios.
- Run 1-generate-data.R to simulate data.
- Run 2-fit-models.R to fit IRT models to simulated data.
- Run 3-compute-fit.R to compute performance measures.
- The code in 4-produce-plots-tables.R can be used to produce plots and tables of simulations results.

