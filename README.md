# EMGS

This repository contains the R package EMGS (currently still developer's version), described in 

Zehang R Li and Tyler H McCormick. _An Expectation Conditional Maximization approach for Gaussian graphical models_, 2018 [arXiv](https://arxiv.org/abs/1709.06970) 

# Using the package
You can install the package with devtools
```
install_github("richardli/EMGS", subdir = "EMGS")
library(EMGS)
```

## First naive example
This example demonstrates the bias reduction of EMGS. It creates a plot (illustraion-emgs.pdf) under the figures/ directory.
```
setwd("codes/")
source("example1.R")
```

## Second example
This example demonstrates the informative priors. It creates a plot (structure.pdf) under the figures/ directory.
```
source("example2.R")
```

## Third example
This example demonstrates the missing data imputation for Burke Gilman Trail example. It creates a plot (burke.pdf) under the figures/ directory. For quicker illustration, this example models only two blocks by combining the pedestrians and bikes from both directions instead of modeling four blocks.
```
source("example3.R")
```



## Simulation studies
Running the full version of the simulation as described in the paper takes a long time and is recommended to be implemented on a cluster. Below is a low dimensional example.
```
source("simulation-sim.R")
source("simulation-example.R")
```


