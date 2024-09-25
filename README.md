# fenrir_paper_code
Code for Scalable Inference for Bayesian Multinomial Logistic-Normal Dynamic Linear Models Paper

## Overview of this repository
This repository consists of all the code (besides our R package [Fenrir](https://github.com/manansaxena/fenrir)) required to produce results for our paper (Link to be added later).

**Directory Structure**
```
.
├── code
│   ├── cpp
│   ├── stan
│   ├── analysis.ipynb
│   ├── generate_runtime_jobs.py
│   ├── generate_simulated_data.R
├── data
│   ├── mallard
└── scripts
├── .gitignore
├── .lintr
├── LICENSE
└── README.md
└── requirements.txt
```

## Guidelines to run experiements

* Search for keyword **random_path** for all the places where you would need to give path's based on your machine.
* All the plots presented in the paper are created using the code provided in analysis.ipynb jupyter notebook.
* To generate simulated data based on the simulation model presented in the paper, use generate_simulated_data.R file.
* To generate results for MAP estimation, Uncertainity Quantification on simulated data and the creation of simulated data itself can be done using generate_runtime_jobs.py. This file has the following two main functions:
  * Create simulated data for different set of values of D (number of multinomial categories) and N (number of time points).
  * Create and submit jobs to HPC cluster for running the Fenrir (cpp) and Stan code.
  * Other finer details are provided in the file. To use this file, you would need to run the following command:
   ```python generate_runtime_jobs.py {vary} {cores}```
* The cpp directory consists of the code to run our Fenrir model. 
  * Optimizer.R file is used to get results on the simulated data.
  * All files starting with **mallard** are for analysis on Artificial Gut Microbiome Data as discussed in the paper which is also provided in this repository under the data directory.
  * mallard.R and mallard_velocity.R is similar to Optimizer.R file in terms of functionality but applied to the microbiome data for MLN-DLM and local trend MLN-DLM model respectively. 
  * Mallard files with **gibbs** at the end, provide the code of our gibbs sampler as described in the paper.
* The stan directory consists of the code to run the our optimized Stan implementation of the model.
  * functions.stan file consists of all the helper and key functions like calculating log probability of the model used across different stan model files.
  * model.stan file is for defining the MLN-DLM model when **W** hyperparamter is known. 
  * model_W.stan and model_W_velocity.stan are used for defining MLN-DLM and local trend MLN-DLM models when we want to infer the **W** hyperparmeter.
  * run_stan.py is used to run our stan model on simulated data where as the files starting with **mallard** are used to run experiments on the microbiome data.
* The scripts directory consists of scripts of jobs submitted to the HPC. For simulated data, they are created automatically by generate_runtime_jobs.py. I have provided two sample scripts to give an idea of how they would look like. For microbiome data results, I would recommend to create new scripts exactly similar to sample ones with correct paths. They can be submitted to the cluster using ```sbatch``` command. For, running Fenrir on microbiome data, you also have the option to run it directly from RStudio for ease of debugging and tracking.    
   
These are general guidelines on how to run code in this repository. For finer details, please refer to individual files. Importantly, please set the correct file/directory paths in the files inplace of placeholders currently there.


## R Package Dependencies

Ensure these packages are installed from CRAN by running the following commands in your R console:

```R
install.packages("tidyverse")
install.packages("lubridate")
install.packages("fido")
install.packages("MCMCpack")
install.packages("fenrir")
install.packages("coda")
install.packages("abind")
install.packages("driver")
install.packages("LaplacesDemon")
install.packages("compositions")
install.packages("doParallel")
install.packages("foreach")
```

## Python Dependencies

```pip install -r requirements.txt``` to install required libraries.

