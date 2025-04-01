# HOGRCDM_V2

## Overview

**HOGRCDM_V2** is an R package accompanying the paper *"Exploratory General-response Cognitive Diagnostic Models with Higher-order Structures."* This package implements the algorithm described in Section 4 of the paper. Although the paper focuses on an exploratory framework, this project has been extended to support confirmatory cases as well. Users can flexibly choose between exploratory and confirmatory settings depending on their analysis requirements.

The package supports multiple measurement models, including:
- **Main-effect**, **All-effect**, and **DINA** models.

Additionally, it supports a variety of bottom layer distributions:
- **Lognormal** and **Gamma** for positive continuous responses. These can also be modified to other normal or transformed normal distributions for continuous responses with different value spaces.
- **Poisson** distribution for count data.
- **Logistic** distribution for binary data.

## Project Structure

The project consists of the following directories:

- **R/**: Contains all the core functions necessary to implement the algorithm.
- **man/**: Reserved for documentation of the primary functions located in the `R` folder. This folder is currently empty but will be populated with help files soon.
- **LoadValue/**: Includes true parameters and other values used in trial examples.
- **MyExample/**: Contains example scripts illustrating the usage of the algorithm for both confirmatory and exploratory cases.

## Script Naming Conventions

- Scripts with `"CC"` in the filename refer to confirmatory cases.
- Scripts with `"CE"` in the filename refer to exploratory cases, as discussed in the paper.

All types of measurement models and distributions are supported, and each script has been successfully tested on the owner's machine.

## Main Functions

The primary functions are located in the **R/** folder. Detailed help files for each function will be included in the **man/** folder soon. Below is a brief description of the key functions:

- **LoadingPackages**: Installs and loads all the necessary R packages.
- **GenerateData**: Generates data for simulation purposes.
- **INI_HOGRCDM_Main**: Initializes the model.
- Functions prefixed with `"FIT"`: These are the main functions for running the algorithm, with filenames following the format:
  
  `"FIT" + "setting (exploratory/confirmatory)" + "distribution"`

### Examples:
  
- For fitting a model with data following a log-normal distribution in a confirmatory setting, use `FIT_CC_Lognormal`.
- For exploratory setting with log-normal distribution, use `FIT_CE_Lognormal`.
- `FIT_CC_Poisson` is used for count data following a Poisson distribution.
- `FIT_CC_LLM` is used for binary data following a Bernoulli distribution.

### Data:
- The RData file named "DataAreBooklet1" contains the dataset used in the Real Data Analysis section of the paper.
- The TIMSS 2019 data is available at: https://timssandpirls.bc.edu/timss2019/
- For more information on the data structure, please see Section 7 of the paper.
