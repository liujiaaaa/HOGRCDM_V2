---
output:
  html_document: default
  pdf_document: default
---
README for HOGRCDM_V2 Project  

This repository contains an R Project named HOGRCDM_V2 for the paper "Exploratory General-response Cognitive Diagnostic Models with Higher-order Structures."  


This project includes the code for Algorithm 1 in Section 4 of the paper. While the paper focuses on an exploratory scenario, we have extended the project to support confirmatory settings as well, allowing users to flexibly choose between these two settings based on their specific needs.  


---Project Structure

--R Folder: This folder contains all the necessary functions required to implement the algorithm.

--man Folder: This folder is intended for the help files of the main functions in folder "R". It is currently empty but will be populated soon.

--LoadValue Folder: This folder contains the true parameters and other values that may be used when running the trial examples.

--MyExample Folder: This folder includes example scripts for running the algorithm in both confirmatory and exploratory cases.

--Script Naming Conventions:

-The R scripts with "CC" in the filename are for confirmatory cases.  

-The scripts with "CE" in the filename are for exploratory cases (as discussed in the paper).  

-All measurement model types and distributions are supported by this project, and all scripts have been successfully tested on the owner's computer.



---Main Functions
The primary functions used in the project are located in the R folder. Since detailed help files will be added to the man folder soon, a brief overview of the main functions is provided below:

--LoadingPackages: This function installs and loads all the required R packages.  


--GenerateData: This function generates data for simulations.  


--INI_HOGRCDM_Main: This function initializes the model.  


--Files with "FIT": These are the main functions for running the algorithms. They are named using the format:  


"FIT" + "setting (exploratory/confirmatory)" + "distribution"

-To fit a model with data following a log-normal distribution in a confirmatory setting, use the function FIT_CC_Lognormal.  

-To fit a model with data following a log-normal distribution in an exploratory setting (as considered in the paper), use FIT_CE_Lognormal.  

-FIT_CC_Poisson is for count data following a Poisson distribution.  

-FIT_CC_LLM is for binary data following a Bernoulli distribution.
