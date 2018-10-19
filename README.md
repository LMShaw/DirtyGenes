# DirtyGenes Readme

R Code and Data for DirtyGenes. Explanation of files given below.

* **Dirichlet_LRT_test.R** R code providing a function to execute all of the tests outlined in the DirtyGenes paper: 
  * Chi-squared likelihood ratio test,
  * Randomization likelihood ratio test,
  * Goodness-of-fit testing to the Dirichlet distribution,
  * Power Analysis.  
* **User Guide.docx** User guide for the Dirichlet_LRT_test.R code
* **DirtyGenes_workspace.RData** R Workspace containing data frames for all of the data used in the DirtyGenes paper. This may be used as example data for testing the Dirichlet_LRT_test.R function and as an example of the data format required to run the code.
* **Raw Sheep Data.zip** Zip file containing all of the raw, previously unpublished sheep hooves data used in the DirtyGenes paper. Folder contains 24 zip files representing samples from 4 sheep over 3 seasons (a,b,c) with each .txt file name describing the sheep, season and whether ARG or bacteria profiles are being considered. Files contain counts for each gene for ARG data and bacteria at genus level for bacteria data as well as outlining the taxonomic breakdown.


