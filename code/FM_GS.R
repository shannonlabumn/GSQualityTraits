library(StageWise)
library(tidyverse)
library(gt)
library(Hmisc)


source("./code/GS_mod.R")
source("./code/CV_func.R")

## Pheno data analysis for each trait and GS model pipeline

traits <- c("yield","sg","roundness")
GS_results <- vector(list,3); names(GS_results) <- c("yield","sg","roundness")

for(trait in traits){
  GS_results <- GS_run(p="../output/Chips_allqcd.csv", g="../data/genoChipsFp21k.csv", model = "", t = trait)

}
