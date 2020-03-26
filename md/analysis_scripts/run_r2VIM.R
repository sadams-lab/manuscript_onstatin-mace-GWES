require(tidyverse)
require(GenABEL)
require(ranger)
require(recipes)
library(sjmisc)
source("functions/r2vim.R")

geno_pheno <- load.gwaa.data(phenofile = "analysis/pheno/phs000963.pheno_MACE.txt", 
                             genofile = "phs000963.remap.af01.nodup.autosomes.genabel")


vars <- read_csv("importance_p.csv") %>% 
  mutate(var = X1) %>% 
  filter(pvalue < 0.01) %>% 
  filter(var != "sex")

geno_pheno_fine <- geno_pheno[,vars$var]

phdata(geno_pheno_fine)$MACE_EVENT <- factor(phdata(geno_pheno_fine)$MACE_EVENT)

Y = phdata(geno_pheno_fine)$MACE_EVENT

w <- 1/table(Y)
w <- w/sum(w)

weights <- rep(0, length(Y))
weights[Y == 0] <- w['0']
weights[Y == 1] <- w['1']

var_r2vim <- var.sel.r2vim(formula = MACE_EVENT ~ .,
                           data = geno_pheno_fine,
                           no.runs = 10,
                           factor = 6,
                           ntree = 10000,
                           mtry.prop = 0.333,
                           nodesize.prop = 0.1,
                           no.threads = 16,
                           method = "ranger",
                           type = "classification",
                           case.weights = weights)

write.csv(var_r2vim$info, "var_r2vim.csv", row.names = TRUE)
