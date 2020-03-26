require(tidyverse)
require(GenABEL)
require(ranger)
library(sjmisc)

geno_pheno <- load.gwaa.data(phenofile = "analysis/pheno/phs000963.pheno_MACE.txt", 
                             genofile = "phs000963.remap.af01.nodup.autosomes.genabel")


phdata(geno_pheno)$MACE_EVENT <- factor(phdata(geno_pheno)$MACE_EVENT)

Y = phdata(geno_pheno)$MACE_EVENT

w <- 1/table(Y)
w <- w/sum(w)

weights <- rep(0, length(Y))
weights[Y == 0] <- w['0']
weights[Y == 1] <- w['1']

rg.mace <- ranger(MACE_EVENT ~ .,
                  case.weights = weights,
                  geno_pheno,
                  importance = 'impurity_corrected',
                  num.trees = 1000,
                  mtry = round(0.333 * geno_pheno@gtdata@nsnps),
                  write.forest = FALSE)
#
mace.importance <- data.frame(as.list(importance(rg.mace)))
mace <- t(mace.importance)
write.csv(data.frame(mace), "importance.csv", row.names = TRUE)
#
mace_p <- importance_pvalues(rg.mace)
#
write.csv(data.frame(mace_p), "importance_p.csv", row.names = TRUE)