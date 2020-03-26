# Generation of figures and some table modifications
# Note that due to data restrictions, some raw data files referenced here are not available to be shared

require(tidyverse)
require(pander)
require(gridExtra)
require(party)
require(ggpubr)

# 1 Manhattan Plots  ----

# * 1.1 Load/Transform Data ----

assoc_data <- read.table("md/data/plink_mace.assoc", header = TRUE) %>%
  mutate(SNP = sapply(SNP, function(x) {str_replace(str_replace(x, "\\.", ""), "-", "")}))

importance_coarse <- read_csv("md/data/importance_p.csv") %>%
  rename(SNP = X1) %>%
  mutate(SNP = sapply(SNP, function(x) {str_replace(str_replace(x, "\\.", ""), "-", "")}))

importance_fine <- read_csv("md/data/var_r2vim.csv") %>%
  rename(SNP = X1) %>%
  mutate(SNP = sapply(SNP, function(x) {str_replace(str_replace(x, "\\.", ""), "-", "")})) %>%
  select(SNP, rel.vim.min)

full_data <- assoc_data %>%
  full_join(importance_coarse, by = "SNP") %>%
  full_join(importance_fine, by = "SNP") %>%
  filter(SNP != "sex")

# * 1.2 Prepare Manhattan Plots ----

nCHR <- length(unique(full_data$CHR))
CHR <- unique(full_data$CHR)

full_data$BPcum <- NA
s <- 0
nbp <- c()

for (i in CHR){
  nbp[i] <- max(full_data[full_data$CHR == i,]$BP)
  full_data[full_data$CHR == i,"BPcum"] <- full_data[full_data$CHR == i,"BP"] + s
  s <- s + nbp[i]
}

axis.set <- full_data %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

# * 1.3 Plink Manhattan Plot  ----

tiff("figures/plink_manhattan.tiff", units = "in", width = 12, height = 4, res = 300)

ggplot(full_data, aes(x = BPcum, y = -log10(P), color = as.factor(CHR))) +
  geom_point(alpha = 0.75, size = 0.25) +
  geom_point(data = full_data[full_data$rel.vim.min > 1,], 
             aes(y = -log10(P), shape = "circle"), color = "red", alpha = 1, size = 2, stroke = 0.5) +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 7)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log10(p)") + 
  theme_minimal() +
  theme( 
    legend.position="bottom",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  ) + 
  scale_shape_manual(name = "", 
                        values = c("circle"), labels = c("Selected Variant")) + 
  guides(color = FALSE)

dev.off()

# * 1.4 Ranger Manhattan Plot  ----

tiff("figures/ranger_manhattan.tiff", units = "in", width = 12, height = 4, res = 300)

ggplot(full_data, aes(x = BPcum, y = -log10(pvalue), color = as.factor(CHR))) +
  geom_point(alpha = 0.75, size = 0.25) +
  geom_point(data = full_data[full_data$rel.vim.min > 1,], 
             aes(y = -log10(pvalue), shape = "circle"), color = "red", alpha = 1, size = 2, stroke = 0.5) +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 7)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log10(p)") + 
  theme_minimal() +
  theme( 
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  ) + 
  guides(color = FALSE, shape = FALSE)

dev.off()

# * 1.5 Combined Manhattan Plots ----

tiff("figures/manhattans.tiff", res = 600, units = "in", height = 2, width = 3)
grid.arrange(
  grid::rasterGrob(tiff::readTIFF("figures/plink_manhattan.tiff")), 
  grid::rasterGrob(tiff::readTIFF("figures/ranger_manhattan.tiff")),
  ncol = 1
  )
dev.off()

png("figures/manhattans.png", res = 300, units = "in", height = 2, width = 3)
grid.arrange(
  grid::rasterGrob(tiff::readTIFF("figures/plink_manhattan.tiff")), 
  grid::rasterGrob(tiff::readTIFF("figures/ranger_manhattan.tiff")),
  ncol = 1
)
dev.off()

# 2 Ensemble pairs  ----

# * 2.1 Load Data ----

ensemble_pairs <- read_csv("md/data/pair_data.csv") %>% 
  select(-X1) %>% 
  mutate(v12_prob_on = v12_prob * 10000)

ensemble_pairs$p_upper <- mapply(FUN = function(x, y) {
  binom.test(c(x, 10000 - x), p = y, alternative = "greater")$p.value 
  },
  x = ensemble_pairs$v12_c,
  y = ensemble_pairs$v12_prob)

ensemble_pairs$p_lower <- mapply(FUN = function(x, y) {
  binom.test(c(x, 10000 - x), p = y, alternative = "less")$p.value 
},
x = ensemble_pairs$v12_c,
y = ensemble_pairs$v12_prob)


# * 2.2 Build Plot  ----

ggplot(ensemble_pairs[ensemble_pairs$p_upper >= 0.05 & ensemble_pairs$p_lower >= 0.05,], 
       aes(x = v12_prob * 10000, y = v12_c)) + 
  geom_point(color = "grey", size = 1.5) + 
  geom_point(data = ensemble_pairs[ensemble_pairs$p_upper < 0.05,],
             aes(x = v12_prob * 10000, y = v12_c, color = "darkred"), size = 1.5) +
  geom_point(data = ensemble_pairs[ensemble_pairs$p_lower < 0.05,],
             aes(x = v12_prob * 10000, y = v12_c, color = "darkblue"), size = 1.5) +
  theme_bw() + 
  geom_abline(slope = 1, intercept = 0, linetype = 2) + 
  xlim(c(0, 80)) + 
  ylim(c(0, 80)) + 
  xlab("Expected Paired Selection Frequency (N Trees)") + 
  ylab("Actual Paired Selection Frequency (N Trees)") +
  scale_colour_manual(name = "Probable Pair\nRelationship", 
                      values =c("darkblue"="darkblue", "darkred"="darkred"), labels = c("LD", "Epistasis"))

ggsave(last_plot(), filename = "figures/paired_selection.png", width = 20, height = 15, dpi = 300, units = "cm")
ggsave(last_plot(), filename = "figures/paired_selection.tiff", width = 20, height = 15, dpi = 600, units = "cm")
  
  
  
# 3 Ensemble results  ----

# * 3.1 Load Data ----



# 4 Decision Trees  ----

# * 4.1 Modified node_terminal function for party tree ----

node_terminal2 <- function(ctreeobj,
                           digits = 3,
                           abbreviate = FALSE,
                           fill = c("lightgray", "white"),
                           id = TRUE,
                           BASE_OR = 1)
{
  getLabel1 <- function(x) {
    if (!x$terminal) return(rep.int("", 2))
    nlab <- paste("n =", sum(x$weights))
    ylab <- paste("OR =", round((x$prediction[[2]] / x$prediction[[1]]) / BASE_OR, 2), sep ="")
    return(c(nlab, ylab))
  }
  
  maxstr <- function(node) {
    lab <- getLabel1(node)
    msl <- ifelse(node$terminal, "", maxstr(node$left))
    msr <- ifelse(node$terminal, "", maxstr(node$right))
    lab <- c(lab, msl, msr)
    return(lab[which.max(nchar(lab))])
  }
  
  nstr <- maxstr(ctreeobj@tree)
  
  ### panel function for simple n, Y terminal node labelling
  rval <- function(node) {
    fill <- rep(fill, length.out = 2)	
    
    node_vp <- viewport(x = unit(0.5, "npc"),   
                        y = unit(0.5, "npc"),   
                        width = unit(1, "strwidth", nstr) * 1.1,
                        height = unit(3, "lines"),
                        name = paste("node_terminal", node$nodeID, sep = ""))
    pushViewport(node_vp)
    
    lab <- getLabel1(node)
    
    grid.rect(gp = gpar(fill = fill[1]))
    grid.text(y = unit(2, "lines"), lab[1])
    grid.text(y = unit(1, "lines"), lab[2])
    
    if (id) {
      nodeIDvp <- viewport(x = unit(0.5, "npc"), y = unit(1, "npc"),
                           width = max(unit(1, "lines"), unit(1.3, "strwidth", as.character(node$nodeID))),
                           height = max(unit(1, "lines"), unit(1.3, "strheight", as.character(node$nodeID))))
      pushViewport(nodeIDvp)
      grid.rect(gp = gpar(fill = fill[2], lty = "solid"))
      grid.text(node$nodeID)
      popViewport()
    }
    upViewport()
  }
  return(rval)
}

# * 4.2 Prep Date ----

fix_name <- list(
  "rs7790976" = "TMEM178B rs7790976",
  "exm1609541" = "APOBEC3H rs139298",
  "rs927629" = "STMND1 rs927629",
  "rs3747945" = "VAV3-AS1 rs3747945",
  "rs17222326" = "TENM3 rs17222326",
  "rs9312547" = "HAND2-AS1 (3’ 157.42kb) rs9312547",
  "sex" = "sex"
)

get_haploreg_gene <- function(rsid, big = FALSE) {
  if (startsWith(rsid, "rs")) {
    data <- queryHaploreg(query = rsid, file = NULL, study = NULL, ldThresh = 1,
                          ldPop = "EUR", epi = "vanilla", cons = "siphy", genetypes = "gencode",
                          url = "https://pubs.broadinstitute.org/mammals/haploreg/haploreg.php",
                          timeout = 10, encoding = "UTF-8", verbose = FALSE)
    var_data <- data[data$query_snp_rsid == rsid & data$is_query_snp == 1,] %>% 
      dplyr::select(chr, pos_hg38, rsID, EUR, RefSeq_name, RefSeq_distance, RefSeq_direction)
    if (big == TRUE) {
      return(data)
    } else {
      if (var_data$RefSeq_distance[[1]] > 0) {
        return(paste(var_data$RefSeq_name[[1]], 
                     paste("(", var_data$RefSeq_direction, "' ", as.character(round((as.numeric(var_data$RefSeq_distance[[1]]) / 1000), 2)), "kb", ")", sep = ""), 
                     sep = " "))
      } else {
        return(var_data$RefSeq_name[[1]])
      }
    }
  } else {
    return("NA")
  }
}

final_networks <- readRDS("md/data/RDS/final_networks.rds")
geno_pheno_min <- readRDS("md/data/RDS/geno_pheno_min.rds")

final_networks[[1]] <- unique(c(final_networks[[1]], final_networks[[5]]))
final_networks[[5]] <- "NA"
final_networks[[7]] <- unique(c(final_networks[[7]], final_networks[[8]]))
final_networks[[8]] <- "NA"
final_networks <- final_networks[final_networks != "NA"]

get_gene <- function(snp, fix_list) {
  if (snp %in% names(fix_list)) {
    return(fix_list[snp][[1]])
  }
  snp <- str_replace(snp, "exm\\.", "")
  return(paste(get_haploreg_gene(snp), snp, sep = " "))
}

geno_pheno_final <- geno_pheno_min %>% 
  dplyr::select(c("MACE_EVENT", "sex", unique(unlist(final_networks)))) %>% 
  rename_at(.vars = vars(-c(1)), ~sapply(., function(x) get_gene(x, fix_name)))

final_networks_gene <- lapply(final_networks, function(x) sapply(x, function(x) {get_gene(x, fix_name)}))

# * 4.3 Make Figures ----

control <- ctree_control(teststat = "quad",
                         testtype = "Univariate",
                         mincriterion = 0.999, 
                         minsplit = 125, 
                         minbucket = 125,
                         stump = FALSE, 
                         maxsurrogate = 0,
                         mtry = 0, 
                         savesplitstats = TRUE, 
                         maxdepth = 0, 
                         remove_weights = FALSE)

BASE_OR <- (sum(geno_pheno_final$MACE_EVENT == 1) / nrow(geno_pheno_final))


make_tree <- function(df, vars) {
  formula <- as.formula(paste("MACE_EVENT~", 
                              paste(sapply(vars, function(x) {paste0("`", x, "`")}), 
                                    collapse = "+"), 
                              sep=""))
  output.tree <- ctree(
    controls = control,
    formula,
    data = df)
  return(output.tree)
}

make_plot <- function(x, index) {
  tiff(paste("figures/plot_", index, ".tiff", sep = ""), units = "cm", width = 30, height = 15, res = 300)
  plot(x, main = paste("Network ", index, sep = ""),
       terminal_panel = node_terminal2(x,
                                       digits = 3,
                                       abbreviate = FALSE,
                                       fill = c("lightgray", "white"),
                                       id = TRUE,
                                       BASE_OR = BASE_OR)
  )
  dev.off()
}

sapply(1:length(final_networks_gene), function(i) {
  tree <- make_tree(geno_pheno_final, final_networks_gene[[i]])
  make_plot(tree, i)
})

# * 4.3 Make Figures ----

tiff("figures/trees.tiff", res = 600, units = "cm", height = 10, width = 30)

grid.arrange(
  grid::rasterGrob(tiff::readTIFF("figures/plot_1.tiff")), 
  grid::rasterGrob(tiff::readTIFF("figures/plot_2.tiff")),
  grid::rasterGrob(tiff::readTIFF("figures/plot_3.tiff")),
  grid::rasterGrob(tiff::readTIFF("figures/plot_4.tiff")),
  grid::rasterGrob(tiff::readTIFF("figures/plot_5.tiff")),
  grid::rasterGrob(tiff::readTIFF("figures/plot_6.tiff")),
  ncol = 3
)

dev.off()


png("figures/trees.png", res = 300, units = "cm", height = 10, width = 30)

grid.arrange(
  grid::rasterGrob(tiff::readTIFF("figures/plot_1.tiff")), 
  grid::rasterGrob(tiff::readTIFF("figures/plot_2.tiff")),
  grid::rasterGrob(tiff::readTIFF("figures/plot_3.tiff")),
  grid::rasterGrob(tiff::readTIFF("figures/plot_4.tiff")),
  grid::rasterGrob(tiff::readTIFF("figures/plot_5.tiff")),
  grid::rasterGrob(tiff::readTIFF("figures/plot_6.tiff")),
  ncol = 3
)

dev.off()


# * 4.4 Ensemble Table ----

read_csv("data/ensemble_final.csv") %>% 
  filter(fdr < 0.05) %>% 
  full_join(read_csv("data/sig_pairs.csv"), by = c("Var1" = "v1", "Var2" = "v2")) %>% 
  mutate(`Variant 1` = sapply(Var1, FUN = function(x) {get_gene(x, fix_name)})) %>% 
  mutate(`Variant 2` = sapply(Var2, FUN = function(x) {get_gene(x, fix_name)})) %>% 
  mutate(fdr.x = ifelse(round(fdr.x, 3) == 0, "<0.001", round(fdr.x, 3))) %>% 
  mutate(fdr.y = ifelse(round(fdr.y, 3) == 0, "<0.001", round(fdr.y, 3))) %>% 
  select(`Variant 1`, `Variant 2`, `$FDR_{ensemble}$` = fdr.x, `$FDR_{interaction}$` = fdr.y) %>% 
  write_csv("tables/ensemble_pairs.csv")

# * 4.5 Ensemble Table ----

# This is made by filtering the bio results by *cardio* and *angio*

read_tsv("data/ingenuity_cvd.tsv", skip = 1) %>% 
  filter(`B-H p-value` < 0.05 & `# Molecules` > 1) %>% 
  select(`Diseases or Functions` = `Diseases or Functions Annotation`, Genes = Molecules, FDR = `B-H p-value`) %>% 
  mutate(Genes = sapply(
                                                      Genes, 
                                                      function(x) {
                                                        str_replace_all(x, ",", " ")
                                                        })) %>% 
  write_csv("tables/pathway.csv")
  

