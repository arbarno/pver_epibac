## Sleuth analysis for RNAseq data ##

# data was first pseud aligned with kallisto (https://pachterlab.github.io/kallisto/)

## https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html
## https://pachterlab.github.io/sleuth_walkthroughs/boj/analysis.html
## https://pachterlab.github.io/sleuth_walkthroughs/pval_agg/analysis.html

# specify where kallisto results are
base_dir <- "kallisto_gene_models"

sample_id <- dir(file.path(base_dir,"results"))

# list of paths to the kallisto in each of the folders
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results", id))

# load auxiliary table with the experimental design
s2c <- read.table(file.path(base_dir, "rnaseq_metadata.txt"), header = TRUE, stringsAsFactors=FALSE)

library(dplyr)

s2c <- select(s2c, sample, temperature, inoculation, group) %>%
  arrange(factor(temperature, levels = c("A", "H")), factor(inoculation, levels = c("ctrl", "bmc", "path", "both")))

# change order of kal_dirs to fit the order of the s2c
kal_dirs <- kal_dirs[c("C10","C12","C5","C9","A10","A13","A7","A9","D10","D12","D3","D9","B11","B12","B13","B7",
                       "G1","G2","G6","G8","F10","F11","F4","F5","H13","H2","H5","H7","E1","E11","E4","E9")]
s2c <- mutate(s2c, path = kal_dirs)

t2g <- read.table("target_id_gene_id.txt", sep = '\t', header = TRUE)

library(sleuth)

s2c$inoculation <- as.factor(s2c$inoculation)
s2c$inoculation <- relevel(s2c$inoculation, ref = "ctrl") 

# Now the “sleuth object” can be constructed. 
# This requires four commands that (1) load the kallisto processed data into the object 
# (2) estimate parameters for the sleuth response error measurement (full) model 
# (3) estimate parameters for the sleuth reduced model, and (4) perform differential analysis (testing). 
# On a laptop the four steps should take about a few minutes altogether.
so <- sleuth_prep(s2c, ~temperature * inoculation, target_mapping = t2g, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE,
                  read_bootstrap_tpm = TRUE, transform_fun_counts = function(x) log2(x + 1), )

# set up the PC values so the plot shows the percents alongside the "PC1"/"PC2"
ppv <- plot_pc_variance(so, units = "tpm")
PCpc <- ppv$data$var
names(PCpc) <- paste0("PC", seq(1:length(PCpc)))
pc1 <- round(PCpc["PC1"], 2)
pc2 <- round(PCpc["PC2"], 2)

library(ggplot2)

# plot PCA
plot_pca(so, units = "tpm", color_by = 'group', point_size = 5, point_alpha = 1) +
  aes(shape = inoculation) + 
  scale_color_manual(name = "Group", values = c("#808080", "#0075DC","#FFCC99","#F0A3FF","#94FFB5","#8F7C00","#9DCC00","#C20088")) +
  scale_shape_manual(name = "Inoculation", values = c(16, 17, 3, 15)) +
  labs(x = paste("PC1 (", pc1, "%)", sep = ""), y = paste("PC2 (", pc2, "%)", sep = "")) +
  geom_vline(linetype="dotted", xintercept = 0) + # create a center dotted line for reference
  geom_hline(linetype="dotted", yintercept = 0) +
  theme_bw() +
  theme(panel.grid.minor=element_blank())+
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),
        axis.title=element_text(size=14, family = "sans", colour = "black"))

# calculate the beta diversity
library(vegan)

abund.matrix <- sleuth:::spread_abundance_by(so$obs_norm_filt, "tpm", so$sample_to_covariates$sample)
abund.matrix <- t(abund.matrix)
meta <- so[["sample_to_covariates"]]

adonis2(dist.matrix ~ temperature * inoculation, data = meta, permutations = 1000, by = "term")
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 1000
# 
# adonis2(formula = abund.matrix ~ temperature * inoculation, data = meta, permutations = 1000, by = "term")
#                           Df SumOfSqs      R2       F   Pr(>F)    
# temperature                1  0.11842 0.32069 20.6739 0.000999 ***
# inoculation                3  0.05638 0.15268  3.2809 0.000999 ***
# temperature:inoculation    3  0.05699 0.15434  3.3166 0.002997 ** 
# Residual                  24  0.13747 0.37229                     
# Total                     31  0.36925 1.00000                     

# run pairwise adonis on all the combinations to see which groups are different from each other
dist.matrix <- vegdist(abund.matrix, method = "bray")
anova(betadisper(dist.matrix, meta$inoculation))
pairwiseAdonis::pairwise.adonis(dist.matrix, meta$inoculation, perm = 1000, p.adjust.m = 'BH')
# no pairwise comparisons were significant!

#quality control metrics to visualize the count distributions
plot_group_density(so, use_filtered = TRUE, units = "tpm",
                   trans = "log", grouping = setdiff(colnames(so$sample_to_covariates), "sample"), offset = 1)

# we can do a full sleuth fit model with both temperature and inoculation as factors
so <- sleuth_fit(so, ~temperature * inoculation, 'full', which_var = 'obs_counts')

# next, we fit the reduced model with just the temperature
so <- sleuth_fit(so, ~temperature, 'no_inoc', which_var = 'obs_counts') # test for main effect of inoc
so <- sleuth_fit(so, ~inoculation, 'no_temp', which_var = 'obs_counts') # test for main effect of temp
so <- sleuth_fit(so, ~temperature + inoculation, 'no_interact', which_var = 'obs_counts') # test if inoc depends on temp

# and finally perform the tests
# by doing this, we will statistically test for genes affect by inoculation, while controlling for differences in temperature
# this will identify genes whose abundances are significantly better explained when inoculation is taken into account
so <- sleuth_lrt(so, 'no_inoc', 'full') # test effect of inoc
so <- sleuth_lrt(so, 'no_temp', 'full') # test effect of temp
so <- sleuth_lrt(so, 'no_interact', 'full') # test effect of interaction              

# In general, we can test models that are nested using the likelihood ratio test. 
# Viewing models which have been fit can be done using the models() function.
models(so)

# get the number of genes affected by each
lrt_inoc <- sleuth_results(so, 'no_inoc:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)
lrt_temp <- sleuth_results(so, 'no_temp:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)
lrt_interact <- sleuth_results(so, 'no_interact:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)

# How many genes are broadly significant for each?
sum(lrt_inoc$qval < 0.05, na.rm = TRUE) # 5255 genes
sum(lrt_temp$qval < 0.05, na.rm = TRUE) # 12165 genes
sum(lrt_interact$qval < 0.05, na.rm = TRUE) # 3690 genes

library(tibble)
library(tidyr)

# Perform hierarchical clustering on columns to visualize tree only
sleuth.heat <- sleuth_to_matrix(so, 'obs_norm', 'tpm')
sleuth.heat <- sleuth.heat %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  pivot_longer(-c('sample'),names_to = 'target_id', values_to = 'tpm') %>%
  
  right_join(trans.table, by = 'target_id', keep = FALSE) %>%
  select(sample, target_id, tpm) %>%
  semi_join(., trans.sig, by = 'target_id') %>%
  mutate(temperature = case_when(startsWith(sample, 'A') ~ 'A', startsWith(sample, 'B') ~ 'A', startsWith(sample, 'C') ~ 'A', startsWith(sample, 'D') ~ 'A',
                                 startsWith(sample, 'E') ~ 'H', startsWith(sample, 'F') ~ 'H', startsWith(sample, 'G') ~ 'H', startsWith(sample, 'H') ~ 'H'),
         inoculation = case_when(startsWith(sample, 'A') ~ 'bmc', startsWith(sample, 'B') ~ 'both', startsWith(sample, 'C') ~ 'ctrl', startsWith(sample, 'D') ~ 'path',
                               startsWith(sample, 'E') ~ 'both', startsWith(sample, 'F') ~ 'bmc', startsWith(sample, 'G') ~ 'ctrl', startsWith(sample, 'H') ~ 'path')) %>%
  select(sample, temperature, inoculation, everything()) %>%
  pivot_wider(id_cols = c('sample','temperature', 'inoculation'), names_from = 'target_id', values_from = 'tpm') %>%
  column_to_rownames('sample')

# plot the heatmap
heatr <- t(sleuth.heat[,3:ncol(sleuth.heat)]) # just the data, adjust the group as needed
heatr <- t(scale(t(heatr))) # scale before pheatmap
heatr.gs <- select(sleuth.heat, temperature, inoculation) # groups

heatr.gs$inoculation <- factor(heatr.gs$inoculation, levels = c("ctrl", "bmc", "path", "both")) # lock treatment levels

library(pheatmap)
  
hc <- hclust(dist(sleuth.heat), method = "ward.D")

# Plot the column dendrogram
plot(hc, main = "", xlab = "", sub = "", hang = -1)
plot(as.dendrogram(hc))

# Get the range of the heights
height_range <- range(hc$height)
# Plot the column dendrogram
plot(hc, main = "", xlab = "", sub = "", hang = -1)
plot(as.dendrogram(hc),horiz=T)
# Get the range of the heights
height_range <- range(hc$height)

#

## Get contrasts between all groups (up or down regulated)

#library(sleuth)

# The unique pairwise inoculation comparisons you actually need
# ref=ctrl -> bmc_vs_ctrl, path_vs_ctrl, both_vs_ctrl
# ref=bmc -> path_vs_bmc, both_vs_bmc (ctrl_vs_bmc is just bmc_vs_ctrl flipped)
# ref=path -> both_vs_path (others already covered)
# ref=both -> nothing new

# So only run 3 reference levels, and filter to the new comparisons each time
comparisons_per_ref <- list(
  ctrl = c("inoculationbmc", "inoculationpath", "inoculationboth", "temperatureH",
           "temperatureH:inoculationbmc", "temperatureH:inoculationpath", "temperatureH:inoculationboth"),
  
  bmc  = c("inoculationpath", "inoculationboth", "temperatureH:inoculationpath",
           "temperatureH:inoculationboth"),
  
  path = c("inoculationboth",
           "temperatureH:inoculationboth"))

relevel_results <- list()

for (ref in names(comparisons_per_ref)) {
  
  message("Running sleuth_prep with ref = ", ref)
  s2c$inoculation <- relevel(factor(s2c$inoculation), ref = ref)
  
  so_tmp <- sleuth_prep(s2c, ~temperature * inoculation,
                        target_mapping = t2g,
                        aggregation_column = 'ens_gene',
                        extra_bootstrap_summary = TRUE,
                        read_bootstrap_tpm = TRUE,
                        transform_fun_counts = function(x) log2(x + 1))
  
  so_tmp <- sleuth_fit(so_tmp, ~temperature * inoculation, 'full')
  
  for (coef in comparisons_per_ref[[ref]]) {
    so_tmp <- sleuth_wt(so_tmp, which_beta = coef, which_model = 'full')
    res    <- sleuth_results(so_tmp, coef, which_model = 'full', show_all = FALSE, pval_aggregate = FALSE)
    key    <- paste0("ref_", ref, "__", coef)
    relevel_results[[key]] <- res
  }
  
  rm(so_tmp); gc()
}

# Save all
for (nm in names(relevel_results)) {
  write.csv(relevel_results[[nm]], 
            file = paste0("DE_", nm, ".csv"), 
            row.names = FALSE)
}

# Helper to count sig DEGs from relevel_results format
count_deg_relevel <- function(tbl) {
  tbl %>% 
    filter(!is.na(qval), qval <= 0.05, abs(b) > 1) %>% 
    nrow()
}

# Build counts dataframe from relevel_results only
deg.counts <- data.frame(
  comparison = c(
    # pairwise
    "H_vs_A",
    "bmc_vs_ctrl",
    "path_vs_ctrl",
    "both_vs_ctrl",
    "path_vs_bmc",
    "both_vs_bmc",
    "both_vs_path",
    # interactions
    "interaction_bmc_vs_ctrl",
    "interaction_path_vs_ctrl",
    "interaction_both_vs_ctrl",
    "interaction_path_vs_bmc",
    "interaction_both_vs_bmc",
    "interaction_both_vs_path"
  ),
  DEG_count = c(
    count_deg_relevel(relevel_results$ref_ctrl__temperatureH),
    count_deg_relevel(relevel_results$ref_ctrl__inoculationbmc),
    count_deg_relevel(relevel_results$ref_ctrl__inoculationpath),
    count_deg_relevel(relevel_results$ref_ctrl__inoculationboth),
    count_deg_relevel(relevel_results$ref_bmc__inoculationpath),
    count_deg_relevel(relevel_results$ref_bmc__inoculationboth),
    count_deg_relevel(relevel_results$ref_path__inoculationboth),
    count_deg_relevel(relevel_results$`ref_ctrl__temperatureH:inoculationbmc`),
    count_deg_relevel(relevel_results$`ref_ctrl__temperatureH:inoculationpath`),
    count_deg_relevel(relevel_results$`ref_ctrl__temperatureH:inoculationboth`),
    count_deg_relevel(relevel_results$`ref_bmc__temperatureH:inoculationpath`),
    count_deg_relevel(relevel_results$`ref_bmc__temperatureH:inoculationboth`),
    count_deg_relevel(relevel_results$`ref_path__temperatureH:inoculationboth`)
  ),
  type = c(
    rep("pairwise", 7),
    rep("interactions", 6)
  )
)

#

## Plot the contrasts ##

library(ggplot2)
library(dplyr)

# You'll need up/down counts, not just total
count_deg_direction <- function(tbl) {
  tbl %>%
    filter(!is.na(qval), qval <= 0.05, abs(b) > 1) %>%
    summarise(
      up   =  sum(b > 1),
      down = -sum(b < -1)  # negative so bar goes downward
    )
}

deg_dir <- data.frame(
  comparison = c("Ambient vs. Heat Stress", "Control vs. Cobetia", "Control vs. Vibrio",
                 "Control vs. Both", "Cobetia vs. Vibrio", "Cobetia vs. Both",
                 "Vibrio vs. Both", "Interaction:\nControl vs. Cobetia",
                 "Interaction:\nControl vs. Vibrio", "Interaction:\nControl vs. Both",
                 "Interaction:\nCobetia vs. Vibrio", "Interaction:\nCobetia vs. Both",
                 "Interaction:\nVibrio vs. Both"),
  rbind(count_deg_direction(relevel_results$ref_ctrl__temperatureH),
    count_deg_direction(relevel_results$ref_ctrl__inoculationbmc),
    count_deg_direction(relevel_results$ref_ctrl__inoculationpath),
    count_deg_direction(relevel_results$ref_ctrl__inoculationboth),
    count_deg_direction(relevel_results$ref_bmc__inoculationpath),
    count_deg_direction(relevel_results$ref_bmc__inoculationboth),
    count_deg_direction(relevel_results$ref_path__inoculationboth),
    count_deg_direction(relevel_results$`ref_ctrl__temperatureH:inoculationbmc`),
    count_deg_direction(relevel_results$`ref_ctrl__temperatureH:inoculationpath`),
    count_deg_direction(relevel_results$`ref_ctrl__temperatureH:inoculationboth`),
    count_deg_direction(relevel_results$`ref_bmc__temperatureH:inoculationpath`),
    count_deg_direction(relevel_results$`ref_bmc__temperatureH:inoculationboth`),
    count_deg_direction(relevel_results$`ref_path__temperatureH:inoculationboth`)))

deg_dir$comparison <- factor(deg_dir$comparison, levels = deg_dir$comparison)
write.table(deg_dir, "deg_numbers_contrasts.txt", quote = FALSE, sep = '\t', row.names=FALSE, col.names = TRUE)

deg_long <- deg_dir %>%
  tidyr::pivot_longer(cols = c(up, down), names_to = "direction", values_to = "count")

p <- ggplot(deg_long, aes(x = comparison, y = count, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("up" = "#2e4057", "down" = "#edae49"), name = NULL,
                    labels = c("up" = "Upregulated", "down" = "Downregulated"),
                    limits = c("up", "down")) +
  scale_y_continuous(labels = function(x) abs(x), breaks = seq(-1250, 1000, 250)) +
  labs(x = NULL, y = "Number of DEGs") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.border = element_blank(),
    axis.text.x = element_text(size = 11, family = "sans", colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, family = "sans", colour = "black"),
    axis.title = element_text(size = 13, family = "sans", colour = "black"),
    legend.text = element_text(size = 12, family = "sans", colour = "black"),
    legend.title = element_text(size = 12, family = "sans", colour = "black"),
    legend.position = "right")
p

ggsave("DEG_counts.pdf", plot = p, height = 6, width = 12)

#
                  
## Pairwise comparisons to be used for volcano plots
                  
library(tidyr)

# above tells you whether there's a significant difference across any of the 8 groups, but it does not perform pairwise comparisons between groups.

# create a contrast object (wald test)
so <- sleuth_wt(so, which_beta = 'temperatureH', which_model = 'full')
so <- sleuth_wt(so, which_beta = 'inoculationbmc', which_model = 'full')
so <- sleuth_wt(so, which_beta = 'inoculationpath', which_model = 'full')
so <- sleuth_wt(so, which_beta = 'inoculationboth', which_model = 'full')

# the results of the test can be examined with sleuth_results
temp.table <- sleuth_results(so,'temperatureH','wt', which_model = 'full', pval_aggregate = FALSE)
temp.sig <- filter(temp.table, qval <= 0.05 & abs(b) > 1)

bmc.table <- sleuth_results(so,'inoculationbmc','wt', which_model = 'full', pval_aggregate = FALSE)
bmc.sig <- filter(bmc.table, qval <= 0.05 & abs(b) > 1)

path.table <- sleuth_results(so,'inoculationpath','wt', which_model = 'full', pval_aggregate = FALSE)
path.sig <- filter(path.table, qval <= 0.05 & abs(b) > 1)

both.table <- sleuth_results(so,'inoculationboth','wt', which_model = 'full', pval_aggregate = FALSE)
both.sig <- filter(both.table, qval <= 0.05 & abs(b) > 1)

genes.up <- temp.sig %>% # count how many genes are significantly upregulated (substitute for each comparison)
  filter(b > 0) %>%
  select(ens_gene)
genes.down <- temp.sig %>% # count how many genes are significantly upregulated (substitute for each comparison)
  filter(b < 0) %>%
  select(ens_gene)
                  
## Volcano plots

library(EnhancedVolcano)

temp.table <- temp.table %>% 
  filter(!is.na(b))

g <- EnhancedVolcano(path.table,
                lab = NA,
                title = 'Control vs. Path',
                subtitle = NULL,
                caption = NULL,
                legendPosition = "right",
                legendLabSize = 14,
                #lab = both_table$target_id,
                #selectLab = goi.vol$target_id,
                #max.overlaps = Inf,
                x = 'b',
                y = 'qval',
                pCutoff = 5e-2,
                FCcutoff = 1,
                #xlim = c(-7, 7),
                ylim = c(0,7.5),
                pointSize = 1.5,
                labSize = 6.0,
                colAlpha = 0.5,
                drawConnectors = TRUE,
                widthConnectors = 0.75) +
  scale_x_continuous(limits = c(-3.1,3.1), breaks=seq(-3,3,1.5))
g

#

## correlation between all methylation levels and expression levels across genome

library(dplyr)
library(stringr)
  
sleuth.matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm')
  
med.meths <- read.table("compiled_median_meths.tsv", sep = '\t', header = TRUE)

med.meths.2 <- med.meths %>%
  filter(grepl('T2', timepoint))

# removes the "_T2" trailing the sample names
med.meths.2$sample <- str_extract(med.meths.2$sample, "^[^_]+")

library(tibble)
library(tidyr)

# rename the sleuth samples to read the same as the median methylation files. i.e. A7 -> A07
sleuth.rna <- sleuth.matrix %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('sample')
sleuth.rna$sample <- ifelse(grepl("^[A-Z]\\d$", sleuth.rna$sample), 
                            paste0(substr(sleuth.rna$sample, 1, 1), str_pad(substr(sleuth.rna$sample, 2, nchar(sleuth.rna$sample)), width = 2, side = "left", pad = "0")),
                            sleuth.rna$sample)

# clean the sleuth rna file
sleuth.rna <- sleuth.rna %>%
  column_to_rownames('sample') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('target_id') %>%
  left_join(t2g, by = 'target_id') %>%
  select(-target_id, ens_gene, everything()) %>%
  column_to_rownames('ens_gene') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  pivot_longer(-c('sample'), names_to = 'gene', values_to = 'tpm')
sleuth.rna$tpm <- as.numeric(as.character(sleuth.rna$tpm))

# combine the sleuth rna table with the median methylation table to perform linear regression
meth.exp <- med.meths.2 %>%
  left_join(., sleuth.rna, by = c("sample", "gene")) %>%
  na.omit(.) %>%
  filter(tpm != 0)

library(ggplot2)

# plot the linear resgression with all of the points and the linear model trendline
k <- ggplot(meth.exp, aes(x = median_meths, y = log2(tpm + 1))) +
  geom_point(alpha = 0.9, color = "darkgrey", stroke = 0) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_classic() +
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),
        axis.title=element_text(size=14, family = "sans", colour = "black")) +
  labs(x = "Median Methylation (%)", 
       y = expression("Log"[2] ~ "(Transcripts per million (TPM) + 1)"))
k

# Fit the linear model to the pooled data
model <- lm(log2(tpm + 1) ~ median_meths, data = meth.exp)

# Extract R² and p-value
summary(model)

#

## Plot the linear individual linear regressions correlating signficant gene from methylation analysis with expression data

# read in compiled median methylation file and curated list of significant genes
med.meths <- read.table("compiled_median_meths.tsv", sep = '\t', header = TRUE)
glm.ctrlpath <- read.table("glm_converged_ctrl_path_T2.tsv", sep = '\t', header = TRUE)

# retain genes with >= 5 methylated positions and is methylated
med.meths.f <- med.meths %>%
  filter(meth_pos >= 5) %>%
  select(-meth_pos) %>%
  filter(if_all(everything(), ~. != 0))

# wrangle the dataframe to include all the grouping variables
med.meths.c <- med.meths.f %>%
  column_to_rownames('gene') %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  mutate(timepoint = case_when(endsWith(sample, '0') ~ 'T0', endsWith(sample, '1') ~ 'T1', endsWith(sample, '2') ~ 'T2'),
         temperature = case_when(startsWith(sample, 'A') ~ 'A', startsWith(sample, 'B') ~ 'A', startsWith(sample, 'C') ~ 'A', startsWith(sample, 'D') ~ 'A',
                                 startsWith(sample, 'E') ~ 'H', startsWith(sample, 'F') ~ 'H', startsWith(sample, 'G') ~ 'H', startsWith(sample, 'H') ~ 'H'),
         treatment = case_when(startsWith(sample, 'A') ~ 'bmc', startsWith(sample, 'B') ~ 'both', startsWith(sample, 'C') ~ 'ctrl', startsWith(sample, 'D') ~ 'path',
                               startsWith(sample, 'E') ~ 'both', startsWith(sample, 'F') ~ 'bmc', startsWith(sample, 'G') ~ 'ctrl', startsWith(sample, 'H') ~ 'path')) %>%
  select(sample, timepoint, temperature, treatment, everything()) %>%
  as.data.frame() %>%
  pivot_longer(-c('sample','timepoint', 'temperature', 'treatment'), names_to = 'gene', values_to = 'median_meths')

# subset the timepoint T2
med.meths.2 <- med.meths.c %>%
  filter(grepl('T2', timepoint))
  
# filter the median methylation file to only include the signficant genes
med.glm <- glm.ctrlpath %>% 
  inner_join(., med.meths.2, by ='gene') %>% 
  pivot_wider(id_cols = c('sample','timepoint', 'temperature', 'treatment'), names_from = 'gene', values_from = 'median_meths') %>%
  column_to_rownames('X') %>%
  rename(inoculation = treatment)
med.glm$sample <- str_extract(med.glm$sample, "^[^_]+")

# remove unnecessary "T2" column and clean data
med.glm <- med.glm %>%
  select(-timepoint) %>%
  pivot_longer(-c('sample', 'temperature', 'inoculation'), names_to = 'gene', values_to = 'median_meths')

# combine the sleuth rna table with the median methylation table to perform linear regression
meth.exp <- med.glm %>%
  left_join(., sleuth.rna, by = c("sample", "gene"))

# run the linear models to determine significant correlations
model.stats <- meth.exp %>%
  group_by(gene) %>%
  do({
    model <- lm(log2(tpm + 1) ~ median_meths, data = .)
    model_summary <- summary(model)
    
    tibble(
      r2 = model_summary$r.squared,
      pval = model_summary$coefficients[2, 4],  # p-value for the slope (median_meths)
      df = model_summary$df[2],                  # Degrees of freedom for the residuals
      f_stat = model_summary$fstatistic[1]       # F-statistic value
    )
  })

# adjust pvalues for mutliple comparisons
pvals <- model.stats %>%
  select(gene, pval) %>%
  column_to_rownames('gene')
qval <- pvals %>%
  mutate(bh = p.adjust(pval, method = 'BH'))

# plot the scatterplot containing the significant genes and linear model trendlines for each
l <- ggplot(meth.exp, aes(x = median_meths, y = log2(tpm + 1))) +
  geom_point(aes(color = gene), alpha = 0.2) +
  geom_smooth(method = "lm", aes(color = gene, fill = gene), se = TRUE, alpha = 0.2) +
  theme_classic() + 
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),
        axis.title=element_text(size=14, family = "sans", colour = "black")) +
  labs(x = "Median Methylation (%)", 
       y = expression("Log"[2] ~ "(Transcripts per million (TPM) + 1)")) +
  scale_y_continuous(expand = c(0,0), limits=c(0,11)) +
  scale_x_continuous(expand = c(0,0), limits=c(0,101))
l

#
