## Different median methylation levels of significant genes for heatmap ##

library(dplyr)
library(tibble)
library(tidyr)
library(broom)
library(pheatmap)

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
  column_to_rownames('sample')
  
# plot the heatmap
heatr <- t(med.glm[,4:ncol(med.glm)]) # just the data, adjust the group as needed
heatr <- t(scale(t(heatr))) # scale before pheatmap because pheatmap scaling doesn't separate color enough
heatr.gs <- select(med.glm, temperature, treatment) # groups

heatr.gs$treatment <- factor(heatr.gs$treatment, levels = c("ctrl", "bmc", "path", "both")) # lock treatment levels

map <- pheatmap(heatr,
         color = rev(colorRampPalette(RColorBrewer::brewer.pal(10, 'RdYlBu'))(256)),
         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'ward.D2', 
         cutree_rows = 1, cutree_cols = 2,
         show_rownames = FALSE,
         angle_col = 45,
         border_color= NA,
         annotation_col = heatr.gs)

# Perform hierarchical clustering on columns to visualize tree only
hc_col <- hclust(dist(t(heatr)), method = "ward.D2")
# Plot the column dendrogram
plot(hc_col, main = "Column Dendrogram Only", xlab = "Columns", sub = "",
     hang = 0.3)
# Get the range of the heights
height_range <- range(hc_col$height)

#
