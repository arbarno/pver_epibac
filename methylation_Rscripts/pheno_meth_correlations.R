#### linear regressions and boxplots of individual genes found to have correlated methylation levels and phenotypic responses ####

## Correlations & linear regressions

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

# read in pct compiled median methylation file
med.meths <- read.table("compiled_median_meths.tsv", sep = '\t', header = TRUE)

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
         treatment = case_when(startsWith(sample, 'A') ~ 'bmc', startsWith(sample, 'B') ~ 'both', startsWith(sample, 'C') ~ 'ctrl', startsWith(sample, 'D') ~ 'path',
                               startsWith(sample, 'E') ~ 'both', startsWith(sample, 'F') ~ 'bmc', startsWith(sample, 'G') ~ 'ctrl', startsWith(sample, 'H') ~ 'path')) %>%
  pivot_longer(-c('sample','timepoint', 'treatment'), names_to = 'gene', values_to = 'median_meths')

# read in phenotypic csv data for left join
pheno <- read.csv("phenotype_analysis.csv", header = TRUE)

pheno.meth <- left_join( med.meths.c, pheno, by = 'sample')

# subset gene for regression (test all genes)
gene <- pheno.meth %>%
  filter(grepl('g8986$', gene))

# create model for plotting
gene.model <- lm(median_meths~color_intensity, data = gene)
summary(gene.model)

# compile all pvals and adjust pval using bh
lm_pvals <- read.csv("linear_model_regressions.csv", header = TRUE)
pval.bh <- lm_pvals %>% 
  mutate(adjusted_pval = p.adjust(pval, method = 'BH'))

gene$treatment <- factor(gene$treatment, levels = c("ctrl", "bmc", "path", "both")) # lock treatment levels

# plot linear regression
l <- ggplot(gene, aes(x = color_intensity, y = median_meths)) +
  geom_smooth(formula = 'y~x', method = 'lm', se = TRUE, level = 0.95,fullrange = T,
              color = 'black', linewidth = 0.8, alpha = 0.3) +
  geom_point(aes(color = timepoint, shape = treatment), size = 3) +
  scale_color_manual(values=c('#8d96a3','#edae49', '#2e4057')) +
  theme_bw() +
  ggtitle('Pver_gene_g8986') + # change based on gene plotted
  theme(plot.title = element_text(size=14, family = 'sans', color = 'black', hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(), panel.border = element_blank())+
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),
        axis.title=element_text(size=14, family = "sans", colour = "black"))+
  theme(legend.text = element_text(size=14, family = "sans", colour = "black"),
        legend.title = element_text(size=14, family = "sans", colour = "black")) +
  scale_shape_discrete(name = 'Treatment', labels=c('Control', 'BMC', 'Pathogen', 'Both')) +
  labs(x = 'Color Intensity',
       y = 'Median Methylation (%)')
l

#

## Median methylation boxplots of specific genes

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

# read in pct compiled median methylation file
med.meths <- read.table("compiled_median_meths.tsv", sep = '\t', header = TRUE)

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

# subset gene for regression (test all genes)
gene <- med.meths.c %>%
  filter(grepl('g8986$', gene)) %>%
  filter(grepl('T2', timepoint))

# lock the order of the x axis
gene$treatment <- factor(gene$treatment, levels = c('ctrl', 'bmc', 'path', 'both'))

# plot the boxplot
b<-ggplot(gene, aes(x = temperature, y = median_meths)) + # for color intensity
  geom_boxplot(aes(fill=treatment), colour= "black", lwd=0.5, fatten = 0.8) +
  geom_point(aes(group = treatment), size=1.5, position = position_dodge(width=0.75)) +
  scale_x_discrete(expand = c(0,0.4), labels=c('Ambient', 'Heat-stressed')) +
  theme_bw() +
  theme(plot.title = element_text(size=14, family = 'sans', color = 'black', hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(), panel.border = element_blank())+
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),
        axis.title=element_text(size=14, family = "sans", colour = "black"))+
  theme(legend.text = element_text(size=14, family = "sans", colour = "black"),
        legend.title = element_text(size=14, family = "sans", colour = "black")) +
  scale_fill_discrete(name = 'Treatment', labels=c('Control', 'BMC', 'Pathogen', 'Both')) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0, hjust = 0.5))+
  theme(strip.background = element_blank(), strip.text.x = element_text(size = 18, colour = 'black'), strip.placement = "outside") +
  labs(x = "", y="Methylation (%)")
b

#
