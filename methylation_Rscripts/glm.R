## Generalized linear model to determine significant genes in pairwise comparisons ##

library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(broom)

# read in compiled median methylation file
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

# subset timepoint T2 to determine genes specifically driving differences between inoculation groups at T2
med.meths.2 <- med.meths.c %>%
  filter(grepl('T2', timepoint)) %>%
  select(-timepoint, -sample) %>%
  relocate('gene')

med.meths.2$treatment <- factor(med.meths.2$treatment, levels = c("ctrl", "bmc", "path", "both"))

# run the glm on each gene
genlm <- med.meths.2 %>%
  nest(data = -gene) %>%
  mutate(
    model = map(data, ~ {
      l <- glm(median_meths ~ temperature + treatment, data = ., family = gaussian())
      s <- tryCatch(step(l, trace = 0), error = function(e) e)
      if (inherits(s, "error") || !s$converged) {
        tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
      } else {
        tidy(s)
      }
    })
  ) %>%
  unnest_wider(model) %>%
  select(-data)

# filter the results to keep only the relevant comparisons
rel.glm <- genlm %>%
  unnest(c(term, estimate, std.error, statistic, p.value))

# pivot the table for better readability
rel.glm.wide <- pivot_wider(rel.glm, names_from = term, values_from = c(estimate, std.error, statistic, p.value))
clean.rel.glm <- arrange(rel.glm.wide, gene)

glm.ctrlpath <- clean.rel.glm %>%
  select(gene, p.value_treatmentpath) %>%
  filter(!is.na(p.value_treatmentpath)) %>%
  mutate(adj_pval = p.adjust(p.value_treatmentpath, method = 'BH')) %>%
  filter(adj_pval < 0.05) %>%
  distinct(gene)
  
write.table(glm.ctrlpath, "glm_converged_ctrl_path_T2.tsv", sep="\t", row.names=FALSE)

#
