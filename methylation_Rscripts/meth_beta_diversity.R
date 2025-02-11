## Beta diversity of median methylation levels in samples ##

library(dplyr)
library(tibble)
library(vegan)
library(ggplot2)

# read in median methylation tsv file to run nMDS
med.meths <- read.table("compiled_median_meths.tsv", sep = '\t', header = TRUE)

# retain genes with >= 5 methylated positions and is methylated
med.meths.f <- med.meths %>%
  filter(meth_pos >= 5) %>%
  select(-meth_pos) %>%
  filter(if_all(everything(), ~. != 0))

# make dataframe including all grouping variables
df.div <- med.meths.f %>%
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
  column_to_rownames('sample')

# run the nMDS
meth.pct.0 <- subset(df.div, timepoint == 'T0')  # substitute t0/t1/t2
meth.data.0 <- meth.pct.0[,4:ncol(meth.pct.0)] # separate out just the data
meth.group.0 <- meth.pct.0[,1:3] %>% # create grouping table with variable
  as_tibble(rownames = 'sample')
set.seed(1)  # because nmds uses random number generator, this locks the number for replication
dist.0 <- vegdist(meth.data.0, method = 'bray') # get the distances between points
nmds.0 <- metaMDS(dist.0, k = 2, trymax = 1000) # run the nMDS
print(nmds.0)
s0 <- scores(nmds.0) %>% # this get the x and y axis points
  as_tibble(rownames = 'sample') %>%
  inner_join(., meth.group.0, by ='sample') # attach the axis points to the grouping variables for ggplot 

s <- bind_rows(s0, s1, s2) %>% # get all conditions in one dataframe for ggplot
  select(sample, timepoint, temperature, treatment, everything())

# lock the order of the x axis
s$treatment <- factor(s$treatment, levels = c('ctrl', 'bmc', 'path', 'both'))

# plot the figure
ggplot(s, aes(x = NMDS1, y = NMDS2, color = temperature)) +
  geom_point(aes(shape = treatment, color = temperature), size = 5) +
  scale_shape_discrete(name = 'Treatment', labels=c('Control','Cobetia','Vibrio','Both')) +
  facet_wrap(vars(timepoint), nrow = 1, ncol = 3, scales = 'free') +
  stat_ellipse(show.legend = FALSE) +
  geom_vline(linetype="dotted", xintercept = 0) +
  geom_hline(linetype="dotted", yintercept = 0) +
  theme_bw() +
  theme(panel.grid.minor=element_blank())+
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),
        axis.title=element_text(size=14, family = "sans", colour = "black"))+
  theme(legend.position="bottom",
       legend.spacing.x = unit(0.3, 'cm'),
        legend.text = element_text(size=12, family = "sans", colour = "black"), 
        legend.title = element_text(size=14, family = "sans", colour = "black")) +
  theme(strip.background = element_blank(), strip.text.x = element_text(size = 18, colour = 'black'),
        strip.placement = "outside")
  
# determining the significance
adonis2(dist.0 ~ meth.group.0$temperature * meth.group.0$treatment, permutations = 1000) # substitute t0/t1/t2
anova(betadisper(dist.0, meth.group.0$temperature))
anova(betadisper(dist.0, meth.group.0$treatment))
pairwiseAdonis::pairwise.adonis(dist.0, meth.group.0$temperature, p.adjust.m = 'BH', perm = 1000)
pairwiseAdonis::pairwise.adonis(dist.0, meth.group.0$treatment, p.adjust.m = 'BH', perm = 1000)

# distance between centroids
usedist::dist_between_centroids(dist.1,1:15, 16:29) #1:15 = A, 16:29 = H, value = [1] 0.02916496
usedist::dist_between_centroids(dist.2,1:16, 17:31) #1:16 = A, 17:31 = H, value = [1] 0.0221721

#
