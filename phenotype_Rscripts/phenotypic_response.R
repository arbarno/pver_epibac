#### Analysis of phenotypic responses of P. verrucosa in reponse to heat stress and bacterial inoculation ####
    
## PAM analysis

pam<-read.csv("pam_timepoints.csv", header=TRUE)

# statistics to compare sample groups 
hist(pam$t0)
shapiro.test(pam$t0) # W = 0.94394, p-value = 0.0969

# statistics to compare sample groups 
hist(pam$t1)
shapiro.test(pam$t1) # W = 0.8702, p-value = 0.001173

# statistics to compare sample groups 
hist(pam$t2)
shapiro.test(pam$t2) # W = 0.56448, p-value = 1.384e-08

# compare ambient to heat stress at all timepoints with wilcoxon (2/3 not normal)
wilcox.test(pam$t0 ~ pam$temp) # W = 112, p-value = 0.5585
wilcox.test(pam$t1 ~ pam$temp) # W = 256, p-value = 1.542e-06
wilcox.test(pam$t2 ~ pam$temp) # W = 248, p-value = 6.635e-06

# subset temperatures to get inoculation comparisons
pam.a <- subset(pam, temp == 'a')
pam.h <- subset(pam, temp == 'h')

kruskal.test(pam.a$t1, pam.a$inoculation) # Kruskal-Wallis chi-squared = 3.5291, df = 3, p-value = 0.317
kruskal.test(pam.a$t2, pam.a$inoculation) # Kruskal-Wallis chi-squared = 4.8673, df = 3, p-value = 0.1818
kruskal.test(pam.h$t1, pam.h$inoculation) # Kruskal-Wallis chi-squared = 2.3382, df = 3, p-value = 0.5052
kruskal.test(pam.h$t2, pam.h$inoculation) # Kruskal-Wallis chi-squared = 3.9709, df = 3, p-value = 0.2646

## Nubbin picture analysis

library(ggplot2)
library(dplyr)

pheno<-read.csv("phenotype_analysis.csv", header=TRUE)

# get just the means of the color and percent tissue for each condition
pheno.means<- pheno %>%
  group_by(timepoint, temperature, group) %>%
  summarize(across(c(color, percent, live_color), list(mean=mean, sd=sd), na.rm=TRUE))

# lock the order of the inoculation groups
pheno.means$group <- factor(pheno.means$group, levels = c('ctrl', 'bmc', 'path', 'both'))

# plot the color intensity
ggplot(pheno.means, aes(x = timepoint, y = color_mean, group = group)) +
  geom_line(aes(y = color_mean, color = group), linewidth = 1) +
  geom_errorbar(aes(ymin = color_mean - color_sd, ymax = color_mean + color_sd, color = group), width = 0.05, linewidth = 1) +
  scale_x_discrete(expand = c(0.01,0.01), limits=c('T0', 'T1', 'T2')) +
  scale_y_continuous(expand = c(0,0), limits = c(25,160)) +
  facet_wrap(vars(temperature), nrow = 2, ncol = 1, scales = 'free') +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(), panel.border = element_blank())+
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),
        axis.title=element_text(size=14, family = "sans", colour = "black"))+
  theme(legend.text = element_text(size=14, family = "sans", colour = "black"),
        legend.title = element_text(size=14, family = "sans", colour = "black")) +
  scale_color_discrete(name = 'Treatment', labels=c('Control', 'BMC', 'Pathogen', 'Both')) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0, hjust = 0.5))+
  theme(axis.title.x = element_blank()) +
  theme(strip.background = element_blank(), strip.text.x = element_text(size = 18, colour = 'black'), strip.placement = "outside") +
  labs(y='Color Intensity')

# plot percent live tissue
ggplot(pheno.means, aes(x = timepoint, y = percent_mean, group = group)) +
  geom_line(aes(y = percent_mean, color = group), linewidth = 1) +
  geom_errorbar(aes(ymin = percent_mean - percent_sd, ymax = percent_mean + percent_sd, color = group), width = 0.05, linewidth = 1) +
  scale_x_discrete(expand = c(0.01,0.01), limits=c('T0', 'T1', 'T2')) +
  scale_y_continuous(expand = c(0,0), limits = c(-5,110)) +
  facet_wrap(vars(temperature), nrow = 2, ncol = 1, scales = 'free') +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(), panel.border = element_blank())+
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),
        axis.title=element_text(size=14, family = "sans", colour = "black"))+
  theme(legend.text = element_text(size=14, family = "sans", colour = "black"),
        legend.title = element_text(size=14, family = "sans", colour = "black")) +
  scale_color_discrete(name = 'Treatment', labels=c('Control', 'BMC', 'Pathogen', 'Both')) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0, hjust = 0.5))+
  theme(axis.title.x = element_blank()) +
  theme(strip.background = element_blank(), strip.text.x = element_text(size = 18, colour = 'black'), strip.placement = "outside") +
  labs(y='Percent live tissue (%)')

# plot color intensity of live tissue
ggplot(pheno.means, aes(x = timepoint, y = live_color_mean, group = group)) +
  geom_line(aes(y = live_color_mean, color = group), linewidth = 1) +
  geom_errorbar(aes(ymin = live_color_mean - live_color_sd, ymax = live_color_mean + live_color_sd, color = group), width = 0.05, size = 1) +
  scale_x_discrete(expand = c(0.01,0.01), limits=c('T0', 'T1', 'T2')) +
  scale_y_continuous(expand = c(0,0), limits = c(25,135)) +
  facet_wrap(vars(temperature), nrow = 2, ncol = 1, scales = 'free') +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(), panel.border = element_blank())+
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),
        axis.title=element_text(size=14, family = "sans", colour = "black"))+
  theme(legend.text = element_text(size=14, family = "sans", colour = "black"),
        legend.title = element_text(size=14, family = "sans", colour = "black")) +
  scale_color_discrete(name = 'Treatment', labels=c('Control', 'BMC', 'Pathogen', 'Both')) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0, hjust = 0.5))+
  theme(axis.title.x = element_blank()) +
  theme(strip.background = element_blank(), strip.text.x = element_text(size = 18, colour = 'black'), strip.placement = "outside") +
  labs(y='Color Intensity of live tissue')

## statistics to compare samples 
pheno.0 <- subset(pheno, timepoint == 'T0')
pheno.1 <- subset(pheno, timepoint == 'T1')
pheno.2 <- subset(pheno, timepoint == 'T2')
hist(pheno.0$color)
shapiro.test(pheno.0$color) # W = 0.97535, p-value = 0.06708
shapiro.test(pheno.1$color) # W = 0.8433, p-value = 0.000298
shapiro.test(pheno.2$color) # W = 0.90991, p-value = 0.01118

# both color and percent are not normal, so run kruskal-wallis -> dunn.test
pheno %>%
  group_split(timepoint, temperature) %>%
  purrr::map(~kruskal.test(.$color, .$group))
# only [[6]] is signficant (which is T2_H): Kruskal-Wallis chi-squared = 8.4926, df = 3, p-value = 0.036860

pheno %>%
  group_split(timepoint, temperature) %>%
  purrr::map(~kruskal.test(.$percent, .$group))
# only [[6]] is signficant (which is T2_H): Kruskal-Wallis chi-squared = 10.081, df = 3, p-value = 0.01789

pheno %>%
  group_split(timepoint, temperature) %>%
  purrr::map(~kruskal.test(.$live_color, .$group))
# no significance = [[6]] (which is T2_H): Kruskal-Wallis chi-squared = 7.4779, df = 3, p-value = 0.05813

# subset T2_H to get the pairwise comparisons for color and percent
pheno.t2.h <- subset(pheno, timepoint == 'T2' & temperature == 'H')
FSA::dunnTest(data = pheno.t2.h, color ~ group, method="bh")
# ctrl - path Z = -2.30209304, p-unadjusted = 0.02132993, p-adjusted = 0.04265985
FSA::dunnTest(data = pheno.t2.h, percent ~ group, method="bh")
# ctrl - path Z = 3.0447037, p-unadjusted = 0.002329097, p-adjusted = 0.01397458

#
