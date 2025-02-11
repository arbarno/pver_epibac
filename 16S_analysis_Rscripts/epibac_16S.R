#### Epibac bacteria community analysis ####

## Preparing the data

# read OTU data, calling this 'dat' 
# use the 'stringsAsFactors' argument to read OTU names as character
dat <- read.delim('ASV_counts_nosingle.txt', stringsAsFactors=F)

dat$timepoint <- factor(dat$timepoint, levels=c('T0','T1','T2'))
dat$temp <- factor(dat$temp, levels=c('A','H'))
dat$treatment <- factor(dat$treatment, levels=c('CTRL', 'BMC', 'PATH', 'BOTH'))
dat$group <- factor(dat$group, levels=c('T0_ACTRL', 'T0_ABMC', 'T0_APATH', 'T0_ABOTH','T0_HCTRL', 'T0_HBMC', 'T0_HPATH', 'T0_HBOTH',
                                        'T1_ACTRL', 'T1_ABMC','T1_APATH', 'T1_ABOTH','T1_HCTRL', 'T1_HBMC','T1_HPATH', 'T1_HBOTH',
                                        'T2_ACTRL', 'T2_ABMC','T2_APATH', 'T2_ABOTH','T2_HCTRL', 'T2_HBMC','T2_HPATH', 'T2_HBOTH'))
str(dat[, 1:5])

# read in the taxonomy, use 'stringsAsFactors=F' to keep text as 'character'
tax <- read.delim('tax_table_16S_nosingle.txt', stringsAsFactors=F)

## Creating phyloseq object

# We need to specify the different components to create a [phyloseq](https://github.com/joey711/phyloseq) object.

# load phyloseq, dplyr(to manipulate data)
library(phyloseq)
library(dplyr)

# extract the OTU data from 'dat' and convert to matrix format
# using 'starts_with' with 'select' allows us to chose all columns containing OTU counts
otus <- as.matrix(select(dat, starts_with('ASV')))

# rows need to be named with the sample names from 'dat', which we can do directly because they are in the same order
rownames(otus) <- dat$sample

# extract the sample data from from 'dat' and keep as dataframe format (mix of character and numeric variables)
samps <- select(dat, sample:group)

# rows need to be named with the sample names from 'dat', which we can do directly because they are in the same order
rownames(samps) <- dat$sample

# extract the taxonomy info for each OTU from 'tax' and convert to matrix format
taxonomy <- as.matrix(select(tax, kingdom:genus))

# rows need to be named with the OTU names from 'tax', which we can do directly because they are in the same order
rownames(taxonomy) <- tax$sample

# merge the three objects into a single phyloseq object using 'merge_phyloseq'
# each function nexted within the call to 'merge_phyloseq' creates a special object for that type of data
# because of how the OTU matrix is oriented, we have to specify 'taxa_are_rows=F' 
phy <- merge_phyloseq(otu_table(otus, taxa_are_rows=F), 
                sample_data(samps),
                tax_table(taxonomy))

# remove extra objects to keep workspace tidy
rm(otus, samps, taxonomy)

# remove archaea/eukarya
phy <- subset_taxa(phy, kingdom=="k_Bacteria")

# remove OTUs with abundance zero
phy <- prune_samples(sample_sums(phy)>0, phy)
                     
# remove OTUs that hit to Chloroplast or Mitochondria (because they are showing eukarya)
phy <- subset_taxa(phy, order != "o_Chloroplast")
phy <- subset_taxa(phy, family != "f_Mitochondria")

## Alpha diversity

library(microbiome)
library(ggpubr)

# Rarefaction using the phyloseq function rarefy_even_depth
set.seed(1)
phy.rarefied<- rarefy_even_depth(phy, sample.size = min(sample_sums(phy)), replace=FALSE, rngseed = TRUE)
phy
phy.rarefied

# we can look at the alpha diversity in each of the samples
set.seed(1)
rich<-estimate_richness(phy.rarefied)
even<-evenness(phy.rarefied)

# subset just the high temps to do 4 plots showing alpha diversity progression over time
phy.h <- subset_samples(phy.rarefied, temp=="H")
rich.h <- estimate_richness(phy.h)
even.h <- evenness(phy.h)

# we can plot at the alpha diversity in each of the samples
l<-plot_richness(phy.h, x="group", measures=c("Shannon")) +
  geom_boxplot(aes(fill=(sample_data(phy.h)$group)), outlier.shape = NA) +
  geom_point(size=1)+
  facet_wrap(~treatment, scales = "free_x", nrow = 1, strip.position = 'top',
             labeller = labeller(treatment = c('CTRL' = 'Control','BMC' = 'Cobetia',
                                           'PATH' = 'Vibrio','BOTH' = 'Both'))) +
  scale_x_discrete(labels=c('T0','T1','T2')) +
  theme_bw() +
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),axis.title=element_text(size=14, family = "sans", colour = "black"))+
  theme(legend.position="none") + 
  theme(axis.text.x=element_text(hjust=0.5,vjust=1)) +
  theme(axis.title.x = element_blank(), strip.text = element_blank()) +
  scale_y_continuous(name="Shannon diversity", limits=c(2,5.7), breaks = seq(2,6,1)) +
  theme(strip.placement = 'outside', strip.background = element_blank(), strip.text = element_text(size=14, family = "sans", colour = "black"))
l

ggsave("shannon_heats.pdf", plot = l, height = 5, width = 8)

# plot the evenness
k<-ggplot(even.h,aes(x=(sample_data(phy.h)$group), y=pielou,))+
  geom_boxplot(aes(fill=(sample_data(phy.h)$group)), outlier.shape = NA) +
  geom_point(size=1)+
  theme_bw()+
  facet_wrap((sample_data(phy.h)$treatment), scales = "free_x", nrow = 1, strip.position = 'top',
             labeller = labeller(treatment = c('CTRL' = 'Control','BMC' = 'Cobetia',
                                               'PATH' = 'Vibrio','BOTH' = 'Both'))) +
  scale_x_discrete(labels=c('T0','T1','T2')) +
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),axis.title=element_text(size=14, family = "sans", colour = "black"))+
  theme(legend.position="none")+ 
  theme(axis.text.x=element_text(hjust = 0.5, vjust = 1))+
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(name="Evenness", limits=c(0.3,0.8), breaks = seq(0.3,0.8,0.1)) +
  theme(strip.placement = 'outside', strip.background = element_blank(), strip.text = element_text(size=14, family = "sans", colour = "black"))
k

ggsave("evenness_heats.pdf", plot = k, height = 5, width = 8)

# statistics for comparing across timepoints
hist(rich.h$Shannon, main="Shannon index", xlab="")
ggqqplot(rich.h$Shannon, main="Shannon index", xlab="")
shapiro.test(rich.h$Shannon) ##Not normal W = 0.88128, p-value = 0.000165

# subset for each treatment
phy.h.ctrl <- subset_samples(phy.h, treatment=="CTRL")
rich.h.ctrl <-estimate_richness(phy.h.ctrl)
even.h.ctrl <- evenness(phy.h.ctrl)
phy.h.bmc <- subset_samples(phy.h, treatment=="BMC")
rich.h.bmc <-estimate_richness(phy.h.bmc)
even.h.bmc <- evenness(phy.h.bmc)
phy.h.path <- subset_samples(phy.h, treatment=="PATH")
rich.h.path <-estimate_richness(phy.h.path)
even.h.path <- evenness(phy.h.path)
phy.h.both <- subset_samples(phy.h, treatment=="BOTH")
rich.h.both <-estimate_richness(phy.h.both)
even.h.both <- evenness(phy.h.both)

kruskal.test(even.h.ctrl$pielou ~ sample_data(phy.h.ctrl)$timepoint) # substitute richness/evenness
FSA::dunnTest(even.h.ctrl$pielou ~ sample_data(phy.h.ctrl)$timepoint, 
              data = sample_data(phy.h.ctrl)$timepoint, method = "bh") # substitute richness/evenness

## Abundance analysis barplots

library(ggplot2)

# create a dummy variable with both the categories to make a relative abundance plot
phy.dum<-phy
sample_data(phy.dum)$dummy <- paste0(sample_data(phy.dum)$timepoint, sample_data(phy.dum)$temp, sample_data(phy.dum)$treatment)
phym = merge_samples(phy.dum, "dummy")

# repair the variables that you just destroyed
sample_data(phym)$timepoint <- levels(sample_data(phy.dum)$timepoint)[get_variable(phym, "timepoint")]
sample_data(phym)$temp <- levels(sample_data(phy.dum)$temp)[get_variable(phym, "temp")]
sample_data(phym)$treatment <- levels(sample_data(phy.dum)$treatment)[get_variable(phym, "treatment")]

# to do relative abundance (out of 100%), we need to also transform the counts to percents
phy.prop <- transform_sample_counts(phym, function(otu) 100 * otu/sum(otu))

# convert phyloseq object to dataframe to make figures in ggplot2
phy.propdf<-psmelt(phy.prop)

# keep the same levels on the x axis
phy.propdf$genus <- factor(phy.propdf$genus,levels=rev(unique(phy.propdf$genus)))

# filter to just the top 20 taxa at the family level.

# agglomerate at family level
phyg <- tax_glom(phy.prop, "family", NArm = FALSE)

# get just the top 20 order in the samples
top20otus = names(sort(taxa_sums(phyg), TRUE)[1:20])
taxtab20 = cbind(tax_table(phyg), phy20 = NA)
taxtab20[top20otus, "fam20"] <- as(tax_table(phyg)[top20otus, "family"], 
    "character")
tax_table(phyg) <- tax_table(taxtab20)

# prune samples
phyg20 = prune_taxa(top20otus, phyg)

# convert phyloseq object to dataframe to make figures in ggplot2
phyg20df<-psmelt(phyg20)

# keep the same levels on the x axis
phyg20df$phy20 <- factor(phyg20df$phy20,levels=rev(unique(phyg20df$phy20)))
phyg20df$Sample <- factor(phyg20df$Sample, 
                          levels = c("T0ACTRL", "T0ABMC", "T0APATH", "T0ABOTH", "T0HCTRL", "T0HBMC", "T0HPATH", "T0HBOTH", 
                                     "T1ACTRL", "T1ABMC", "T1APATH", "T1ABOTH", "T1HCTRL", "T1HBMC", "T1HPATH", "T1HBOTH", 
                                     "T2ACTRL", "T2ABMC", "T2APATH", "T2ABOTH", "T2HCTRL", "T2HBMC", "T2HPATH", "T2HBOTH"))

# you can define the colors
colours <- c("#808080", "#0075DC","#FFCC99","#F0A3FF","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405",
             "#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FF8000", "#993F00","#4C005C","#2BCE48")

# plot
n <- ggplot(phyg20df, aes(x=Sample, fill = phy20, y=Abundance)) +
  geom_bar(aes(color=phy20), stat = 'identity', position = 'stack', colour = NA) +
  facet_wrap(~timepoint, nrow = 1, scales = "free_x", strip.position = 'top') +
  scale_x_discrete(labels=c('A_Control','A_Cobetia','A_Vibrio','A_Both',
                            'H_Control','H_Cobetia','H_Vibrio','H_Both')) +
  theme_bw() +
  theme(axis.text=element_text(size=12, family = "sans", colour = "black"),axis.title=element_text(size=14, family = "sans", colour = "black"))+
  theme(legend.text = element_text(size=12, family = "sans", colour = "black"), legend.title = element_text(size=12, family = "sans", colour = "black")) +
  scale_fill_manual(name = "Top 20 families", values = colours) +
  scale_y_continuous(name= 'Abundance (%)', expand = c(0,0), limits=c(0,101)) +
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))+
  theme(axis.title.x = element_blank()) +
  theme(strip.placement = 'outside', strip.text = element_text(size=14, family = "sans", colour = "black"))
n

ggsave("top20_families.pdf", plot = n, height = 6, width = 9)

-------------------------------------------------------------------------------------------------------------------

## Beta diversity

# Now we want to check the beta diversity by making a nMDS or PCoA (Eslam did PCoA), and run some tests.

# first load the ggplot2 library
library(ggplot2)

# we can just use a subset of the samples so that we don't clog the plot
phy.t0 = subset_samples(phy, timepoint == "T0")
#phy.t0 = subset_samples(phy.t0, sample != "H12")

phy.t1 = subset_samples(phy, timepoint == "T1")

phy.t2 = subset_samples(phy, timepoint == "T2")

# this is to normalize the abundances of the samples (here used as a percent)
phy.t0.norm = transform_sample_counts(phy.t0, function(x) 100 * x/sum(x))
phy.t1.norm = transform_sample_counts(phy.t1, function(x) 100 * x/sum(x))
phy.t2.norm = transform_sample_counts(phy.t2, function(x) 100 * x/sum(x))

# calculate the distance matrix
phy.dist.t0 = phyloseq::distance(phy.t0.norm, method="bray")
phy.dist.t1 = phyloseq::distance(phy.t1.norm, method="bray")
phy.dist.t2 = phyloseq::distance(phy.t2.norm, method="bray")

# calculate the ordination (using NMDS, PCoA, PCA) using distance (bray, euclidean)
ord0 <- ordinate(phy.t0.norm, method="PCoA", distance=phy.dist.t0)
ord1 <- ordinate(phy.t1.norm, method="PCoA", distance=phy.dist.t1)
ord2 <- ordinate(phy.t2.norm, method="PCoA", distance=phy.dist.t2)

# see the eigen value of the axes
plot_scree(ord0)
plot_scree(ord1)
plot_scree(ord2)

# plot
p <- plot_ordination(phy.t2.norm, ord2, color = "temp", shape = "treatment")+
  #facet_wrap('timepoint', nrow = 1, ncol = 3, scales = 'free') +
  geom_point(size=2, alpha=0.8)+
  #scale_shape_manual(values=c(15, 16, 17,18))+
  theme_bw()+
  #scale_color_manual(values = c('#d7191c','#fdae61','#abdda4','#2b83ba'))+
  stat_ellipse(aes(color=phy.t2.norm@sam_data[["temp"]], group = phy.t2.norm@sam_data[["temp"]]))
p

# calculate the beta diversity
library(vegan)
adonis2(phy.dist.t0 ~ sample_data(phy.t0.norm)$temp * sample_data(phy.t0.norm)$treatment, permutations = 1000)
adonis2(phy.dist.t1 ~ sample_data(phy.t1.norm)$temp * sample_data(phy.t1.norm)$treatment, permutations = 1000)
adonis2(phy.dist.t2 ~ sample_data(phy.t2.norm)$temp * sample_data(phy.t2.norm)$treatment, permutations = 1000)

# run pairwise adonis on all the combinations to see which groups are different from each other
anova(betadisper(phy.dist.t0, sample_data(phy.t0.norm)$temp))
pairwiseAdonis::pairwise.adonis(phy.dist.t0, sample_data(phy.t0.norm)$group, perm = 1000, p.adjust.m = 'BH')

anova(betadisper(phy.dist.t1, sample_data(phy.t1.norm)$temp))
pairwiseAdonis::pairwise.adonis(phy.dist.t1, sample_data(phy.t1.norm)$group, perm = 1000, p.adjust.m = 'BH')

anova(betadisper(phy.dist.t2, sample_data(phy.t2.norm)$temp))
pairwiseAdonis::pairwise.adonis(phy.dist.t2, sample_data(phy.t2.norm)$group, perm = 1000, p.adjust.m = 'BH')

    ## OTU differential abundance testing with edgeR ----

library(edgeR)

phy.trans = transform_sample_counts(phy, function(x){x/sum(x)})
hist(log10(apply(otu_table(phy.trans), 2, var)),
     xlab="log10(variance)", breaks=50)
#we will set up a variance threshold to remove large fraction of OTUs with low variance
#Here we’ve used an arbitrary but not-unreasonable variance threshold of 10-5. 
#It is important to keep in mind that this filtering is independent of our downstream test. 
#The sample classifications were not used.
varianceThreshold = 1e-5
keepOTUs = names(which(apply(otu_table(phy.trans), 2, var) > varianceThreshold))
phy = prune_taxa(keepOTUs, phy)
phy
head(otu_table(phy))

# create function phyloseq to edgeR

phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}

    ## edgeR ----

# subset all the combos to get the pairwise edgeRs
#phy.t1.actrl.abmc = subset_samples(phy, timepoint == "T1" & treatment %in% c("ACTRL", "ABMC"))
#phy.t1.actrl.apath = subset_samples(phy, timepoint == "T1" & treatment %in% c("ACTRL", "ABMC"))
#phy.t1.actrl.aboth = subset_samples(phy, timepoint == "T1" & treatment %in% c("ACTRL", "ABMC"))
#phy.t1.actrl.hctrl = subset_samples(phy, timepoint == "T1" & treatment %in% c("ACTRL", "ABMC"))
#phy.t1.actrl.hbmc = subset_samples(phy, timepoint == "T1" & treatment %in% c("ACTRL", "ABMC"))
phy.t2.hctrl.hbmc = subset_samples(phy, timepoint == "T2" & group %in% c("T2_HCTRL", "T2_HBMC"))
#I think no genes with qualifying differential abundances
phy.t2.hctrl.hpath = subset_samples(phy, timepoint == "T2" & group %in% c("T2_HCTRL", "T2_HPATH"))
phy.t2.hctrl.hboth = subset_samples(phy, timepoint == "T2" & group %in% c("T2_HCTRL", "T2_HBOTH"))
phy.t2.hbmc.hpath = subset_samples(phy, timepoint == "T2" & group %in% c("T2_HBMC", "T2_HPATH"))
phy.t2.hbmc.hboth = subset_samples(phy, timepoint == "T2" & group %in% c("T2_HBMC", "T2_HBOTH"))
phy.t2.hpath.hboth = subset_samples(phy, timepoint == "T2" & group %in% c("T2_HPATH", "T2_HBOTH"))

phy.t2.temp = subset_samples(phy, timepoint == "T2")

# Now let’s use our newly-defined function to convert the phyloseq data object phyloseq_transformed2 into an edgeR “DGE” data object, called dge.
dge.t2.temp = phyloseq_to_edgeR(phy.t2.temp, group="group")

# Perform binary test
et = exactTest(dge.t2.temp)

# Extract values from test results
tt = topTags(et, n=nrow(dge.t2.temp$table), adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.001
sigtab = res[(res$FDR < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy.t2.temp)[rownames(sigtab), ], "matrix")) 
# Error in dimnames(x) <- dn : 
# length of 'dimnames' [1] not equal to array extent
dim(sigtab)
#[1]
head(sigtab)
write.csv(sigtab, file = "table_foldchange_t2_temp.csv", na = "NA") 

    ## Fold-change barplots ----

# plot Here is a bar plot showing the log-fold-change, showing Genus and Phylum. Uses some ggplot2 commands.
fcp <- read.delim('foldchange_t2.txt', stringsAsFactors=F)

colours <- c("#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FF8000", "#993F00","#4C005C","#2BCE48",
              "#808080", "#0075DC","#FFCC99","#F0A3FF","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405",
             "#E495A5", "#ABB065", "#39BEB1", "#ACA4E2", "#009999")

p<-ggplot(
  data = fcp,
  aes(x = comparison, y = logFC, group = reorder(row.names(fcp), as.numeric(genus)))) + 
  geom_bar(stat = "identity", aes(fill = genus), position = position_dodge(preserve = "single")) +
  scale_fill_manual(values = colours) +
  ylab("Log Fold Change") + 
  coord_flip()
p

    ## Indicspecies ----

library(indicspecies)
library(stats)

# the data should be the same as above, but if not there, run command
dat.piv <- read.delim('genus_counts_nosingle_nochlor_nomito.txt', stringsAsFactors=F)

ind.sp = subset(dat.piv, timepoint == 'T1' & temp == 'H')
ind.sp.abund = ind.sp[,6:ncol(ind.sp)]
treat = ind.sp$treatment

# comparing GROUPS of potential indicators instead of individual species
#ind.sp.combo = combinespecies(ind.sp)
#dim(ind.sp.combo)
#indvalspcomb = multipatt(ind.sp.combo, groups, duleg = TRUE, control = how(nperm=999))
#summary(indvalspcomb, indvalcomp = TRUE)
#p.adjust(c$p, method = "BH")

library(data.table)
indisp<- multipatt(ind.sp.abund, treat, duleg = TRUE, func = "r.g", control = how(nperm=1000))
#extract table of stats
indisp.sign<-as.data.table(indisp$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
indtable <- indisp.sign[p.value.bh<=0.05 & stat>=0.75, ]

write.csv(indtable, file = "indicator_species_t1h.csv", na = "NA")  #export table

    ## multi_boxplot_w/_inoculums ----

# to do relative abundance (out of 100%), we need to also transform the counts to percents
phy.rel <- transform_sample_counts(phy, function(otu) 100 * otu/sum(otu))

#phy.rel.rueg <- subset_taxa(phy.rel, genus == "g_Ruegeria")
#phy.rel.inoc <- subset_taxa(phy.rel, genus == "g_Cobetia" | genus == "g_Vibrio")
#phy.rel.inoc <- subset_taxa(phy.rel, family == "f_Endozoicomonadaceae")

Taxa = c("ASV176", "ASV319")
# "ASV36", "ASV37", "ASV57" vibrio
# "ASV176", "ASV319" cobetia
allTaxa = taxa_names(phy.rel)
myTaxa <- allTaxa[(allTaxa %in% Taxa)]
phy.rel.inoc = prune_taxa(myTaxa, phy.rel)

phy.rel.inoc<-tax_glom(phy.rel.inoc, "genus", NArm = FALSE)

phy.reldf.inoc<-psmelt(phy.rel.inoc)

phy.reldf.inoc<-phy.reldf.inoc %>%
  filter(., timepoint != "T2")

# plot the boxplot
q<-ggplot(phy.reldf.inoc,aes(x=group, fill = genus, y=Abundance))+
  geom_boxplot(aes(fill=genus), outlier.shape = NA) +
  #scale_fill_manual(values=c("plum")) +
  geom_point(size=0.8,position=position_dodge(width=0.75)) +
  facet_wrap(~timepoint, scales = "free_x") +
  theme_bw()+
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))
  #scale_y_continuous(name="Abundance (%)", limits=c(0,4), breaks = seq(0,4,1))
q

# we can just use a subset of the samples so we can do stats within timepoints
phy.t0.rel.vib = subset_samples(phy.rel.inoc, timepoint == "T0")
phy.t0.rel.vib = subset_taxa(phy.t0.rel.vib, genus == "g_Vibrio")
phy.t0.rel.vib.a = subset_samples(phy.t0.rel.vib, temp == "A")
phy.t0.rel.vib.h = subset_samples(phy.t0.rel.vib, temp == "H")
phy.t1.rel.vib = subset_samples(phy.rel.inoc, timepoint == "T1")
phy.t1.rel.vib = subset_taxa(phy.t1.rel.vib, genus == "g_Vibrio")
phy.t1.rel.vib.a = subset_samples(phy.t1.rel.vib, temp == "A")
phy.t1.rel.vib.h = subset_samples(phy.t1.rel.vib, temp == "H")
phy.t2.rel.vib = subset_samples(phy.rel.inoc, timepoint == "T2")
phy.t2.rel.vib = subset_taxa(phy.t2.rel.vib, genus == "g_Vibrio")
phy.t2.rel.vib.a = subset_samples(phy.t2.rel.vib, temp == "A")
phy.t2.rel.vib.h = subset_samples(phy.t2.rel.vib, temp == "H")

phy.t0.reldf.vib<-psmelt(phy.t0.rel.vib)
phy.t0.reldf.vib.a<-psmelt(phy.t0.rel.vib.a)
phy.t0.reldf.vib.h<-psmelt(phy.t0.rel.vib.h)
phy.t1.reldf.vib<-psmelt(phy.t1.rel.vib)
phy.t1.reldf.vib.a<-psmelt(phy.t1.rel.vib.a)
phy.t1.reldf.vib.h<-psmelt(phy.t1.rel.vib.h)
phy.t2.reldf.vib<-psmelt(phy.t2.rel.vib)
phy.t2.reldf.vib.a<-psmelt(phy.t2.rel.vib.a)
phy.t2.reldf.vib.h<-psmelt(phy.t2.rel.vib.h)

# we can just use a subset of the samples so we can do stats within timepoints
phy.t0.rel.cob = subset_samples(phy.rel.inoc, timepoint == "T0")
phy.t0.rel.cob = subset_taxa(phy.t0.rel.cob, genus == "g_Cobetia")
phy.t0.rel.cob.a = subset_samples(phy.t0.rel.cob, temp == "A")
phy.t0.rel.cob.h = subset_samples(phy.t0.rel.cob, temp == "H")
phy.t1.rel.cob = subset_samples(phy.rel.inoc, timepoint == "T1")
phy.t1.rel.cob = subset_taxa(phy.t1.rel.cob, genus == "g_Cobetia")
phy.t1.rel.cob.a = subset_samples(phy.t1.rel.cob, temp == "A")
phy.t1.rel.cob.h = subset_samples(phy.t1.rel.cob, temp == "H")
phy.t2.rel.cob = subset_samples(phy.rel.inoc, timepoint == "T2")
phy.t2.rel.cob = subset_taxa(phy.t2.rel.cob, genus == "g_Cobetia")
phy.t2.rel.cob.a = subset_samples(phy.t2.rel.cob, temp == "A")
phy.t2.rel.cob.h = subset_samples(phy.t2.rel.cob, temp == "H")

phy.t0.reldf.cob<-psmelt(phy.t0.rel.cob)
phy.t0.reldf.cob.a<-psmelt(phy.t0.rel.cob.a)
phy.t0.reldf.cob.h<-psmelt(phy.t0.rel.cob.h)
phy.t1.reldf.cob<-psmelt(phy.t1.rel.cob)
phy.t1.reldf.cob.a<-psmelt(phy.t1.rel.cob.a)
phy.t1.reldf.cob.h<-psmelt(phy.t1.rel.cob.h)
phy.t2.reldf.cob<-psmelt(phy.t2.rel.cob)
phy.t2.reldf.cob.a<-psmelt(phy.t2.rel.cob.a)
phy.t2.reldf.cob.h<-psmelt(phy.t2.rel.cob.h)

# statistics within each timepoint
hist(phy.t1.reldf.vib$Abundance)
shapiro.test(phy.t1.reldf.vib$Abundance) # p>0.05 is normal, can run anova

kruskal.test(phy.t0.reldf.cob.h$Abundance,phy.t0.reldf.cob.h$treatment) # for non-normal data
FSA::dunnTest(phy.t0.reldf.cob.h$Abundance~phy.t0.reldf.cob.h$treatment, data = phy.t0.reldf.cob.h, method='bh')
#pairwise.wilcox.test(phy.t2.reldf.vib$Abundance,phy.t2.reldf.vib$treatment, p.adj='fdr')
