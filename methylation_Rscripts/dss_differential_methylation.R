## DSS differential methylation analysis ##

library(DSS)
require(bsseq)

path = file.path("dss_analysis/")

# read in the design file
design <- read.table("design.txt", sep = '\t', header = TRUE)
design$Temperature <- factor(design$Temperature,
                              levels = c("ambient", "heat"))
design$Inoculation <- factor(design$Inoculation,
                             levels = c("ctrl", "bmc", "both", "path"))

# read in all the methylation data from the t2 files
meth_list <- lapply(seq_len(nrow(design)), function(i) {
  sample_id <- design$X[i]
  inoc <- design$Inoculation[i]
  read.table(file.path(path, inoc, paste0(sample_id, "_dss_input.txt")), sep = "\t", header = TRUE)
})

names(meth_list) <- design$X

# make the bs object with the files
BSobj <- makeBSseqData(meth_list, sampleNames = design$X)

# observe the model matrix
X = model.matrix(~Temperature+Inoculation, design)

# fit the model (similar to multiple regression); model only needs to be fit once
DMLfit = DMLfit.multiFactor(BSobj, design=design, 
                            formula=~Temperature+Inoculation)

# check to see which column refers to which comparison
colnames(DMLfit$X)

# look to see the positions affected by each factor 
DMLtest.temp = DMLtest.multiFactor(DMLfit, coef=2) #ambient v heat-stressed
DMLtest.path = DMLtest.multiFactor(DMLfit, coef=5) #control v pathogen
dml = DMLtest.multiFactor(DMLfit, coef=1)

# now we can call different methylated regions
temp.dmrs = callDMR(DMLtest.temp, p.threshold=0.001) # 96 DMRs
path.dmrs = callDMR(DMLtest.path, p.threshold=0.001) # 8 DMRs
dmrs = callDMR(dml, p.threshold=0.001) # no DMRs from bmc/both

# save to files
write.table(temp.dmrs, "t2_temperature_dmrs.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(path.dmrs, "t2_pathogen_dmrs.txt", sep="\t", row.names=FALSE, quote = FALSE)

# create a Granges object from results to get the sequences associated with the DMRs
library(GenomicRanges)
dmr_gr <- GRanges(seqnames = path.dmrs$chr,
                  ranges = IRanges(start = path.dmrs$start, end = path.dmrs$end), strand = "*")

mcols(dmr_gr)$length <- path.dmrs$length
mcols(dmr_gr)$nCG <- path.dmrs$nCG
mcols(dmr_gr)$areaStat <- path.dmrs$areaStat

library(Rsamtools)
library(Biostrings)

genome <- FaFile("Pver_genome_assembly_v1.0.fasta")
open(genome)
 
dmr_seqs <- getSeq(genome, dmr_gr)

writeXStringSet(dmr_seqs, "path_dmr_sequences.fa")
#

## Plot the DMRs ##

library(bsseq)
library(rtracklayer)
library(GenomicRanges)

# dmr table
dmrs <- read.table("t2_pathogen_dmrs.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# gff file
gff <- import("Pver_genome_assembly_v1.0.gff")

# colors
inoc_colors <- c(bmc = "steelblue", both = "purple", ctrl = "gray50", path = "tomato")

exon_col <- rgb(0.85, 0.92, 1.0, 0.6) # pale blue for exons
intron_col <- rgb(0.95, 0.95, 0.95, 0.6) # light grey for introns

# smoothing on 100bp sliding window
smooth_meth <- function(positions, values, bandwidth = 100) {
  smoothed <- numeric(length(positions))
  for (i in seq_along(positions)) {
    in_window <- abs(positions - positions[i]) <= bandwidth / 2
    smoothed[i] <- mean(values[in_window], na.rm = TRUE)
  }
  smoothed
}

# plot
genes <- unique(dmrs$gene)
n_cols <- 3
n_rows <- 2

pdf("dmrs_by_gene.pdf", width = 18, height = n_rows * 6)
par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 2), oma = c(0, 0, 3, 8))

for (gene_id in genes) {
  gene_dmrs <- dmrs[dmrs$gene == gene_id, ]
  chr <- gene_dmrs$chr[1]
  
  # get gene coordinates from gff
  gene_hits <- gff[!is.na(gff$ID) & grepl(gene_id, gff$ID, fixed = TRUE) &
                     as.character(seqnames(gff)) == chr]
  
  if (length(gene_hits) == 0) {
    gene_start <- min(gene_dmrs$start)
    gene_end <- max(gene_dmrs$end)
    exons <- NULL
  } else {
    gene_start <- min(start(gene_hits))
    gene_end <- max(end(gene_hits))
    
    # pull exon features for gene
    exons <- gene_hits[!is.na(gene_hits$type) & gene_hits$type == "exon"]
  }
  
  # subset to just visualize gene body
  roi <- GRanges(chr, IRanges(gene_start, gene_end))
  bsub <- subsetByOverlaps(BSobj, roi)
  
  if (length(bsub) == 0) {
    plot.new()
    title(main = paste0(gene_id, "\n(no data)"))
    next
  }
  
  meth <- getMeth(bsub, type = "raw", what = "perBase")
  positions <- start(bsub)
  
  # zoom x axis to just the region within the gene body with methylation data
  has_data <- rowSums(!is.na(meth) & meth > 0) > 0
  if (sum(has_data) == 0) has_data <- !is.na(rowMeans(meth, na.rm = TRUE))
  
  data_min <- min(positions[has_data])
  data_max <- max(positions[has_data])
  pad_zoom <- 100 # +/- 100 bp around methylation regions
  xlim_abs <- c(data_min - pad_zoom, data_max + pad_zoom)
  
  # relative positions (to gene start, for readable axis)
  pos_rel <- positions - gene_start
  xlim_rel <- xlim_abs - gene_start
  
  # start with a blank canvas
  plot(NULL,
       xlim = xlim_rel,
       ylim = c(0, 1),
       xlab = paste0("Position in ", gene_id, " (bp)"),
       ylab = "Methylation",
       main = gene_id)
  
  # shade intron background across whole plot window
  rect(xlim_rel[1], 0, xlim_rel[2], 1,
       col = intron_col, border = NA)
  
  # shade exon regions
  if (!is.null(exons) && length(exons) > 0) {
    exon_starts <- start(exons) - gene_start
    exon_ends <- end(exons) - gene_start
    for (k in seq_along(exon_starts)) {
      rect(exon_starts[k], 0, exon_ends[k], 1,
           col = exon_col, border = NA)
    }
  }
  
  # add legend for exon/intron shading
  legend("topleft",
         legend = c("Exon", "Intron"),
         fill = c(exon_col, intron_col),
         border = "gray60",
         bty = "n", cex = 0.6)
  
  # shade each DMR
  for (j in seq_len(nrow(gene_dmrs))) {
    dmr_start_rel <- gene_dmrs$start[j] - gene_start
    dmr_end_rel <- gene_dmrs$end[j] - gene_start
    rect(dmr_start_rel, 0, dmr_end_rel, 1,
         col = rgb(1, 0, 0, 0.15), border = "red", lty = 2, lwd = 1.2)
    text(x = mean(c(dmr_start_rel, dmr_end_rel)),
         y = 0.97,
         labels = gene_dmrs$genetic_region[j],
         cex = 0.55, col = "red", font = 2)
  }
  
  # smoothed group mean + 95% CI
  for (grp in names(inoc_colors)) {
    idx <- which(design$Inoculation == grp)
    grp_meth <- meth[, idx, drop = FALSE]
    
    grp_mean_raw <- rowMeans(grp_meth, na.rm = TRUE)
    grp_se_raw <- apply(grp_meth, 1,
                          function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
    
    grp_mean <- smooth_meth(positions, grp_mean_raw, bandwidth = 100)
    ci_upper <- pmin(smooth_meth(positions, grp_mean_raw + 1.96 * grp_se_raw, bandwidth = 100), 1)
    ci_lower <- pmax(smooth_meth(positions, grp_mean_raw - 1.96 * grp_se_raw, bandwidth = 100), 0)
    
    polygon(c(pos_rel, rev(pos_rel)),
            c(ci_upper, rev(ci_lower)),
            col = adjustcolor(inoc_colors[grp], alpha.f = 0.15),
            border = NA)
    
    lines(pos_rel, grp_mean, col = inoc_colors[grp], lwd = 2)
  }
}

# add legend
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("right",
       legend = names(inoc_colors),
       col = inoc_colors,
       lty = 1, lwd = 2,
       title = "Inoculation",
       bty = "n", cex = 0.9,
       inset = 0.01)

title(main = "Differentially methylated regions due to Vibrio inoculation",
      outer = TRUE, cex.main = 1.4)
