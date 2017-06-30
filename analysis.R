# Title     : Analysis
# Objective : Process the read count data as differential expression
# Created by: winni
# Created on: 6/23/17

output_dir = file.path("results", "featureCounts")

# read in the data
counts = read.table("results/featureCounts/all_sample_counts.csv", sep = " ", row.names=1,header=TRUE)
names(counts) = gsub('\\.\\d+\\.bam', '', gsub('.*PAV', 'PAV', names(counts)))
head(counts)

# read in sample metadata
metadata = read.table("data/sample_meta_data.tsv", sep="\t", header=TRUE)

cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
sample.def = names(counts)
colors = rep(cbPalette, each = 2)
# need to scale because sample 1 has higher coverage
myPca <- prcomp(t(scale(log2(counts + 1))))

# Adding metadata
pcs = as.data.frame(myPca$x)
pcs$Pavitra.Tube.number = as.numeric(gsub("PAV", "", rownames(pcs)))
pcs = merge(metadata, pcs, by="Pavitra.Tube.number")

# first two principle components
library(ggplot2)
g = ggplot(aes(x=PC1, y=PC2, color=Sample, shape=Time.point), data=pcs) + geom_point()
ggsave(file.path(output_dir, 'pca_plot.ggplot.pdf'))

# pdf(file.path(output_dir, 'pca_plot.pdf'))
# plot(myPca$x[, 1], myPca$x[, 2], col = colors, pch = 1)
# legend("topright", sample.def, pch = 1, col = cbPalette)
# dev.off()

# pairwise correlations between samples
nSamples = ncol(counts)
C = cor(log2(counts + 1), method = "pearson")

d = as.dist(1 - C)
h = hclust(d, method = "ward.D2")
dendro = as.dendrogram(h)

pdf(file.path(output_dir, 'correlation_heatmap.pdf'))
heatmap(C, Colv = dendro, Rowv = dendro)
dev.off()
