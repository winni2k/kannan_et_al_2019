# Title     : Analysis
# Objective : Process the read count data as differential expression
# Created by: winni
# Created on: 6/23/17

source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)
library(ggplot2)

output_dir = file.path("results", "featureCounts", "r_analysis")

# read in the data
cts = read.table("results/featureCounts/all_sample_counts.csv", sep = " ", row.names = 1, header = TRUE)
names(cts) = gsub('\\.\\d+\\.bam', '', gsub('.*PAV', 'PAV', names(cts)))
head(cts)

# read in sample coldata
coldata = read.table("data/sample_meta_data.tsv", sep = "\t", header = TRUE)
rownames(coldata) = paste0("PAV", coldata$Pavitra.Tube.number)

# reorder counts
cts <- cts[, rownames(coldata)]
stopifnot(all(rownames(coldata) == colnames(cts)))
coldata$CellType = ifelse(grepl("Resistant", coldata$Sample), "Resistant", "Control")
coldata$CultureType = ifelse(grepl("culture", coldata$Sample), "Coculture", "Separate")
levels(coldata$Time.point) = c("120h", "24h")
coldata$group = factor(
paste(
substr(coldata$CultureType, 1, 3),
substr(coldata$CellType, 1, 3),
coldata$Time.point,
sep = ""
)
)

dds <- DESeqDataSetFromMatrix(
countData = cts,
colData = coldata,
design = ~ Time.point + CellType + CultureType
)

# Remove rows that have almost no reads
dds <- dds[rowSums(counts(dds)) > 1,]
# re-level the factor data
dds$CultureType = relevel(dds$CultureType, ref = "Separate")


cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
sample.def = names(cts)
colors = rep(cbPalette, each = 2)

# need to scale because sample 1 has higher coverage
myPca <- prcomp(t(scale(log2(cts + 1))))

# also try scaling with regularized transform
rld = rlog(dds)
myPCA.rld = prcomp(t(assay(rld)))

# Adding coldata
pcs = as.data.frame(myPca$x)
# pcs$Pavitra.Tube.number = as.numeric(gsub("PAV", "", rownames(pcs)))
pcs = merge(coldata, pcs, by = "row.names", sort = FALSE)

# first two principle components
library(ggplot2)
g = ggplot(aes(x = PC1, y = PC2, color = Sample, shape = Time.point), data = pcs) + geom_point()
ggsave(file.path(output_dir, 'pca_plot.ggplot.pdf'))

# same thing for reg log
pcs.rld = merge(coldata, as.data.frame(myPCA.rld$x), by = "row.names", sort = FALSE)
g = ggplot(aes(x = PC1, y = PC2, color = Sample, shape = Time.point), data = pcs.rld) + geom_point()
ggsave(file.path(output_dir, 'pca_plot.rld.ggplot.pdf'), width=5, height=5)

# pairwise correlations between samples
data = list(log2 = log2(cts + 1), rld = assay(rld))
for (dat_name in c("log2", "rld")) {
    dat = data[[dat_name]]
    C = cor(dat, method = "pearson")

    d = as.dist(1 - C)
    h = hclust(d, method = "ward.D2")
    dendro = as.dendrogram(h)

    pdf(file.path(output_dir, paste(dat_name, 'correlation_heatmap.pdf', sep = ".")), height=6, width=6)
    heatmap(C, Colv = dendro, Rowv = dendro)
    dev.off()
}

# differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res

# order results by p-value
resOrdered <- res[order(res$padj),]

# Nothing to be found
# let's try interaction effects

design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)

# 1. What genes are differentially expressed between control (separate) and resistant (separate)
# cells at 24 hours?
# 2. What genes are differentially expressed between control (separate) and control (co-culture) at 24 hours?
# 3. What genes are differentially expressed between control (separate) and resistant (separate) cells at 120 hrs?
# 4. What genes are differentially expressed between control (separate) and control (co-culture) at 120 hrs?

questions = list(
q1 = list(design = c("SepCon24h", "SepRes24h")),
q2 = list(design = c("SepCon24h", "CocCon24h")),
q3 = list(design = c("SepCon120h", "SepRes120h")),
q4 = list(design = c("SepCon120h", "CocCon120h")),
q6a = list(design = c("CocCon24h", "SepRes24h")),
q6b = list(design = c("CocCon120h", "SepRes120h"))
)
for (question_name in names(questions)) {
    print(question_name)
    design = questions[[question_name]]$design
    res = results(dds, contrast = c("group", design[1], design[2]))
    res = lfcShrink(dds, coef = 2, res = res)
    questions[[question_name]]$res = res
}

for (question_name in names(questions)) {
    print(question_name)
    design = questions[[question_name]]$design
    res = questions[[question_name]]$res
    res_fdr_pc10 = res[(! is.na(res$padj)) & res$padj < 0.1,]
    write.table(
    res_fdr_pc10,
    file.path(output_dir, "hits", paste0(design[1], "_", design[2], ".fdr_10pc.tsv")),
    sep = "\t"
    )
    questions[[question_name]]$res_fdr_pc10 = res_fdr_pc10
}

library(VennDiagram)
gene_lists = lapply(lapply(questions, "[[", "res_fdr_pc10"), rownames)
v1 <- venn.diagram(gene_lists[1 : 4], filename = NULL, fill = rainbow(4))
pdf(file.path(output_dir, "q1-q4.venn.pdf"), width = 5, height = 5)
grid.newpage()
grid.draw(v1)
dev.off()

names(gene_lists) = sapply(lapply(questions, "[[", "design"), paste, collapse = "--")
v1 <- venn.diagram(gene_lists[1 : 4], filename = NULL, fill = rainbow(4))
pdf(file.path(output_dir, "q1-q4.venn.+names.pdf"), width = 12, height = 12)
grid.newpage()
grid.draw(v1)
dev.off()



