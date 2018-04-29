# Title     : Analysis
# Objective : Process the read count data as differential expression
# Created by: winni
# Created on: 6/23/17

## source("http://bioconductor.org/biocLite.R")
## biocLite("DESeq2")
## biocLite("goseq")
## biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
## biocLite("org.Hs.eg.db")
## biocLite("KEGGREST")
## biocLite("clusterProfiler")
## install.packages('UpSetR')
## biocLite("biomaRt")
## install.packages('VennDiagram')
library(org.Hs.eg.db)
library(DESeq2)
library(ggplot2)
library(goseq)
library(KEGGREST)
library(clusterProfiler)
library("biomaRt")
library(VennDiagram)
library(ggplot2)


output_dir = file.path("results", "featureCounts", "r_analysis")
go_term_enrichment_output_dir = file.path(output_dir, "go_term_enrichment")
kegg_dir = file.path(output_dir, "kegg_analysis")
dir.create(kegg_dir, showWarnings=FALSE)

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
g = ggplot(aes(x = PC1, y = PC2, color = Sample, shape = Time.point), data = pcs) + geom_point()
ggsave(file.path(output_dir, 'pca_plot.ggplot.pdf'))

# same thing for reg log
pcs.rld = merge(coldata, as.data.frame(myPCA.rld$x), by = "row.names", sort = FALSE)
g = ggplot(aes(x = PC1, y = PC2, color = Sample, shape = Time.point), data = pcs.rld) + geom_point()
ggsave(file.path(output_dir, 'pca_plot.rld.ggplot.pdf'), width = 5, height = 5)

# pairwise correlations between samples
data = list(log2 = log2(cts + 1), rld = assay(rld))
for (dat_name in c("log2", "rld")) {
    dat = data[[dat_name]]
    C = cor(dat, method = "pearson")

    d = as.dist(1 - C)
    h = hclust(d, method = "ward.D2")
    dendro = as.dendrogram(h)

    pdf(file.path(output_dir, paste(dat_name, 'correlation_heatmap.pdf', sep = ".")), height = 6, width = 6)
    heatmap(C, Colv = dendro, Rowv = dendro)
    dev.off()
}

# differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# order results by p-value
resOrdered <- res[order(res$padj),]

# Nothing to be found
# let's try interaction effects

design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)

# retrieve gene annotations from biomart 
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# head(listAttributes(ensembl))
gene_annotations = getBM(
    attributes=c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene',
                 'chromosome_name', 'start_position', 'end_position'),
    mart = ensembl
)

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
question_name = 'q2'
for (question_name in names(questions)) {
    print(question_name)
    design = questions[[question_name]]$design
    contrasts=c("group", design[1], design[2])

    res_unshrunk = results(dds, contrast = contrasts)

    res_unfiltered = results(dds, contrast = contrasts, independentFiltering=FALSE)
    
    res = lfcShrink(dds, contrast = contrasts, res = res_unshrunk)

    questions[[question_name]]$res = res
    questions[[question_name]]$res_unshrunk = res_unshrunk
    questions[[question_name]]$res_unfiltered = res_unfiltered
}

for (question_name in names(questions)) {
    print(question_name)
    design = questions[[question_name]]$design
    res = questions[[question_name]]$res
    res_unfiltered =  questions[[question_name]]$res_unfiltered

    print(head(res))
    write.table(
        res,
        file.path(output_dir, "hits", paste0(design[1], "_", design[2], ".padj.tsv")),
        sep = "\t"
    )
    write.table(
        res_unfiltered,
        file.path(output_dir, "hits", paste0(design[1], "_", design[2], ".padj.unfiltered.tsv")),
        sep = "\t"
    )
    res_fdr_pc10 = res[(! is.na(res$padj)) & res$padj < 0.1,]
    write.table(
        res_fdr_pc10,
        file.path(output_dir, "hits", paste0(design[1], "_", design[2], ".fdr_10pc.tsv")),
        sep = "\t"
    )
    questions[[question_name]]$res_fdr_pc10 = res_fdr_pc10

    unshrunk = questions[[question_name]]$res_unshrunk
    unshrunk_fdr_pc10 = unshrunk[(! is.na(unshrunk$padj)) & unshrunk$padj < 0.1,]
    write.table(
        unshrunk_fdr_pc10,
        file.path(output_dir, "hits", paste0(design[1], "_", design[2], ".unshrunk.fdr_10pc.tsv")),
        sep = "\t"
    )
}

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

# GO analysis
gene_lengths = read.table("results/featureCounts/gene_lengths.tsv", header = FALSE, row.names = 1)
tmp = gene_lengths[, 1]
names(tmp) = rownames(gene_lengths)
all_gene_lengths = tmp
convert_to_enrichment_map_gen_enrich_res_format = function(goseq_result){
    over_represented_pvalue_is_smaller = goseq_result$over_represented_pvalue < goseq_result$under_represented_pvalue
    phenotype = ifelse(over_represented_pvalue_is_smaller, '+1', '-1')
    pval = ifelse(over_represented_pvalue_is_smaller,
    goseq_result$over_represented_pvalue,
    goseq_result$under_represented_pvalue)

    enrich_df = data.frame(
    GO.ID = goseq_result$category,
    Description = c(""),
    p.val = pval,
    FDR = c(""),
    Phenotype = phenotype
    )
    return(enrich_df)
}
for(q in questions){
    print(q$design)
    if(nrow(q$res_fdr_pc10) < 2){
        next
    }
    genes = as.integer(rownames(q$res) %in% rownames(q$res_fdr_pc10))
    names(genes) = rownames(q$res)
    gene_lengths = all_gene_lengths
    gene_lengths = gene_lengths[names(all_gene_lengths) %in% names(genes)]
    gene_lengths[match(names(gene_lengths), names(genes))]
    stopifnot(all(match(names(gene_lengths), names(genes)) == match(names(genes), names(gene_lengths))))
    pwf = nullp(genes, bias.data = gene_lengths)
    GO.wall = goseq(pwf, "hg38", "ensGene")
    kegg.wall = goseq(pwf, "hg38", "ensGene", test.cats = c("KEGG"))
    em_enrichment = convert_to_enrichment_map_gen_enrich_res_format(GO.wall)

    write.table(
    em_enrichment,
    file = file.path(go_term_enrichment_output_dir, paste(paste(q$design, collapse = "-"), "de", "go-term-enrichment", "tsv", sep = ".")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
    )

    em_kegg_enrichment = convert_to_enrichment_map_gen_enrich_res_format(kegg.wall)
    em_kegg_enrichment$kegg.id = sub("^", "hsa", as.character(em_kegg_enrichment$GO.ID))

    # let's try complete clusterProfiler KEGG analysis
    eg = select(org.Hs.eg.db, keys=rownames(q$res_fdr_pc10), columns=c("ENTREZID"), keytype="ENSEMBL")
    eg = eg[!duplicated(eg$ENTREZID), ]

    kk <- enrichKEGG(gene = eg$ENTREZID,
    organism = 'hsa',
    pvalueCutoff = 0.005)
    geneList = q$res_fdr_pc10$stat
    names(geneList) = rownames(q$res_fdr_pc10)
    geneList = geneList[names(geneList) %in% eg$ENSEMBL]
    names(geneList) = eg$ENTREZID[match(names(geneList), eg$ENSEMBL)]

    design_str = paste(q$design, collapse="-")

    pdf(file.path(kegg_dir, paste0(design_str, ".cnetplot.pdf")), width = 30, height = 30)
    cnetplot(kk, showCategory=8, categorySize="pvalue", foldChange=geneList)
    dev.off()

    ## pdf(file.path(kegg_dir, paste0(design_str,".upsetplot.pdf")))
    ## upsetplot(kk)
    ## dev.off()

}
