### Goal: to 

## source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade")
## biocLite("KEGGREST")
## biocLite("DESeq2")

library(readr)
library(KEGGREST)
library(tidyverse)
library(DESeq2)
library("biomaRt")
library(tidyr)
library(ggplot2)
library(scales)
library(gplots)
library(dendextend)
library(RColorBrewer)

### KEGG
### listDatabases()
### kegg_human = keggList("hsa")
### head(kegg_human)

### let's grab the kegg pathways in cancer
listDatabases()
keggList("hsa")
pathways = keggLink("hsa", "pathway")
## 05162  Measles
## 05164  Influenza A
## 05161  Hepatitis B
## 05160  Hepatitis C
## 05168  Herpes simplex infection
## 05167  Kaposi's sarcoma-associated herpesvirus infection
## 05169  Epstein-Barr virus infection
## 05165  Human papillomavirus infection
pathways_of_interest_nums =c('05162', '05164','05160','05168')
pathways_of_interest = paste0('path:hsa', pathways_of_interest_nums)

pathway_genes = pathways[names(pathways) %in% pathways_of_interest]
length(pathway_genes)

get_table_as_tibble = function(file){
    res = read.table(file)
    res$gene=rownames(res)
    as.tibble(res)
}

## load up annotations
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

gene_annotations = getBM(
    attributes=c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene',
                 'chromosome_name', 'start_position', 'end_position'),
    mart = ensembl
)

pathway_entrez_ids = as.integer(gsub('hsa:','',pathway_genes))
pathway_gene_annotations = as.tibble(gene_annotations[which(gene_annotations$entrezgene %in% pathway_entrez_ids),])
pathway_gene_annotations
### load files
results_unfiltered = c(
'./results/featureCounts/r_analysis/hits/SepCon120h_SepRes120h.padj.unfiltered.tsv',
'./results/featureCounts/r_analysis/hits/SepCon24h_SepRes24h.padj.unfiltered.tsv'
)
results_shrunk = c(
    './results/featureCounts/r_analysis/hits/SepCon120h_SepRes120h.padj.tsv',
    './results/featureCounts/r_analysis/hits/SepCon24h_SepRes24h.padj.tsv'
)

p_tables = tibble(
    exp1_file=c(results_unfiltered[1], results_shrunk[1]),
    exp2_file=c(results_unfiltered[2], results_shrunk[2]),
    exp1_name = 'SepCon120h_SepRes120h',
    exp2_name = 'SepCon24h_SepRes24h',
    lfc_type = c('unfiltered_shrunk', 'unshrunk')
)
nrow(p_tables)
row_num=1

sort_types = c('sort_none', 'sort_size')
for(row_num in 1:nrow(p_tables)){
    row = p_tables[row_num,]
    res1 = get_table_as_tibble(row$exp1_file)
    res2 = get_table_as_tibble(row$exp2_file)
    res1$experiment = row$exp1_name
    res2$experiment = row$exp2_name
    res = bind_rows(res1, res2)
    res_kcp = res[res$gene %in% pathway_gene_annotations$ensembl_gene_id,]
    dim(res_kcp)
    for(sort_type in sort_types){
        res_table = res_kcp
        if(sort_type == 'sort_size'){
            res_table$gene = factor(res_table$gene, levels=unique(res1$gene[order(res1$log2FoldChange)]))
        }
        g = ggplot(res_table,
               aes(experiment, gene)) +
            geom_tile(aes(fill = log2FoldChange)) + 
            scale_fill_gradient2(midpoint = 1)
        g = g + labs(title=paste("Spearman correlation:", cor(res1$log2FoldChange, res2$log2FoldChange, method='spearman')))
        out_file = file.path('analyses/heatmap',paste(row$exp1_name, 'vs', row$exp2_name, row$lfc_type ,sort_type, 'heatmap.png', sep='.'))
        ggsave(out_file)
    }
}


## plot data

### for later
## lets load up the raw counts and regularize them for the heatmap
cts = read_delim("results/featureCounts/all_sample_counts.csv", " ")
cts_names = gsub('.*PAV', 'PAV', names(cts))
cts_names = gsub('-\\d+.bam', '', cts_names)
names(cts) = cts_names

coldata = read_tsv("data/sample_meta_data.tsv")

coldata$sample_name = paste0("PAV", coldata$"Pavitra Tube number")
cts = cts[, coldata$sample_name]

stopifnot(all(coldata$sample_name == colnames(cts)))
coldata$CellType = ifelse(grepl("Resistant", coldata$Sample), "Resistant", "Control")
coldata$CultureType = ifelse(grepl("culture", coldata$Sample), "Coculture", "Separate")
levels(coldata$"Time point") = c("120h", "24h")
coldata$group = factor(
    paste(
        substr(coldata$CultureType, 1, 3),
        substr(coldata$CellType, 1, 3),
        coldata$"Time point",
        sep = ""
    )
)
coldata$Time.point = coldata$"Time point"
coldata$full_sample_name = gsub(' ', '', paste(coldata$group, coldata$Repeat, sep='.'))

dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~ Time.point + CellType + CultureType
)


gene_lengths = read_tsv("results/featureCounts/gene_lengths.tsv", col_names=FALSE)
gene_lengths = rename(gene_lengths,X1='ensembl_id', X2='gene_length')
rld = rlog(dds)
rld = as_tibble(assay(rld))

rld = bind_cols(gene_lengths, rld)

# filter out genes with low counts
rld = rld[rowSums(counts(dds)) > 1, ]
rld
pathway_rld = inner_join(pathway_gene_annotations, rld, c('ensembl_gene_id' = 'ensembl_id'))
pathway_rld


remove_cols = which(names(pathway_rld) %in% filter(coldata, Sample != 'Control' & Sample != 'Resistant')$sample_name)
sep_culture_pathway_rld = pathway_rld[, -remove_cols]

sep_culture_pathway_rld = sep_culture_pathway_rld[!duplicated(sep_culture_pathway_rld$ensembl_gene_id),]

rld_dist = dist(scale(sep_culture_pathway_rld[,8:15]), method='euclidean')

clustering = hclust( rld_dist, method = "ward.D" )
plot(clustering)
ord <- clustering$order
sep_culture_pathway_rld[ord < 5,] %>% print(width=1000)
sep_culture_pathway_rld$order = ord
sep_culture_pathway_rld

sep_culture_pathway_rld$ensembl_gene_id= factor(
    sep_culture_pathway_rld$ensembl_gene_id,
    levels = sep_culture_pathway_rld$ensembl_gene_id[sep_culture_pathway_rld$order]
)


sep_culture_pathway_rld_melt = gather(sep_culture_pathway_rld, sample_name, rlog, PAV1:PAV19)
sep_culture_pathway_rld_melt = left_join(sep_culture_pathway_rld_melt, coldata)


names(sep_culture_pathway_rld_melt)
ggplot( sep_culture_pathway_rld_melt, aes(full_sample_name, hgnc_symbol) ) +
  geom_tile(aes(fill = rlog)) +
    scale_fill_gradient(low="white", high=muted("red"))

## scaled_iris2 <- iris2 %>% as.matrix %>% scale
sep_culture_pathway_rld_mat = as.matrix(sep_culture_pathway_rld[,grep('PAV', names(sep_culture_pathway_rld))])
some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
names(coldata)
colnames(sep_culture_pathway_rld_mat) = coldata$full_sample_name[match(colnames(sep_culture_pathway_rld_mat), coldata$sample_name)]

hist(sep_culture_pathway_rld_mat)
min(sep_culture_pathway_rld_mat)
out_file = file.path('analyses/heatmap/separate_some-viral-pathways_heatmap.pdf')

pdf(out_file, height=30, width=8)
heatmap.2(sep_culture_pathway_rld_mat,
          labRow=sep_culture_pathway_rld$hgnc_symbol,
          col=brewer.pal(9,"Reds"),
          scale='row',
          trace='none',
          margins=c(10,5),
          xlab='Biological sample',
          ylab='HGNC gene symbol',
          main='Row-scaled expression of genes in\nMeasles, Influenza A, Hepatitis C,\nor Herpes simplex KEGG pathway')
dev.off()
