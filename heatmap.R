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
### KEGG
### listDatabases()
### kegg_human = keggList("hsa")
### head(kegg_human)

### let's grab the kegg pathways in cancer
listDatabases()
keggList("hsa")
pathways = keggLink("hsa", "pathway")
kegg_cancer_pathway = "path:hsa05200"
kegg_cancer_pathway_genes = pathways[names(pathways) == kegg_cancer_pathway]

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

kegg_cancer_pathway_entrez_ids = as.integer(gsub('hsa:','',kegg_cancer_pathway_genes))
kcp_gene_annotations = gene_annotations[which(gene_annotations$entrezgene %in% kegg_cancer_pathway_entrez_ids),]


### load files

results_unfiltered = c(
'./results/featureCounts/r_analysis/hits/SepCon120h_SepRes120h.padj.unfiltered.tsv',
'./results/featureCounts/r_analysis/hits/SepCon24h_SepRes24h.padj.unfiltered.tsv',
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
    row$exp1_name
    res1$experiment = row$exp1_name
    res2$experiment = row$exp2_name
    res = bind_rows(res1, res2)
    res_kcp = res[res$gene %in% kcp_gene_annotations$ensembl_gene_id,]
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
        out_file = file.path('analyses/heatmap',paste(exp1, 'vs',exp2, row$lfc_type ,sort_type, 'heatmap.png', sep='.'))
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
coldata_df = coldata %>% as.data.frame
names(coldata_df)

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

coldata$


