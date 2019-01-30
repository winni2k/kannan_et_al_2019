### Created on: 2017-10-17
### Created by: winni!

### source("https://bioconductor.org/biocLite.R")
### biocLite("biomaRt")
### install.packages("ggrepel")


library("biomaRt")
library(ggplot2)
library(ggrepel)
q4_tsv = '../results/featureCounts/r_analysis/hits/SepCon120h_CocCon120h.padj.tsv'
q4 = read.csv(q4_tsv, sep="\t")
q4$ensembl_gene_id = rownames(q4)
head(q4)
dim(q4)

## Only keep genes with non missing padj
q4 = q4[complete.cases(q4),]

## Let's retrieve gene symbols for ensemble ids
## listMarts()
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
attributes
attributes$name[grep('symbol',attributes$name)]
gene_symbols = getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),                     
                     mart=ensembl)

nrow(q4)
nrow(gene_symbols)
head(gene_symbols)

q4 = merge(q4, gene_symbols, by='ensembl_gene_id')
names(q4)
nrow(q4)
tail(q4[order(q4$padj),])

q4$Significant <- ifelse(q4$padj < 0.05, "FDR < 0.1", "Not Sig")

g = ggplot(q4, aes(x = log2FoldChange, y = -log10(pvalue)))
g = g + geom_point(aes(color = Significant), size=0.5, alpha=0.7) 
g = g + scale_color_manual(values = c("red", "grey"))
g = g + theme_bw(base_size = 12) + theme(legend.position = "bottom")
g = g + geom_text_repel(
            data = subset(q4, padj < 0.05),
            aes(label = hgnc_symbol),
            size = 5,
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines")
        )
g


### Let's take a look at the count data: What is OAS3 expression level?
# read in the data
cts = read.table("../results/featureCounts/all_sample_counts.csv", sep = " ", row.names = 1, header = TRUE)
names(cts) = gsub('\\.\\d+\\.bam', '', gsub('.*PAV', 'PAV', names(cts)))
head(cts)

# read in sample coldata
coldata = read.table("../data/sample_meta_data.tsv", sep = "\t", header = TRUE)
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
head(coldata)

coldata$group == 'SepCon120h'
hgnc_gene = 'OAS3'
oas3_ensembl_gene_id = gene_symbols[gene_symbols$hgnc_symbol == hgnc_gene,]$ensembl_gene_id
oas3_ensembl_gene_id
cts[rownames(cts) == oas3_ensembl_gene_id, coldata$group == 'SepCon120h']
cts[rownames(cts) == oas3_ensembl_gene_id, coldata$group == 'CocCon120h']
q4[q4$hgnc_symbol == hgnc_gene,]

MANSC4 = gene_symbols[gene_symbols$hgnc_symbol == 'MANSC4',]$ensembl_gene_id
cts[rownames(cts) == MANSC4,]
