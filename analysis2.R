# Title     : TODO
# Objective : TODO
# Created by: winni
# Created on: 9/28/17
# Creating an intersection of SepCon120h-CocCon120h vs SepCon120h-SepRes120h

library(ggplot2)

q4_tsv = 'results/featureCounts/r_analysis/hits/SepCon120h_CocCon120h.fdr_10pc.tsv.gz'
q3_tsv = 'results/featureCounts/r_analysis/hits/SepCon120h_SepRes120h.fdr_10pc.tsv.gz'
output = 'results/featureCounts/r_analysis/hits/analysis_2_isec_genes.tsv'

q4 = read.csv(gzfile(q4_tsv), sep="\t")
q3 = read.csv(gzfile(q3_tsv), sep="\t")

isec_genes = rownames(q3[rownames(q3) %in% rownames(q4),])

comp_table = data.frame(gene_name = isec_genes, stat_q3=q3[isec_genes,]$stat, stat_q4=q4[isec_genes,]$stat)

ggplot(aes(x=stat_q3, y=stat_q4), data=comp_table) + geom_point()

isec_genes = data.frame(q3_q4_isec_genes = isec_genes)
write.table(isec_genes, output, sep="\t", row.names=FALSE)
