# Install packages for differential expression analysis

options(repos=structure(c(CRAN="https://ftp.acc.umu.se/mirror/CRAN/")))

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")

BiocManager::install(
                 c("DESeq2", "goseq", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "KEGGREST", "clusterProfiler", "biomaRt"),
                 version = "3.8"
             )

install.packages('UpSetR')
install.packages('VennDiagram')





