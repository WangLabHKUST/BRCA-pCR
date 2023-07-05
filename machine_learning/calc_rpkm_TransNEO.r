library(edgeR)
library(data.table)
library("GenomicFeatures")

gtf_txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh37.87.gtf.gz") # Downloaded from http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
gene_list <- genes(gtf_txdb)
gene_list <- as.data.frame(gene_list)


transneo.counts <- data.frame(fread("transneo-diagnosis-RNAseq-rawcounts.tsv.gz",header=T, sep="\t",stringsAsFactors = F),row.names = 1) # Downloaded from https://github.com/cclab-brca/neoadjuvant-therapy-response-predictor/raw/master/data/transneo-diagnosis-RNAseq-rawcounts.tsv.gz

y <- DGEList(transneo.counts, genes = data.frame(Length=gene_list$width))
minCPM      <- 1
minNoToKeep <- 10
keep        <- rowSums(cpm(y)>minCPM)>=minNoToKeep
y <- y[keep , , keep.lib.sizes=FALSE]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y, method = "TMM")
y_rpkm <- rpkm(y)
write.csv(y_rpkm,'./data/transneo.rpkm.csv')