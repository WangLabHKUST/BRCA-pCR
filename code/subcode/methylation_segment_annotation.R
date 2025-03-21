setwd("/Users/zmoad/OneDrive - HKUST Connect/HKUST/BreastGGH/12.Submission/03.ScienceAdvance/Code/subcode")

require("TxDb.Hsapiens.UCSC.hg38.knownGene")
require("ChIPseeker")

zsegment_qcBed = "/Users/zmoad/OneDrive - HKUST Connect/HKUST/BreastGGH/22.Methylation/03.Block_Ranksum/00.BlockMatrix/zsegment.qc.bed"

wgbstools = fread(zsegment_qcBed)
peakAnno <- annotatePeak(paste(zsegment_qcBed,sep = ""), tssRegion=c(-5000, 5000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb='org.Hs.eg.db')

peakAnnodf = as.data.frame(as.GRanges(peakAnno))
peakAnnodf$key = paste(peakAnnodf$seqnames,peakAnnodf$start-1,peakAnnodf$end,sep = "_")
write.table(peakAnnodf, file = paste("table.3-2.wgbstools_block.anno.txt",sep = ""),col.names = T,row.names = F,sep = "\t",quote = F)

