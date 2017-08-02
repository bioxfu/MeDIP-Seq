suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(ChIPseeker))
# setting the java.parameters option before the rJava package is loaded.
options(java.parameters = "-Xmx8000m")
library(xlsx)

argv <- commandArgs(T)
input_rdata <- argv[1]
sqlite <- argv[2]
gene2GO <- argv[3]
output_rdata <- argv[4]
output_xlsx <- argv[5]

load(input_rdata)
txdb <- loadDb(sqlite)
GO <- read.table(gene2GO, sep='\t', quote="", stringsAsFactors = F)
colnames(GO) <- c('geneId', 'GOterm')

gene_anno <- function(x) {
  peakAnno <- annotatePeak(x, TxDb=txdb)
  anno <- as.data.frame(peakAnno)
  anno <- merge(anno, GO, all.x=T)
  return(anno)
}
DMR_sig_Grange <- lapply(DMR_sig, makeGRangesFromDataFrame, keep.extra.columns=TRUE)
DMR_sig_anno <- lapply(DMR_sig_Grange, gene_anno)
save(list = c('MeDIP_Set', 'DMR_sig_anno'), file = output_rdata)

wb <- createWorkbook()
for (i in 1:length(DMR_sig_anno)) {
  sheet <- createSheet(wb, sheetName=names(DMR_sig_anno)[i])
  addDataFrame(DMR_sig_anno[[i]], sheet, row.names=FALSE)
}
saveWorkbook(wb, output_xlsx)
