suppressPackageStartupMessages(library(MEDIPS))
library(BSgenome.Bdistachyon.Ensembl.Bdi1)
library(yaml)

BSgenome <- 'BSgenome.Bdistachyon.Ensembl.Bdi1'
uniq <- 1e-3
extend <- 300
shift <- 0
ws <- 100

output <- commandArgs(T)[1]
config <- yaml.load_file('config.yaml')

## createSet ##
bam_files <- paste0('bam/', config$samples, '.bam')

MeDIP_Set <- list()
for (i in 1:length(bam_files)) {
  MeDIP_Set[[i]] <- MEDIPS.createSet(file=bam_files[i], BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, window_size=ws)
}
names(MeDIP_Set) <- config$samples
CS <- MEDIPS.couplingVector(pattern = 'CG', refObj = MeDIP_Set[[1]])

VS <- t(combn(config$treat, 2))
DMR_sig <- list()
for (i in 1:nrow(VS)) {
  cat(VS[i, 1], 'VS', VS[i, 2])
  mr.edgeR <- MEDIPS.meth(MSet1=MeDIP_Set[[VS[i, 1]]], MSet2=MeDIP_Set[[VS[i, 2]]], CSet=CS, ISet1=MeDIP_Set$Input, ISet2=MeDIP_Set$Input, p.adj='fdr', diff.method='edgeR', MeDIP=T, CNV=F, minRowSum = 10)
  DMR_sig[[i]] <- MEDIPS.selectSig(results=mr.edgeR, p.value=0.1, adj=T, ratio=NULL, bg.counts=NULL, CNV=F)
}
names(DMR_sig) <- paste0(VS[,1], "_vs_", VS[,2])

save(list = c('MeDIP_Set', 'DMR_sig'), file = output)

