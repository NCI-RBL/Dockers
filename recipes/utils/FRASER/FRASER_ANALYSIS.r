#!/usr/bin/env Rscript

library(FRASER)

sampleTable <- fread("/data2/test.table.tsv")
bamFiles <- sampleTable[,bamFile]

settings <- FraserDataSet(colData=sampleTable, workingDir="/data2/FRASER_output")
fds <- countRNAData(settings)
fds <- calculatePSIValues(fds)
fds <- filterExpressionAndVariability(fds, minExpressionInOneSample=20,
minDeltaPsi=0.0, filter=TRUE)
fds <- FRASER(fds, q=c(psi5=2, psi3=3, theta=3))

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgDb <- org.Hs.eg.db
fds <- annotateRangesWithTxDb(fds, txdb=txdb, orgDb=orgDb)
res <- results(fds, zScoreCutoff=NA, padjCutoff=NA, deltaPsiCutoff=NA)
write.table(res, '/data2/FRASER_RES')
plotVolcano(fds, sampleID="sample1", type="psi5", aggregate=TRUE)
