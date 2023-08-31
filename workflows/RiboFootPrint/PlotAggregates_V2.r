# PlotAggregates_V2.r

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

DATA <- read.table(args[1])
DATA <- data.frame(DATA)


OUTDIR <- args[2]
OUTNAME <- args[3]
NORM_FACTOR <- read.table(args[4])[1:1]

utr5 <- DATA[which(DATA$V1 == '5utr'), ]
cds  <- DATA[which(DATA$V1 == 'cds' ), ]
utr3 <- DATA[which(DATA$V1 == '3utr'), ]

utr5$transcript_pos <- utr5$V2 + 0
cds$transcript_pos  <- cds$V2  + 1 + 0.0001
utr3$transcript_pos <- utr3$V2 + 2 + 0.0002

#utr5$V3 <- utr5$V3 / (19006 - 958)
#cds$V3 <- cds$V3 / (19006)
#utr3$V3 <- utr3$V3 / (19006 - 648)

DATA <- rbind.data.frame(utr5,cds,utr3)

DATA <- DATA[order( DATA[,1], DATA[,4] ),]


#write.csv(DATA,file=paste0(OUTDIR,'test_aggregate.csv'), row.names=FALSE)

#g=ggplot() +
#        geom_point(data = utr5, 
#                  aes(x=V2, y= V3 / NORM_FACTOR$V1 * 1e6),
#                  size = 0.05) +

#        theme_classic(base_size = 20) +
#        labs(y = 'RPM',
#             x = '5utr Position',
#             title = OUTNAME) +
#        guides(color = guide_legend(title = "")) +
#        coord_cartesian(ylim = c(0,2e3))

#ggsave(g, file=paste0(OUTDIR,OUTNAME,".5utr.aggregates.png"), width = 8, height = 5)

#g=ggplot() +
#        geom_point(data = cds, 
#                  aes(x=V2, y= V3 / NORM_FACTOR$V1 * 1e6),
#                  size = 0.05) +

#        theme_classic(base_size = 20) +
#        labs(y = 'RPM',
#             x = 'cds Position',
#             title = OUTNAME) +
#        guides(color = guide_legend(title = "")) +
#        coord_cartesian(ylim = c(0,2e3))


#ggsave(g, file=paste0(OUTDIR,OUTNAME,".cds.aggregates.png"), width = 8, height = 5)

#g=ggplot() +
#        geom_point(data = utr3, 
#                  aes(x=V2, y= V3 / NORM_FACTOR$V1 * 1e6),
#                  size = 0.05) +

#        theme_classic(base_size = 20) +
#        labs(y = 'RPM',
#             x = '3utr Position',
#             title = OUTNAME) +
#        guides(color = guide_legend(title = "")) +
#        coord_cartesian(ylim = c(0,2e3))


#ggsave(g, file=paste0(OUTDIR,OUTNAME,".3utr.aggregates.png"), width = 8, height = 5)

g=ggplot() +
        geom_point(data = DATA,
                  aes(x=transcript_pos, y= V3 / NORM_FACTOR$V1 * 1e6),
                  size = 0.05) +

        theme_classic(base_size = 20) +
        labs(y = 'RPM',
             x = 'Transcript Position',
             title = OUTNAME) +
        guides(color = guide_legend(title = "")) +
        coord_cartesian(ylim = c(0,4e3))#, xlim = c(0.99, 1.02))

ggsave(g, file=paste0(OUTDIR,OUTNAME,".transcript.aggregates.png"), width = 8, height = 5)
