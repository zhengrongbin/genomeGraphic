## draw igv tracks
library(data.table)
library(ggplot2)
library(Gviz)
library(RColorBrewer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#Fig3-C igv track to show the actual peak
txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
grt <- GeneRegionTrack(txdb_hg38, genome="hg38",showId=TRUE, geneSymbol=TRUE, name="UCSC")
z <- mapIds(org.Hs.eg.db, gene(grt), "SYMBOL", "ENTREZID", multiVals = "first")
zz <- sapply(z, is.null)
z[zz] <- gene(grt)[zz]
gr <- ranges(grt)
mcols(gr)$symbol <- z
grt@range <- gr

bwInfo<-read.table("easy_input.txt",header=F,row.names=1,as.is=T)
head(bwInfo)
gloci<-read.table("loci.bed",header=F,as.is=T)
head(gloci)
genefold<-as.numeric("1")

colnames(gloci)<-c("chr","start","end","strand")
chr<-gloci[rownames(gloci),]$chr
gloci$width<-with(gloci,end-start)
startpoint<-gloci[rownames(gloci),]$start-genefold*gloci[rownames(gloci),]$width
endpoint<-gloci[rownames(gloci),]$end+genefold*gloci[rownames(gloci),]$width

tracklist<-list()
itrack <- IdeogramTrack(genome = "hg38", chromosome = chr,outline=T)
tracklist[["itrack"]]<-itrack

scalebar <- GenomeAxisTrack(scale=0.25,col="black",fontcolor="black",name="Scale",labelPos="above",showTitle=TRUE)
tracklist[["scalebar"]]<-scalebar

axisTrack <- GenomeAxisTrack(labelPos="above",col="black",fontcolor="black",name=paste(chr,":",sep=""),exponent=0,showTitle=TRUE)
tracklist[["axisTrack"]]<-axisTrack

colpal<-rep(brewer.pal(12,"Paired"),20)
coldf<-data.frame(col=colpal[1:nrow(bwInfo)],row.names = rownames(bwInfo),stringsAsFactors = F)

for(index in rownames(bwInfo)){
  bgFile<-file.path("./bw/", bwInfo[index,])
  tracklist[[index]]<-DataTrack(range = bgFile,genome="hg38",type="histogram",
                                name=chartr("_","\n",index),
                                col.histogram=COLORS[which(rownames(bwInfo)==index)])
}

tracklist[["grt"]]<-grt

plotTracks(tracklist, from = startpoint, to = endpoint,
           chromosome=chr,background.panel = "white", background.title = "white",
           col.title="black",col.axis="black",
           rot.title=0,cex.title=0.9,margin=38,title.width=1.5,
           collapseTranscripts = "longest")

pdf("loci.pdf",height=8,width=13)
plotTracks(tracklist, from = startpoint, to = endpoint,
           chromosome=chr,background.panel = "white", background.title = "white",
           col.title="black",col.axis="black",
           rot.title=0,cex.title=0.9,margin=38,title.width=1.5,
           collapseTranscripts = "longest")
dev.off()
sessionInfo()