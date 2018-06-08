# --- Author: Rongbin Zheng ---
# install packages:
# source("https://bioconductor.org/biocLite.R")
# biocLite("Gviz", 'GenomicRanges')
library(Gviz)
library(GenomicRanges)

hg38 = read.csv('./db/hg38.refGene', header = F, sep = '\t')

format = function(gene, mat){
  ## gene = offical name, mat is the sub matrix of hg38 object
  mat = as.matrix(mat)
  minSite = min(as.numeric(do.call(c, strsplit(mat[,5], '\\,'))))
  maxSite = max(as.numeric(do.call(c, strsplit(mat[,6], '\\,'))))
  for (i in 1:nrow(mat)){
    x = mat[i,]
    t = paste(x[3], x[5], x[6], x[2], x[13], sep = ':')
    t = gsub(' ', '', t)
    transcripts = as.character(mat[i,2])
    chromosome = as.character(mat[i,3])
    strand = as.character(mat[i,4])
    exons = t(do.call(rbind, strsplit(mat[i,10:11], '\\,')))
    transcripts = rep(transcripts, nrow(exons))
    chromosome = rep(chromosome, nrow(exons))
    strand = rep(strand, nrow(exons))
    subdf = data.frame('chr' = chromosome, 'start' = as.numeric(exons[,1]),
                       'end' = as.numeric(exons[,2]), 'strand' = strand, 'transcript' = transcripts)
    chr = chromosome[1]
    gen = c(chr='hg38')
    grtrack <- GeneRegionTrack(subdf, genome = gen, 
                               chromosome = chr,
                               strand = subdf$strand, 
                               background.panel = "#FFFEDB",
                               background.title = 'darkblue',
                               fill = 'blue')
    png(paste0('./result/', t, '.png'), res = 100, height = 20, width = 500)
    plotTracks(list(grtrack), col = 'blue', panel.only=T, from = minSite, to = maxSite)
    dev.off()
  } 
}

chrs = paste0('chr', c(1:22, 'X', 'Y'))
hg38 = subset(hg38, V3 %in% chrs)
for (gene in unique(as.vector(hg38$V13))[1:5]){
  print(gene)
  mat = subset(hg38, V13 == gene)
  if (nrow(mat) > 1){
    try(format(gene, mat))
  }
}
