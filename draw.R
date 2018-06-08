# --- Author: Rongbin Zheng ---
# --- Bioinformatics & Epigenetics, Tongji University ---

# install packages:
# source("https://bioconductor.org/biocLite.R")
# biocLite("Gviz", 'GenomicRanges')
library(Gviz)
library(GenomicRanges)

hg38 = read.csv('./db/hg38.refGene', header = F, sep = '\t')

## plot transcript in separate figure but keep gene range 
format1 = function(gene, mat){
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

## plot gene transcripts in one figure
format2 = function(gene, mat){
  ## gene = offical name, mat is the sub matrix of hg38 object
  mat = as.matrix(mat)
  minSite = min(as.numeric(do.call(c, strsplit(mat[,5], '\\,'))))
  maxSite = max(as.numeric(do.call(c, strsplit(mat[,6], '\\,'))))
  df = data.frame()
  for (i in 1:nrow(mat)){
    transcripts = as.character(mat[i,2])
    chromosome = as.character(mat[i,3])
    strand = as.character(mat[,4])
    exons = t(do.call(rbind, strsplit(mat[i,10:11], '\\,')))
    transcripts = rep(transcripts, nrow(exons))
    chromosome = rep(chromosome, nrow(exons))
    strand = rep(strand, nrow(exons))
    df = rbind(df, data.frame('chr' = chromosome, 'start' = as.numeric(exons[,1]),
                    'end' = as.numeric(exons[,2]), 'strand' = strand, 'transcript' = transcripts))
  }
  chr = chromosome[1]
  gen = c(chr='hg38') 
  gtrack <- GenomeAxisTrack(name="Axis",
                            range <- IRanges(start=df$start, end=df$end))
  grtrack <- GeneRegionTrack(df, genome = gen,
                           chromosome = chr, name = gene,
                           strand = df$strand,
                           transcriptAnnotation = "transcript",
                           background.panel = "#FFFEDB",
                           background.title = 'darkblue')
  png(paste0('./result/hg38_', gene, '.png'), res = 300, height = 280*nrow(mat), width = 1800)
  plotTracks(list(itrackList[[chr]], gtrack, grtrack), col = 'red',
             margin = 10, innerMargin = 10, fontcolor = 'black', fontface = 'bold',
             stackHeight = 0.1, extend.left = 50)
  dev.off()
}

chrs = paste0('chr', c(1:22, 'X', 'Y'))
hg38 = subset(hg38, V3 %in% chrs)
for (gene in unique(as.vector(hg38$V13))[1:5]){
  print(gene)
  mat = subset(hg38, V13 == gene)
  if (nrow(mat) > 1){
    try(format1(gene, mat))
  }
}









