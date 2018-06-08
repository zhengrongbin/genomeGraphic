# genomeGraphic
In terms of biologists need to draw the genome region graph sometimes, here is a simple script for drawing beautiful transcripts.
This is depandant on R packages called <a href='https://bioconductor.org/packages/release/bioc/html/Gviz.html'>Gviz</a> and <a href='https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html'>GenomicRanges</a>.

The main script is written in draw.R.
The input or required file is genome annotation file.
```{r setup, include=FALSE}
hg38 = read.csv('./db/hg38.refGene')
head(hg38, n = 2)
```
