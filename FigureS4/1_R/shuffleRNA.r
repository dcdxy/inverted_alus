library("regioneR")
library("glue")

i = commandArgs(trailingOnly = TRUE)

circRNA <- read.table("human_circrna_sorted.bed", header=FALSE, sep="\t",stringsAsFactors=FALSE)
mask <- read.table("intergenic.bed", header=FALSE, sep="\t",stringsAsFactors=FALSE)

result <- randomizeRegions(circRNA, mask=mask, genome="hg38", per.chromosome=TRUE)

df <- data.frame(chr=seqnames(result), start=start(result)-1, end=end(result))
write.table(format(df, scientific=FALSE), file=glue("temp_{i}.bed"), quote=F, sep="\t", row.names=F, col.names=F)

