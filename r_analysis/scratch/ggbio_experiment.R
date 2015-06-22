# Experiment with ggbio

# install from bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite("ggbio")
install.packages("ggbio")

library(ggbio)

# 1 ideogram and track
p.ideo <- Ideogram(genome = "hg19")
p.ideo

library(GenomicRanges)
p.ideo + xlim(GRanges("chr2", IRanges(1e8, 1e8+10000000)))

# 2 

biocLite("Homo.sapiens")
library(Homo.sapiens)
data("CRC", package = "biovizBase")
head(hg19sub)
autoplot(hg19sub, layout = "circle", fill = "gray70")

# circular plot with ordinate
p <- ggbio() + circle(hg19sub, geom = "ideo", fill = "gray70") + circle(hg19sub, geom = "scale", size = 2) +
  circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
p

# circular wider radius
p <- ggbio(trackWidth = 10, buffer = 0, radius = 10) 
p <- p + circle(hg19sub, geom = "ideo", fill = "grey70")
p <- p + circle(hg19sub, geom = "scale", size = 2) 
p <- p + circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
p

head(mut.gr)

p <- ggbio() + circle(mut.gr, geom = "rect", color = "steelblue") +
  circle(hg19sub, geom = "ideo", fill = "gray70") +
  circle(hg19sub, geom = "scale", size = 2) + circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
p

# add some links
head(crc.gr)
gr.crc1 <- crc.gr[values(crc.gr)$individual == "CRC-1"]
gr.crc1

p <- p + circle(gr.crc1, geom = "point", aes(y = score, size = tumreads), color = "red", grid = TRUE, radius = 30) + scale_size(range = c(1, 2.5))
p


#
p <- p + circle(gr.crc1, geom = "link", linked.to = "to.gr", aes(color = rearrangements), radius = 23)
p

# do the whole thing
p <- ggbio() + circle(gr.crc1, geom = "link", linked.to = "to.gr", aes(color = rearrangements)) +
  circle(gr.crc1, geom = "point", aes(y = score, size = tumreads), color = "red", grid = TRUE) + scale_size(range = c(1, 2.5)) +
  circle(mut.gr, geom = "rect", color = "steelblue") + circle(hg19sub, geom = "ideo", fill = "gray70") + circle(hg19sub, geom = "scale", size = 2) +
  circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
p

snp <- read.table(system.file("extdata", "plink.assoc.sub.txt", package = "biovizBase"), header = TRUE)
require(biovizBase)
gr.snp <- transformDfToGr(snp, seqnames = "CHR", start = "BP", width = 1)
require(GenomicRanges)
gr.snp <- keepSeqlevels(gr.snp, as.character(1:22))
seqlengths(gr.snp)
data(ideoCyto, package = "biovizBase")
seqlengths(gr.snp) <- as.numeric(seqlengths(ideoCyto$hg18)[1:22])
gr.snp <- gr.snp[!is.na(gr.snp$P)]
values(gr.snp)$pvalue <- -log10(values(gr.snp)$P)
head(gr.snp)
autoplot(gr.snp, geom = "point", coord = "genome", aes(y = pvalue))
plotGrandLinear(gr.snp, aes(y = pvalue), color = c("#7fc97f", "#fdc086"))

plotGrandLinear(gr.snp, aes(y = pvalue), color = c("#7fc97f", "#fdc086"), cutoff = 3, cutoff.color = "blue", cutoff.size = 0.2)

plotGrandLinear(gr.snp, aes(y = pvalue, color = OR), spaceline = TRUE, legend = TRUE)



# Stacked karyogram
data(ideoCyto, package = "biovizBase")
autoplot(seqinfo(ideoCyto$hg19), layout = "karyogram")
biovizBase::isIdeogram(ideoCyto$hg19)
autoplot(ideoCyto$hg19, layout = "karyogram", cytoband = TRUE)


# add data to karyogram
data(darned_hg19_subset500, package = "biovizBase") 
dn <- darned_hg19_subset500 
library(GenomicRanges) 
seqlengths(dn)
## add seqlengths ## we have seqlegnths information in another data set 
seqlengths(dn) <- seqlengths(ideoCyto$hg19)[names(seqlengths(dn))] 
## then we change order 
dn <- keepSeqlevels(dn, paste0("chr", c(1:22, "X"))) 
seqlengths(dn) 
autoplot(dn, layout = "karyogram")

## since default is geom rectangle, even though it's looks like segment 
## we still use both fill/color to map colors 
autoplot(dn, layout = "karyogram", aes(color = exReg, fill = exReg))

## since default is geom rectangle, even though it's looks like segment 
## we still use both fill/color to map colors 
autoplot(dn, layout = "karyogram", aes(color = exReg, fill = exReg), alpha = 0.5) +
scale_color_discrete(na.value = "brown")



## 5

# Link ranges to your data
library(TxDb.Hsapiens.UCSC.hg19.knownGene) 
library(ggbio) 
data(genesymbol, package = "biovizBase") 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
model <- exonsBy(txdb, by = "tx")
model17 <- subsetByOverlaps(model, genesymbol["RBM17"]) 
exons <- exons(txdb) 
exon17 <- subsetByOverlaps(exons, genesymbol["RBM17"]) 
## reduce to make sure there is no overlap
## just for example
exon.new <- reduce(exon17) 
## suppose 
values(exon.new)$sample1 <- rnorm(length(exon.new), 10, 3) 
values(exon.new)$sample2 <- rnorm(length(exon.new), 10, 10) 
values(exon.new)$score <- rnorm(length(exon.new)) 
values(exon.new)$significant <- sample(c(TRUE,FALSE), size = length(exon.new),replace = TRUE) 
## data ready 
exon.new

p17 <- autoplot(txdb, genesymbol["RBM17"])
plotRangesLinkedToData(exon.new, stat.y = c("sample1", "sample2"), annotation = list(p17))


## Themes
p
p + theme_alignment() 
p + theme_clear() 
p + theme_null()


# theme alignment
library(GenomicRanges) 
set.seed(1) 
N <- 100 
gr <- GRanges(seqnames = sample(c("chr1", "chr2", "chr3"), 
                                size = N, replace = TRUE), 
              IRanges(start = sample(1:300, size = N, replace = TRUE),
              width = sample(70:75, size = N,replace = TRUE)), 
              strand = sample(c("+", "-"), size = N, replace = TRUE), 
              value = rnorm(N, 10, 3), 
              score = rnorm(N, 100, 30), 
              sample = sample(c("Normal", "Tumor"),
                              size = N, replace = TRUE),
              pair = sample(letters, size = N, replace = TRUE))

seqlengths(gr) <- c(400, 1000, 500)
autoplot(gr)
autoplot(gr) + theme_tracks_sunset()

setwd("~/Dropbox/Consulting/Rick/recruitment_plot")
df <- read.table("gr_frame.txt", sep="\t", header=T)
gr <- makeGRangesFromDataFrame(df, keep.extra.columns=T)
autoplot(gr)

p <- ggbio() + circle(gr, geom = "ideo", fill = "gray70")

df2 <- read.table("gr_frame2.txt", sep="\t", header=T)
gr2 <- makeGRangesFromDataFrame(df2, keep.extra.columns=T)
seqlengths(gr2) <- unique(df2$contig_l)

p <- ggbio() + circle(gr2, geom = "ideo", fill = "gray70")
p

setwd("/Users/nielsh/Dropbox/manuscripts/SOFA/r_analysis/")
ecoli_df <- read.table("NC_000913.3_ggbio_gr.txt", sep="\t", header=T, quote = "")
ecoli_gr <- makeGRangesFromDataFrame(ecoli_df, keep.extra.columns=T)
seqlengths(ecoli_gr) <- 4640542
p <- ggbio() + circle(ecoli_gr, geom = "rect", color = "steelblue") +circle(ecoli_gr, geom = "ideo", fill = "gray70") + circle(ecoli_gr, geom = "scale", size = 2)
quartz()
p

