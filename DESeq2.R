#Purpose of this exercise:  will learn how to analyse a count table, such as arising from a summarised RNA-Seq experiment, for differentially expressed genes.
# purpose of the experiment was to investigate the role of the estrogen receptor in parathy- roid tumors.
library( "DESeq2" )
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library( "parathyroidSE" )
library( "parathyroidSE" )
library( "parathyroidSE" )
library("parathyroidSE")
source("http://bioconductor.org/biocLite.R")
biocLite("parathyroidSE")
data("parathyroidGenesSE")
source("http://bioconductor.org/biocLite.R")
biocLite("parathyroidSE")
data("parathyroidGenesSE")
head( assay( parathyroidGenesSE ) )
library("parathyroidGenesSE")
source("http://bioconductor.org/biocLite.R")
biocLite("parathyroiGenesdSE")
source("http://bioconductor.org/biocLite.R")
biocLite("parathyroiGenesSE")
source("http://bioconductor.org/biocLite.R")
biocLite("parathyroidSE")
library("parathyroidGenesSE")
library("parathyroidSE")
data("parathyroidGenesSE")
head( assay( parathyroidGenesSE ) )
# in the table above: each row represents an Ensembl gene, each column a sequenced RNA library
nrow(parathyroidGenesSE)
rowdata(parathyroidGenesSE)
rowData(parathyroidGenesSE)
as.data.frame( colData(parathyroidGenesSE)[,c("sample","patient","treatment","time")] )
# as.data.frame forces R to show us the full list within eachcolumn
# The following: an example of a typical preparatory data manipulation task done with elementary R functions
allColSamples <- colData(parathyroidGenesSE)$sample
sp <- split( seq(along=allColSamples), allColSamples )
# use the function split to see which columns need to be collapsed
countdata <- sapply(sp, function(columns)
rowSums( assay(parathyroidGenesSE)[,columns,drop=FALSE] ) )
head(countdata)
coldata <- colData(parathyroidGenesSE)[sapply(sp, `[`, 1),]
rownames(coldata) <- coldata$sample
coldata
coldata <- coldata[ , c( "patient", "treatment", "time" ) ]
head( coldata )
rowdata <- rowData(parathyroidGenesSE)
rowdata
#coldata: only keep the columns needed for your analysis
# ∼ patient + treatment, which means that we want to test for the effect of treatment (the last factor), controlling for the effect of patient (the first factor)
ddsFull <- DESeqDataSetFromMatrix(
countData = countdata,
colData = coldata,
design = ~ patient + treatment,
rowData = rowdata)
ddsFull
library( "DESeq2" )
ddsFull <- DESeqDataSetFromMatrix(
countData = countdata,
colData = coldata,
design = ~ patient + treatment,
rowData = rowdata)
ddsFull
# Next: analyze a subset of the samples, taken after 48 hours, with either control or DPN treatment, taking into account the multifactor design
# First: subset the relevant columns from the full dataset
ddsFull <- DESeqDataSetFromMatrix(
countData = countdata,
colData = coldata,
design = ~ patient + treatment,
rowData = rowdata)
ddsFull
stopifnot(sum(assay(parathyroidGenesSE)) == sum(counts(ddsFull)))
dds <- ddsFull[ , colData(ddsFull)$treatment %in% c("Control","DPN") &
colData(ddsFull)$time == "48h" ]
dds$patient <- factor(dds$patient)
dds$treatment <- factor(dds$treatment)
dds$treatment <- relevel( dds$treatment, "Control" )
colData(dds)
# The function reveal: make sure that Control is the first level in the treatment factor, so that the log2 fold changes are calculated as treatment over control
dds$treatment <- relevel( dds$treatment, "Control" )
colData(dds)
dds <- DESeq(dds)
res <- results(dds)
mcols(res)
# column log2FoldChange is the effect size estimate --> tells us how much the gene’s expression seems to have changed due to treatment with DPN in comparison to control
all.equal(res$baseMean, rowMeans(counts(dds)))
resSig <- res[ which(res$padj < 0.1 ), ]
head( resSig[ order( resSig$log2FoldChange ), ] )
sum( res$padj < 0.1, na.rm=TRUE )
tail( resSig[ order( resSig$log2FoldChange ), ] )
table(sign(resSig$log2FoldChange))
# visualize dispersion values
plotDispEsts( dds )
plotMA(dds, ylim = c( -1.5, 1.5) )
hist( res$pvalue, breaks=100 )
# MA plot suggests that for genes with less than one or two counts per sample, averaged over all samples, there is no real inferential power
filterThreshold <- 2.0
keep <- rowMeans( counts( dds, normalized=TRUE ) ) > filterThreshold
table( keep )
min( res$padj[!keep], na.rm=TRUE )
table( p.adjust( res$pvalue, method="BH" ) < .1 )
table( p.adjust( res$pvalue[keep], method="BH" ) < .1 )
library( "org.Hs.eg.db" )
columns(org.Hs.eg.db)
convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
stopifnot( inherits( db, "AnnotationDb" ) )
ifMultiple <- match.arg( ifMultiple )
suppressWarnings( selRes <- AnnotationDbi::select(
db, keys=ids, keytype=fromKey, cols=c(fromKey,toKey) ) )
if( ifMultiple == "putNA" ) {
duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
return( selRes[ match( ids, selRes[,1] ), 2 ] )
}
res$symbol <- convertIDs( row.names(res), "ENSEMBL", "SYMBOL", org.Hs.eg.db )
res
write.csv( as.data.frame(res), file="results.csv" )
library( "reactome.db" )
source("http://bioconductor.org/biocLite.R")
biocLite("reactome.db")
