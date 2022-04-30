###################################################################################################

##This is an example capsule for CoBE(https://www.pmcobe.ca) executables.

##The capsule is made using software, source code and data coming from public sources and other authors.

##Please find the original links for code/data in capsule description section.

###################################################################################################
args = commandArgs(trailingOnly=TRUE)
## ----setup, echo=FALSE, results="hide"----------------------------------------
knitr::opts_chunk$set(tidy = FALSE,
                      cache = FALSE,
                      dev = "png",
                      message = FALSE, error = FALSE, warning = TRUE)

## ----quickStart, eval=FALSE---------------------------------------------------
#  dds <- DESeqDataSetFromMatrix(countData = cts,
#                                colData = coldata,
#                                design= ~ batch + condition)
#  dds <- DESeq(dds)
#  resultsNames(dds) # lists the coefficients
#  res <- results(dds, name="condition_trt_vs_untrt")
#  # or to shrink log fold changes association with condition:
#  res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")

## ----txiSetup-----------------------------------------------------------------
library("tximport")
library("readr")
library("tximportData")
#dir <- system.file("extdata", package="tximportData")
#samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
#samples$condition <- factor(rep(c("A","B"),each=3))
#rownames(samples) <- samples$run
#samples[,c("pop","center","run","condition")]

## ----loadData--------------------------------------------------------------
cts <- as.matrix(read.table(args[1]))
coldata <- read.table(args[2])

## ----reorderPasila------------------------------------------------------------
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

## ----matrixInput--------------------------------------------------------------
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
## ----prefilter----------------------------------------------------------------
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## ----factorlvl----------------------------------------------------------------
dds$condition <- factor(dds$condition, levels = c("untreated","treated"))

## ----deseq--------------------------------------------------------------------
dds <- DESeq(dds)
res <- results(dds)

## ----MAidentify, eval=FALSE---------------------------------------------------
#  idx <- identify(res$baseMean, res$log2FoldChange)
#  rownames(res)[idx]

## ----warning=FALSE------------------------------------------------------------
resultsNames(dds)
# because we are interested in treated vs untreated, we set 'coef=2'
resNorm <- lfcShrink(dds, coef=2, type="normal")
log_fold <- data.frame(lfc=resNorm@listData$log2FoldChange,row.names = resNorm@rownames)
p_value <- data.frame(p_value=resNorm@listData$pvalue,adj=resNorm@listData$padj,row.names=resNorm@rownames)
write.table(log_fold,"/results/log_fold.txt")
write.table(p_value,"/results/p_values.txt")

## ----lfcThresh----------------------------------------------------------------
par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-2.5,2.5)
resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
pdf("/results/log_folds.pdf")
plotMA(resGA, ylim=ylim,main="LFC <= 0.5 or LFC >= 0.5"); drawLines()
plotMA(resLA, ylim=ylim,main="-0.5 <= LFC <= 0.5"); drawLines()
plotMA(resG, ylim=ylim,main="LFC > 0.5"); drawLines()
plotMA(resL, ylim=ylim,main="LFC < -0.5"); drawLines()
dev.off()


