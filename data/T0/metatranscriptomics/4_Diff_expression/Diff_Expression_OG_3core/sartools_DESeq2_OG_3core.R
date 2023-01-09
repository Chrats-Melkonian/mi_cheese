
#install.packages("devtools")
#library(devtools)
#install_github("PF2-pasteur-fr/SARTools", build_opts="--no-resave-data")

#SARTools DESeq2 pipeline



library(devtools)
library(DESeq2)
library(SARTools)

################################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Hugo Varet
### April 20th, 2015
### designed to be executed with SARTools 1.1.0
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
rm(list=ls())                                        # remove all the objects from the R session

workDir <- "/Users/dkinkj/OneDrive - Chr Hansen/InKj/Analysis/metaT/Diff_Expression_OG_3core"      # working directory for the R session

projectName <- "FoodTrans_metaT_milk_OG_3core"                         # name of the project
author <- "InKj"                                # author of the statistical analysis/report

targetFile <- "target.txt"                           # path to the design/target file
rawDir <- "OG_count_samples_core"                                      # path to the directory containing raw counts files
featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                      "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual")

varInt <- "condition"                                    # factor of interest
condRef <- "All"                                      # reference biological condition
batch <- NULL                                        # blocking factor: NULL (default) or "batch" for example

fitType <- "parametric"                              # mean-variance relationship: "parametric" (default) or "local"
cooksCutoff <- TRUE                                # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- 0.01                                        # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors

#colors <- c("cornflowerblue","blue", "lightcoral","red3", "palegreen3", "forestgreen", "orangered", "orangered4")              # vector of colors of each biological condition on the plots

library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(10, "Paired")) (6)

colors <- c("#8B535D", "#BE939B", "#456C6F", "#85B0B3", "#B08B53", "#DAC19B", "#59596E", "#9898AC", "#245776", "#5BA2CD")

################################################################################
###                             running script                               ###
################################################################################
setwd(workDir)


# checking parameters
checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                       rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                       condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                       independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                       typeTrans=typeTrans,locfunc=locfunc,colors=colors)

# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDir, skip='1')

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# analysis with DESeq2
out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)

# PCA + clustering
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)

# summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
                                          independentFiltering=independentFiltering, 
                                          cooksCutoff=cooksCutoff, alpha=alpha)

# save image of the R session
save.image(file=paste0(projectName, ".RData"))

# generating HTML report
writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                   targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                   condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
                   independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
                   typeTrans=typeTrans, locfunc=locfunc, colors=colors)







################################################
#           Additional commands !!             #
################################################

#get normalized and transformed count matrix
counts.trans.repl <- assay(varianceStabilizingTransformation(out.DESeq2$dds))

#remove rows with an SD of 0
counts.trans.noSD <- as.matrix(counts.trans.repl[!(apply(counts.trans.repl,1,sd)==0),])

#merge samples by columns, (every third column is averaged)
counts.trans.merged <- sapply(split.default(as.data.frame(counts.trans.noSD), ((1:18)-1)%/%3 + 1), rowMeans)



#####################################
##Heatmap (Hierarchical clustering)##
#####################################
#BiocManager::install("bioDist")
library(bioDist)
library(RColorBrewer)
#install.packages("gplots")
library(gplots)
library(viridis)

#colors
hmcol <- colorRampPalette(brewer.pal(11, "PRGn"))(256) ## try 
hmcol <- colorRampPalette(viridis(11))(256)
hmcol <- colorRampPalette(c("#000099", "#0066CC", "white", "#CC66FF", "#990066"))(n = 299) # cyan-like -> magenta-like
brewer.pal(11, "RdGy")
hmcol <- colorRampPalette(c("turquoise4", "turquoise", "white", "orange", "darkorange3"))(n = 100) # cyan-like -> magenta-like


# BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral 
mydist <- function(x)cor.dist(x,abs=FALSE) #pearson
mydist2 <- function(x)dist(x,method = "euclidean")
mydist3 <- function(x)spearman.dist(x)
mydist4 <- function(x)cor.dist(x,abs=FALSE)^2 #pearson
mydist5 <- function(x)(dist(x,method = "euclidean")^2)
myhclust <- function(x){hclust(x,method='average')}
myhclust2 <- function(x){hclust(x,method='ward')} 
myhclust3 <- function(x){hclust(x,method='ward.D')} #for euclidian squared
myhclust4 <- function(x){hclust(x,method='ward.D2')} 

pdf(file="./heatmap2.pdf")
heatmap.2(counts.trans.noSD, Colv=TRUE,
          dendrogram="both", distfun=mydist,
          hclustfun=myhclust4,
          margins=c(5,1),
          sepwidth=c(0.1,0.1),colsep=c(),
          key = TRUE,
          keysize=1.5, density.info=c("none"),
          scale="row", trace="none", 
          cexCol=1, col=hmcol,cexRow = 0.1,
          lmat=rbind(c(0,3,0),c(2,1,0),c(4,0,0)),lhei=c(0.6,4,1),
          lwid=c(1.5,5,0.5),na.rm=F)
dev.off()



#PAM clustering

require(cluster)

dist.cor<-as.dist((cor(t(counts.trans.noSD))*-0.5)+0.5)
k<-20  # change this number as you wish ...
pam.cl<-pam(dist.cor, k)
par(mfrow=n2mfrow(k))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

pdf(file="pam.cluster.pdf")
for(knr in 1:k){
  given.cl<-apply(counts.trans.noSD[pam.cl$clustering==knr,],1,function(v){x<-v-mean(v); return(x/sd(x))})
  cl.size<-length(given.cl[1,])
  plot(given.cl[,1], type="l", ylim=range(given.cl),xaxt="n", xlab="", main=knr, ylab="Gene expression intensity (Z)",col=rainbow(cl.size)[1])
  axis(1, c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),c('1A', '1B', '1C', '2A', '2B', '2C', '3A', '3B', '3C', '4A', '4B', '4C', '5A', '5B', '5C', '6A', '6B', '6C'), cex.axis=0.7,las=2) # you may need to change numbers...
  for(i in 2:cl.size){
    points(given.cl[,i],type="l", col=gg_color_hue(cl.size)[i])
  }}
#dev.print(file="k8.ps")
dev.off()

write.table(pam.cl$cluster, file="pam.cluster.txt")

clusplot(x=dist.cor, k=10, diss=TRUE)
plot(pam.cl) 



#####################################
## Make overview file              ##
#####################################


minus_5614vsAll <- read.csv(file = 'tables/minus 5614vsAll.complete.txt', sep="\t")
minus_STvsAll <- read.csv(file = 'tables/minus STvsAll.complete.txt', sep="\t")

minus_6086vsAll <- read.csv(file = 'tables/minus 6086vsAll.complete.txt', sep="\t")
minus_10675vsAll <- read.csv(file = 'tables/minus 10675vsAll.complete.txt', sep="\t")
minus_SICOvsAll <- read.csv(file = 'tables/minus SICOvsAll.complete.txt', sep="\t")

OG_func_len <- read.csv(file="OG_func_len.csv")

median_geneL <- median(OG_func_len[,'Length'])


# Create overview file with norm value based on gene length and values from the various comparisons
minus_5614vsAll <-merge(minus_5614vsAll, OG_func_len, by.x='Id', by.y="X")

Total_overview <- minus_5614vsAll[,c('Id', 'Description', 'baseMean')]
Total_overview[,c('All_norm', 'minus.ST_norm', 'minus.5614_norm', 'minus.6086_norm', 'minus.10675_norm', 'minus.SICO_norm')] <- minus_5614vsAll[,c('All', 'minus.ST', 'minus.5614', 'minus.6086', 'minus.10675', 'minus.SICO')]/minus_5614vsAll[,c('Length')]*median_geneL

Total_overview[,c('FoldChange.minusSTvsAll', 'log2FoldChange.minusSTvsAll', 'stat.minusSTvsAll', 'pvalue.minusSTvsAll', 'padj.minusSTvsAll')] <- minus_STvsAll[,c('FoldChange', 'log2FoldChange', 'stat', 'pvalue', 'padj')]
Total_overview[,c('FoldChange.minus5614vsAll', 'log2FoldChange.minus5614vsAll', 'stat.minus5614vsAll', 'pvalue.minus5614vsAll', 'padj.minus5614vsAll')] <- minus_5614vsAll[,c('FoldChange', 'log2FoldChange', 'stat', 'pvalue', 'padj')]
Total_overview[,c('FoldChange.minus6086vsAll', 'log2FoldChange.minus6086vsAll', 'stat.minus6086vsAll', 'pvalue.minus6086vsAll', 'padj.minus6086vsAll')] <- minus_6086vsAll[,c('FoldChange', 'log2FoldChange', 'stat', 'pvalue', 'padj')]
Total_overview[,c('FoldChange.minus10675vsAll', 'log2FoldChange.minus10675vsAll', 'stat.minus10675vsAll', 'pvalue.minus10675vsAll', 'padj.minus10675vsAll')] <- minus_10675vsAll[,c('FoldChange', 'log2FoldChange', 'stat', 'pvalue', 'padj')]
Total_overview[,c('FoldChange.minusSICOvsAll', 'log2FoldChange.minusSICOvsAll', 'stat.minusSICOvsAll', 'pvalue.minusSICOvsAll', 'padj.minusSICOvsAll')] <- minus_SICOvsAll[,c('FoldChange', 'log2FoldChange', 'stat', 'pvalue', 'padj')]

write.csv(Total_overview,'OGcore_Total_overview_Lgene_norm.csv', row.names = FALSE)


# Create overview file with TPM and values from the various comparisons
# Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
# Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
# Divide the RPK values by the “per million” scaling factor. This gives you TPM.

perM_scalingFactor_sAll <- sum(minus_5614vsAll[,'All'])/1000000
perM_scalingFactor_sminus.ST <- sum(minus_5614vsAll[,'minus.ST'])/1000000
perM_scalingFactor_sminus.5614 <- sum(minus_5614vsAll[,'minus.5614'])/1000000
perM_scalingFactor_sminus.6086 <- sum(minus_5614vsAll[,'minus.6086'])/1000000
perM_scalingFactor_sminus.10675 <- sum(minus_5614vsAll[,'minus.10675'])/1000000
perM_scalingFactor_sminus.SICO <- sum(minus_5614vsAll[,'minus.SICO'])/1000000

#minus_5614vsAll <-merge(minus_5614vsAll, OG_func_len, by.x='Id', by.y="X")

Total_overview_TPM <- minus_5614vsAll[,c('Id', 'Description', 'baseMean')]
Total_overview_TPM[,c('All_TMP')] <- minus_5614vsAll[,c('All')]/minus_5614vsAll[,c('Length')]/perM_scalingFactor_sAll
Total_overview_TPM[,c('minus.ST_TMP')] <- minus_5614vsAll[,c('minus.ST')]/minus_5614vsAll[,c('Length')]/perM_scalingFactor_sminus.ST
Total_overview_TPM[,c('minus.5614_TMP')] <- minus_5614vsAll[,c('minus.5614')]/minus_5614vsAll[,c('Length')]/perM_scalingFactor_sminus.5614
Total_overview_TPM[,c('minus.6086_TMP')] <- minus_5614vsAll[,c('minus.6086')]/minus_5614vsAll[,c('Length')]/perM_scalingFactor_sminus.6086
Total_overview_TPM[,c('minus.10675_TMP')] <- minus_5614vsAll[,c('minus.10675')]/minus_5614vsAll[,c('Length')]/perM_scalingFactor_sminus.10675
Total_overview_TPM[,c('minus.SICO_TMP')] <- minus_5614vsAll[,c('minus.SICO')]/minus_5614vsAll[,c('Length')]/perM_scalingFactor_sminus.SICO



Total_overview_TPM[,c('FoldChange.minus5614vsAll', 'log2FoldChange.minus5614vsAll', 'stat.minus5614vsAll', 'pvalue.minus5614vsAll', 'padj.minus5614vsAll')] <- minus_5614vsAll[,c('FoldChange', 'log2FoldChange', 'stat', 'pvalue', 'padj')]
Total_overview_TPM[,c('FoldChange.minusSTvsAll', 'log2FoldChange.minusSTvsAll', 'stat.minusSTvsAll', 'pvalue.minusSTvsAll', 'padj.minusSTvsAll')] <- minus_STvsAll[,c('FoldChange', 'log2FoldChange', 'stat', 'pvalue', 'padj')]
Total_overview_TPM[,c('FoldChange.minus6086vsAll', 'log2FoldChange.minus6086vsAll', 'stat.minus6086vsAll', 'pvalue.minus6086vsAll', 'padj.minus6086vsAll')] <- minus_6086vsAll[,c('FoldChange', 'log2FoldChange', 'stat', 'pvalue', 'padj')]
Total_overview_TPM[,c('FoldChange.minus10675vsAll', 'log2FoldChange.minus10675vsAll', 'stat.minus10675vsAll', 'pvalue.minus10675vsAll', 'padj.minus10675vsAll')] <- minus_10675vsAll[,c('FoldChange', 'log2FoldChange', 'stat', 'pvalue', 'padj')]
Total_overview_TPM[,c('FoldChange.minusSICOvsAll', 'log2FoldChange.minusSICOvsAll', 'stat.minusSICOvsAll', 'pvalue.minusSICOvsAll', 'padj.minusSICOvsAll')] <- minus_SICOvsAll[,c('FoldChange', 'log2FoldChange', 'stat', 'pvalue', 'padj')]

write.csv(Total_overview_TPM,'OGcore_Total_overview_TPM_norm.csv', row.names = FALSE)



##################

Total_overview[]

Total_overview[Total_overview$padj.minus5614vsAll < 0.01,]







