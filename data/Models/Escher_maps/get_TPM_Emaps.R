#############################
#metatranscriptomics
#####################
setwd("/home/chrats/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/3_count_tables")
htseq_count_minusCHCC3053 <- read.csv("~/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/3_count_tables/htseq_count_minusCHCC3053.csv")
OG_func_len <- read.csv("~/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/3_count_tables/OG_func_len.csv")
ST_func_len <- read.csv("~/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/3_count_tables/ST_func_len.csv")
OG_geneid_clean <- read.csv("~/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/3_count_tables/OG_geneid_clean.csv", header=FALSE)

library(Biostrings)
genes<-readDNAStringSet("~/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/3_count_tables/DB/CDS.fa")


names_genes<-unlist(lapply(strsplit(names(genes),"ID=|;"), function(x) x[2]))
names(genes)<-names_genes


genomes<-htseq_count_minusCHCC3053[,1]
htseq_count_minusCHCC3053<-apply(as.matrix(htseq_count_minusCHCC3053[,-1]), 2, as.numeric)
rownames(htseq_count_minusCHCC3053)<-genomes


names_genes_ord<-names_genes[match(rownames(htseq_count_minusCHCC3053),names_genes)]
w_genes_ord<-genes@ranges@width[match(rownames(htseq_count_minusCHCC3053),names_genes)]
# metaT.fullgenome_counts_norm<-sweep(metaT.fullgenome_counts, 2, colSums(metaT.fullgenome_counts), FUN="/")
# metaT.fullgenome_counts_norm<-metaT.fullgenome_counts_norm*100

#metaT.fullgenome_counts colSums(metaT.fullgenome_counts)

#colnames(metaT.fullgenome_counts_norm)<-gsub("HTTFFBGXC_||_19s00.*","",colnames(metaT.fullgenome_counts_norm))
colnames(htseq_count_minusCHCC3053)<-gsub("sample_||_mappedReads","",colnames(htseq_count_minusCHCC3053))

htseq_count_minusCHCC3053<-htseq_count_minusCHCC3053[,sort(colnames(htseq_count_minusCHCC3053))]


# Create overview file with TPM and values from the various comparisons
# Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
# Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
# Divide the RPK values by the “per million” scaling factor. This gives you TPM.
htseq_count_scalingFactor <- apply(htseq_count_minusCHCC3053, 2, function(x)sum(x)/1000000)
htseq_count_minusCHCC3053_TPM<-htseq_count_minusCHCC3053
for (i in 1:ncol(htseq_count_minusCHCC3053)){
  htseq_count_minusCHCC3053_TPM[,i]<-htseq_count_minusCHCC3053_TPM[,i]/w_genes_ord/htseq_count_scalingFactor[i]
}

#########PCA all




lactis_st<-c("CHCC5614","CHCC10675","CHCC6086")
out<-"/home/chrats/Desktop/Git/cheese-ft/Data/Models/Escher_maps/lactis/DataEscher2/"
for (i in lactis_st){
  temp<-htseq_count_minusCHCC3053_TPM[grep(i,rownames(htseq_count_minusCHCC3053_TPM)),]
  temp_all<-rowMeans(temp[,1:3])
  
  temp_all_og<-lapply(names(temp_all), function(x) unlist(apply(OG_geneid_clean, 1, function(y) if (x%in%gsub(" ","",y)==T){y[1]})))
  temp_all_og[unlist(lapply(temp_all_og , is.null))] <- NA
  
  names(temp_all)<-unlist(temp_all_og)
  data_all<-data.frame(OGs=names(temp_all),TPM=temp_all)
  write.csv(data_all, paste(out,i,"_all_TMP.csv",sep = ""),row.names = F,quote = F)

  temp_minus_cremoris<-rowMeans(temp[,7:9])
  names(temp_minus_cremoris)<-unlist(temp_all_og)
  data_minus_cremoris<-data.frame(OGs=names(temp_minus_cremoris),TPM=temp_minus_cremoris)
  write.csv(data_minus_cremoris, paste(out,i,"_minus_cremoris_TMP.csv",sep = ""),row.names = F,quote = F)
  
  temp_minus_ST<-rowMeans(temp[,4:6])
  names(temp_minus_ST)<-unlist(temp_all_og)
  data_minus_ST<-data.frame(OGs=names(temp_minus_ST),TPM=temp_minus_ST)
  write.csv(data_minus_ST, paste(out,i,"_minus_ST_TMP.csv",sep = ""),row.names = F,quote = F)
  
}
test2<-temp_all_og
test2[unlist(lapply(test2 , is.null))] <- NA
unlist(apply(OG_geneid_clean, 1, function(y) if (test%in%y==T){y[1]}))

length(table(unlist(apply(OG_geneid_clean,1, function(x) grep("CHCC10675",x,value = T)))))
length(table(unlist(apply(OG_geneid_clean,1, function(x) grep("CHCC5614",x,value = T)))))
