OG_minusST_high <- read.csv("~/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/4_Diff_expression/Diff_Expression_OG_3/OG_minusST_high.csv")
OG_minusST_high[,3:ncol(OG_minusST_high)]<-round(OG_minusST_high[,-c(1,2)],digits = 3)
write.csv(OG_minusST_high[,c("Id","Description","baseMean","log2FoldChange.minusSTvsAll","padj.minusSTvsAll")],"~/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/4_Diff_expression/Diff_Expression_OG_3/OG_minusST_high_latex.csv",row.names = F)

OG_minusST_down <- read.csv("~/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/4_Diff_expression/Diff_Expression_OG_3/OG_minusST_down.csv")
OG_minusST_down[,3:ncol(OG_minusST_down)]<-round(OG_minusST_down[,-c(1,2)],digits = 3)
write.csv(OG_minusST_down[,c("Id","Description","baseMean","log2FoldChange.minusSTvsAll","padj.minusSTvsAll")],"~/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/4_Diff_expression/Diff_Expression_OG_3/OG_minusST_down_latex.csv",row.names = F)

 singletons<-c(107 ,121, 38 ,61 ,21 ,21 ,42, 50 ,78 ,64, 55, 17, 33, 27, 33, 133, 29, 35)

main_sing<-c( 281 ,79 ,134)
 
 
setwd("/media/chrats/Elements/CHEESE_foodtranscriptomics/DKINKJ/milk_exp/Differential_expression/4_Diff_expresseion/Diff_Expression_ST_3/tables")
ST_func_len <- read.csv("/media/chrats/Elements/CHEESE_foodtranscriptomics/DKINKJ/milk_exp/Differential_expression/4_Diff_expresseion/Diff_Expression_ST_3/ST_func_len.csv")
flist <- list.files(pattern = "All.down")

for(j in 1:length(flist)){
  temp <- read.delim(flist[j])
  ST_func_len
  temp1<-merge(temp,ST_func_len,by.x="Id", by.y="X")
  temp2<-temp1[,c("Id","Description","baseMean","log2FoldChange","padj")]
  temp2$padj<-round(temp2$padj, digits = 3)
  temp3<- temp2[order(temp2$log2FoldChange),]
  write.csv(temp3,paste("Annotated",flist[j],sep = ""))
}

flist <- list.files(pattern = "All.up")

for(j in 1:length(flist)){
  temp <- read.delim(flist[j])
  ST_func_len
  temp1<-merge(temp,ST_func_len,by.x="Id", by.y="X")
  temp2<-temp1[,c("Id","Description","baseMean","log2FoldChange","padj")]
  temp2$padj<-round(temp2$padj, digits = 3)
  temp3<- temp2[order(temp2$log2FoldChange,decreasing = T),]
  write.csv(temp3,paste("Annotated",flist[j],sep = ""))
}

|