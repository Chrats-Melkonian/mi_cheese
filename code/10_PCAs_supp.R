setwd("/home/chrats/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/4_Diff_expression/Diff_Expression_OG_3")
minus_5614vsAll <- read.csv(file = 'tables/minus 5614vsAll.complete.txt', sep="\t")

OG_func_len <- read.csv(file="OG_func_len.csv")

median_geneL <- median(OG_func_len[,'Length'])


# Create overview file with norm value based on gene length and values from the various comparisons
minus_5614vsAll <-merge(minus_5614vsAll, OG_func_len, by.x='Id', by.y="X")
Total<-minus_5614vsAll[,grep(pattern = "norm.",colnames(minus_5614vsAll))]


Total_TPM<-Total/minus_5614vsAll$Length*median_geneL

colnames(Total_TPM)<-gsub("norm.sample_","",colnames(Total_TPM))
df_f <- Total_TPM[apply(na.omit(Total_TPM), 1, var, na.rm=TRUE) != 0,]
df_f<- na.omit(df_f)
pca = prcomp(t(df_f), center = T, scale. = T)
pca
Conditions<-mapvalues(substr(colnames(df_f),1,1), 
                  from=c("1","2","3","4","5","6"), 
                  to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
colors_t0<-rev(c('#a50026','#1b7837','#abd9e9','#4575b4','#313695','#fdae61'))
x1<- -3.2
x2<-3.2
LL_OG_PCA<-ggbiplot(pca,var.axes = F,groups=factor(Conditions),ellipse=T)+
  scale_color_manual(name="Conditions",values=colors_t0)+
  geom_point(aes(colour=factor(Conditions)), size = 3) +
  theme(legend.direction ="horizontal", 
        legend.position = "bottom")+ ggtitle(expression(paste(italic("Lactis"), " OGs")))+ylim(x1,x2)+xlim(x1,x2)
LL_OG_PCA
#########################

################# S. thermophilus 

####################
setwd("/home/chrats/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/4_Diff_expression/Diff_Expression_ST_3")

minus_5614vsAll <- read.csv(file = 'tables/minus 5614vsAll.complete.txt', sep="\t")

ST_func_len <- read.csv(file="ST_func_len.csv")
median_geneL <- median(ST_func_len[,'Length'])


# Create overview file with norm value based on gene length and values from the various comparisons
minus_5614vsAll <-merge(minus_5614vsAll, ST_func_len, by.x='Id', by.y="X")

Total_ST<-minus_5614vsAll[,grep(pattern = "norm.",colnames(minus_5614vsAll))]


Total_STn<-Total_ST/minus_5614vsAll$Length*median_geneL

colnames(Total_STn)<-gsub("norm.sample_","",colnames(Total_STn))
df_f <- Total_STn[apply(na.omit(Total_STn), 1, var, na.rm=TRUE) != 0,]
df_f<- na.omit(df_f)
pca = prcomp(t(df_f), center = T, scale. = T)
pca
Conditions1<-mapvalues(substr(colnames(df_f),1,1), 
                  from=c("1","2","3","4","5","6"), 
                  to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
colors_t1<-rev(c('#a50026','#abd9e9','#4575b4','#313695','#fdae61'))
ST_OG_PCA<-ggbiplot(pca,var.axes = F,groups=factor(Conditions1),ellipse=T)+
  scale_color_manual(name="Conditions", values=colors_t1)+
  geom_point(aes(colour=factor(Conditions1)), size = 3) +  
  theme(legend.direction ="horizontal", legend.position = "bottom")+ 
  ggtitle(expression(paste(italic("S. thermophilus"))))+ylim(x1,x2)+xlim(x1,x2)
ST_OG_PCA
#########################

##########################33 all together 

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
names(w_genes_ord)<-names_genes_ord
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
htseq_count_minusCHCC3053<-htseq_count_minusCHCC3053[-c((nrow(htseq_count_minusCHCC3053)-4):nrow(htseq_count_minusCHCC3053)),]
w_genes_ord<-w_genes_ord[-c((length(w_genes_ord)-4):length(w_genes_ord))]
htseq_count_minusCHCC3053_TPM<-htseq_count_minusCHCC3053
htseq_count_minusCHCC3053_TPM<-apply(htseq_count_minusCHCC3053_TPM, 2, function(x) x/w_genes_ord)
htseq_count_scalingFactor <- apply(htseq_count_minusCHCC3053_TPM, 2, function(x)sum(x)/1000000)
for (i in 1:ncol(htseq_count_minusCHCC3053)){
  htseq_count_minusCHCC3053_TPM[,i]<-htseq_count_minusCHCC3053_TPM[,i]/htseq_count_scalingFactor[i]
}



#########PCA all
df_f <- htseq_count_minusCHCC3053_TPM[apply(na.omit(htseq_count_minusCHCC3053_TPM), 1, var, na.rm=TRUE) != 0,]
df_f<- na.omit(df_f)
pca = prcomp(t(df_f), center = T, scale. = T)
pca
groups<-mapvalues(substr(colnames(df_f),1,1), 
                  from=c("1","2","3","4","5","6"), 
                  to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
colors_t0<-rev(c('#a50026','#1b7837','#abd9e9','#4575b4','#313695','#fdae61'))
all_PCA<-ggbiplot(pca,var.axes = F,groups=factor(groups),ellipse=T)+
  scale_color_manual(values=colors_t0)+
  geom_point(aes(colour=factor(groups)), size = 3) + ggtitle("The whole SLAB culture")+ylim(x1,x2)+xlim(x1,x2)
all_PCA


#verify ST pattern

clean_ids<-gsub("-.*|_.*","", rownames(htseq_count_minusCHCC3053))
ST_ids<-which(clean_ids%in%"NLANBGKA"==T)

w_genes_ST<-w_genes_ord[ST_ids]
htseq_count_ST_TPM<-htseq_count_minusCHCC3053[ST_ids,]

htseq_count_ST_TPM<-apply(htseq_count_ST_TPM, 2, function(x) x/w_genes_ST)
htseq_count_scalingFactor <- apply(htseq_count_ST_TPM, 2, function(x)sum(x)/1000000)
for (i in 1:ncol(htseq_count_ST_TPM)){
  htseq_count_ST_TPM[,i]<-htseq_count_ST_TPM[,i]/htseq_count_scalingFactor[i]
}

df_f_ST <- htseq_count_ST_TPM[apply(na.omit(htseq_count_ST_TPM), 1, var, na.rm=TRUE) != 0,]
df_f_ST<- na.omit(df_f_ST)[,-c(4:6)]

pca_ST = prcomp(t(df_f_ST), center = T, scale. = T)
pca_ST
groups<-mapvalues(substr(colnames(df_f_ST),1,1), 
                  from=c("1","3","4","5","6"), 
                  to=c("All","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
colors_t_st<-rev(c('#a50026','#abd9e9','#4575b4','#313695','#fdae61'))
st_test<-ggbiplot(pca_ST,var.axes = F,groups=factor(groups),ellipse=T)+
  scale_color_manual(values=colors_t_st)+
  geom_point(aes(colour=factor(groups)), size = 3) +ylim(x1,x2)+xlim(x1,x2)
st_test


########################3 CHCC10675 - CHCC6086 CHCC5614

clean_ids<-gsub("-.*|_.*","", rownames(htseq_count_minusCHCC3053))
id<-which(clean_ids%in%"CHCC5614"==T)

w_genes_l<-w_genes_ord[id]
htseq_count_l_TPM<-htseq_count_minusCHCC3053[id,]

htseq_count_l_TPM<-apply(htseq_count_l_TPM, 2, function(x) x/w_genes_l)
htseq_count_scalingFactor <- apply(htseq_count_l_TPM, 2, function(x)sum(x)/1000000)
for (i in 1:ncol(htseq_count_l_TPM)){
  htseq_count_l_TPM[,i]<-htseq_count_l_TPM[,i]/htseq_count_scalingFactor[i]
}

df_f_l <- htseq_count_l_TPM[apply(na.omit(htseq_count_l_TPM), 1, var, na.rm=TRUE) != 0,]
df_f_l<- na.omit(df_f_l)[,-c(7:9)]

pca = prcomp(t(df_f_l), center = T, scale. = T)
pca
groups<-mapvalues(substr(colnames(df_f_l),1,1), 
                  from=c("1","2","3","4","5","6"), 
                  to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
colors_t0<-rev(c('#a50026','#1b7837','#abd9e9','#4575b4','#313695'))

cremoris_PCA<-ggbiplot(pca,var.axes = F,groups=factor(groups),ellipse=T)+
  scale_color_manual(values=colors_t0)+ylim(x1,x2)+xlim(x1,x2)+
  geom_point(aes(colour=factor(groups)), size = 3) + ggtitle(expression(paste(italic("L. lactis subsp. cremoris"), " 5614")))
cremoris_PCA

########################3 CHCC10675 - CHCC6086 CHCC5614

clean_ids<-gsub("-.*|_.*","", rownames(htseq_count_minusCHCC3053))
id<-which(clean_ids%in%"CHCC6086"==T)

w_genes_l<-w_genes_ord[id]
htseq_count_l_TPM<-htseq_count_minusCHCC3053[id,]

htseq_count_l_TPM<-apply(htseq_count_l_TPM, 2, function(x) x/w_genes_l)
htseq_count_scalingFactor <- apply(htseq_count_l_TPM, 2, function(x)sum(x)/1000000)
for (i in 1:ncol(htseq_count_l_TPM)){
  htseq_count_l_TPM[,i]<-htseq_count_l_TPM[,i]/htseq_count_scalingFactor[i]
}


df_f_l <- htseq_count_l_TPM[apply(na.omit(htseq_count_l_TPM), 1, var, na.rm=TRUE) != 0,]
df_f_l<- na.omit(df_f_l)[,-c(10:12)]

df_f_l<-df_f_l[!apply(df_f_l, 1, function(x) all(x==0)),]

pca = prcomp(t(df_f_l), center = T, scale. = T)
pca
groups<-mapvalues(substr(colnames(df_f_l),1,1), 
                  from=c("1","2","3","4","5","6"), 
                  to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
colors_t0<-rev(c('#a50026','#1b7837','#abd9e9','#313695','#fdae61'))

lactis6_PCA<-ggbiplot(pca,var.axes = F,groups=factor(groups),ellipse=T)+
  scale_color_manual(values=colors_t0)+ylim(x1,x2)+xlim(x1,x2)+
  geom_point(aes(colour=factor(groups)), size = 3)  + ggtitle(expression(paste(italic("L. lactis subsp. lactis"), " 6086")))
lactis6_PCA

########################3 CHCC10675 - CHCC6086 CHCC5614

clean_ids<-gsub("-.*|_.*","", rownames(htseq_count_minusCHCC3053))
id<-which(clean_ids%in%"CHCC10675"==T)

w_genes_l<-w_genes_ord[id]
htseq_count_l_TPM<-htseq_count_minusCHCC3053[id,]

htseq_count_l_TPM<-apply(htseq_count_l_TPM, 2, function(x) x/w_genes_l)
htseq_count_scalingFactor <- apply(htseq_count_l_TPM, 2, function(x)sum(x)/1000000)
for (i in 1:ncol(htseq_count_l_TPM)){
  htseq_count_l_TPM[,i]<-htseq_count_l_TPM[,i]/htseq_count_scalingFactor[i]
}


df_f_l <- htseq_count_l_TPM[apply(na.omit(htseq_count_l_TPM), 1, var, na.rm=TRUE) != 0,]
df_f_l<- na.omit(df_f_l)[,-c(13:15)]

df_f_l<-df_f_l[!apply(df_f_l, 1, function(x) all(x==0)),]

pca = prcomp(t(df_f_l), center = T, scale. = T)
pca
groups<-mapvalues(substr(colnames(df_f_l),1,1), 
                  from=c("1","2","3","4","5","6"), 
                  to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
colors_t0<-rev(c('#a50026','#1b7837','#abd9e9','#4575b4','#fdae61'))

lactis10_PCA<-ggbiplot(pca,var.axes = F,groups=factor(groups),ellipse=T)+
  scale_color_manual(values=colors_t0)+ylim(x1,x2)+xlim(x1,x2)+
  geom_point(aes(colour=factor(groups)), size = 3) + ggtitle(expression(paste(italic("L. lactis subsp. lactis"), " 10675")))

lactis10_PCA


########################3 CHCC10675 - CHCC6086 CHCC5614

clean_ids<-gsub("-.*|_.*","", rownames(htseq_count_minusCHCC3053))
id<-which(clean_ids%in%c("CHCC10675","NLANBGKA","CHCC6086","CHCC5614")==F)

w_genes_l<-w_genes_ord[id]
htseq_count_l_TPM<-htseq_count_minusCHCC3053[id,]

htseq_count_l_TPM<-apply(htseq_count_l_TPM, 2, function(x) x/w_genes_l)
htseq_count_scalingFactor <- apply(htseq_count_l_TPM, 2, function(x)sum(x)/1000000)
for (i in 1:ncol(htseq_count_l_TPM)){
  htseq_count_l_TPM[,i]<-htseq_count_l_TPM[,i]/htseq_count_scalingFactor[i]
}


df_f_l <- htseq_count_l_TPM[apply(na.omit(htseq_count_l_TPM), 1, var, na.rm=TRUE) != 0,]
df_f_l<- na.omit(df_f_l)[,-c(16:18)]

df_f_l<-df_f_l[!apply(df_f_l, 1, function(x) all(x==0)),]

pca = prcomp(t(df_f_l), center = T, scale. = T)
pca
groups<-mapvalues(substr(colnames(df_f_l),1,1), 
                  from=c("1","2","3","4","5","6"), 
                  to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
colors_t0<-rev(c('#a50026','#1b7837','#4575b4','#313695','#fdae61'))

lactis_sico_PCA<-ggbiplot(pca,var.axes = F,groups=factor(groups),ellipse=T)+
  scale_color_manual(values=colors_t0)+ylim(x1,x2)+xlim(x1,x2)+
  geom_point(aes(colour=factor(groups)), size = 3) + ggtitle(expression(paste(italic("L. lactis"), " mixture of strains (SICO)")))

lactis_sico_PCA


###############3
legend1 <- get_legend(all_PCA)
prow <- plot_grid(all_PCA + theme(legend.position="none"),
                  LL_OG_PCA  + theme(legend.position="none"),
                  cremoris_PCA + theme(legend.position="none"),
                  lactis6_PCA + theme(legend.position="none") ,
                  lactis10_PCA  + theme(legend.position="none"),
                  lactis_sico_PCA+ theme(legend.position="none") ,
                  ST_OG_PCA + theme(legend.position="none"),
                  legend1,
                  labels=c("a", "b","c","d","e","f","g",""), ncol = 4, nrow = 2, label_fontface = "bold",label_size=20,align = "hv")

prow
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/SUPP/PCAs.pdf",prow, units="in", width=16, height=8, dpi=600)


