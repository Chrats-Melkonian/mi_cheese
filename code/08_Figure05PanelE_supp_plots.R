####################
#metatranscriptomics
####################

htseq_count_minusCHCC3053 <- read.csv("./mi_cheese/data/T0/metatranscriptomics/3_count_tables/htseq_count_minusCHCC3053.csv")
OG_func_len <- read.csv("./mi_cheese/data/T0/metatranscriptomics/3_count_tables/OG_func_len.csv")
ST_func_len <- read.csv("./mi_cheese/data/T0/metatranscriptomics/3_count_tables/ST_func_len.csv")
OG_geneid_clean <- read.csv("./mi_cheese/data/T0/metatranscriptomics/3_count_tables/OG_geneid_clean.csv", header=FALSE)


library(Biostrings)
genes<-readDNAStringSet("./mi_cheese/data/T0/metatranscriptomics/3_count_tables/DB/CDS.fa")


names_genes<-unlist(lapply(strsplit(names(genes),"ID=|;"), function(x) x[2]))
names(genes)<-names_genes


genomes<-htseq_count_minusCHCC3053[,1]
htseq_count_minusCHCC3053<-apply(as.matrix(htseq_count_minusCHCC3053[,-1]), 2, as.numeric)
rownames(htseq_count_minusCHCC3053)<-genomes


names_genes_ord<-names_genes[match(rownames(htseq_count_minusCHCC3053),names_genes)]
w_genes_ord<-genes@ranges@width[match(rownames(htseq_count_minusCHCC3053),names_genes)]

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
library(reshape2)
library(ggplot2)
library(plyr)




#OG1270 OG257 NLANBGKA_01275 Alpha-acetolactate decarboxylase
selection_og<-c("OG1270","OG257","OG5779")
data_ids<-unlist(lapply(selection_og, function(x) as.character(OG_geneid_clean[OG_geneid_clean$V1%in%x,])))
data_ids<-data_ids[-which(data_ids=="")]
data_ids<-gsub(" ","", data_ids)


OG_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%data_ids[-grep("^OG.*",data_ids)],]
ST_data1<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_01275",]
data_og<-rbind(OG_data,ST_data1)
rownames(data_og)[match(c("ST_data1"),rownames(data_og))]<-c("NLANBGKA_01275")
################## 
data<-melt(data_og)
data$clean_ids<-gsub("-.*|_01275","", data$Var1)

data$clean_ids[!data$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"LB"
data$condition<-mapvalues(substr(data$Var2,1,1), 
                          from=c("1","2","3","4","5","6"), 
                          to=c("All","-ST","-LC","-LLm2","-LLm1","-LB"))
data$clean_ids_names<-mapvalues(data$clean_ids, 
                                from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                to=c("LC","LLm2","LLm1","LB","ST"))
#colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')
colors_t0s<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837')
size_all<-25
data$condition<-factor(data$condition,levels = c( "All", "-LLm1","-LC","-ST","-LLm2","-LB"))
data$clean_ids_names<-factor(data$clean_ids_names,levels = c("LC","LLm1","LLm2","LB","ST"))

at1<-ggplot(data,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="Alpha-acetolactate decarboxylase (ACLDC)")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)+ theme(legend.position = "none")
at1
# ggsave(".../PLOTS/diacetyl_acetoin/acetolactate_decarboxylase.pdf",at, units="in", width=10, height=6, dpi=600)

###################################################3
#######################################
#OG935 OG1437 OG2435 OG6820 NLANBGKA_01274 NLANBGKA_00958 NLANBGKA_00959 acetolactate synthase
selection_og<-c("OG935","OG1437","OG2435","OG6820")
data_ids<-unlist(lapply(selection_og, function(x) as.character(OG_geneid_clean[OG_geneid_clean$V1%in%x,])))
data_ids<-data_ids[-which(data_ids=="")]
data_ids<-gsub(" ","", data_ids)

OG_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%data_ids[-grep("^OG.*",data_ids)],]
ST_data1<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_01274",]
ST_data2<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00958",]
ST_data3<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00959",]
data_og<-rbind(OG_data,ST_data1,ST_data2,ST_data3)
rownames(data_og)[match(c("ST_data1","ST_data2","ST_data3"),rownames(data_og))]<-c("NLANBGKA_01274","NLANBGKA_00958","NLANBGKA_00959")
################## 
data<-melt(data_og)
data$clean_ids<-gsub("-.*|_01274|_00958|_00959","", data$Var1)


data$clean_ids[!data$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"LB"
data$condition<-mapvalues(substr(data$Var2,1,1), 
                          from=c("1","2","3","4","5","6"), 
                          to=c("All","-ST","-LC","-LLm2","-LLm1","-LB"))
data$clean_ids_names<-mapvalues(data$clean_ids, 
                                from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                to=c("LC","LLm2","LLm1","LB","ST"))
#colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')
colors_t0s<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837')
size_all<-25
data$condition<-factor(data$condition,levels = c( "All", "-LLm1","-LC","-ST","-LLm2","-LB"))
data$clean_ids_names<-factor(data$clean_ids_names,levels = c("LC","LLm1","LLm2","LB","ST"))

at2<-ggplot(data,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="Acetolactate synthase (ACLS)")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)+ theme(legend.position = "none")
at2
#ggsave(".../acetolactate_synthase.pdf",at, units="in", width=10, height=6, dpi=600)



# OG2171,Diacetyl reductase [(S)-acetoin forming],253
# OG2232,Diacetyl reductase [(S)-acetoin forming],256
# OG3096 OG3097 OG4062
# 
# NLANBGKA_01265
# NLANBGKA_01266
########################### 
# Diacetyl reductase
selection_og<-c("OG2171","OG2232","OG3096","OG3097","OG4062")
#OG_geneid_clean[OG_geneid_clean$V1%in%Amo_transporter_og,]
data_ids<-unlist(lapply(selection_og, function(x) as.character(OG_geneid_clean[OG_geneid_clean$V1%in%x,])))
data_ids<-data_ids[-which(data_ids=="")]
data_ids<-gsub(" ","", data_ids)

OG_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%data_ids[-grep("^OG.*",data_ids)],]
ST_data1<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_01265",]
ST_data2<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_01266",]
data_og<-rbind(OG_data,ST_data1,ST_data2)
rownames(data_og)[match(c("ST_data1","ST_data2"),rownames(data_og))]<-c("NLANBGKA_01265","NLANBGKA_01265")
################## 
data<-melt(data_og)
data$clean_ids<-gsub("-.*|_01265|_01265","", data$Var1)

data$clean_ids[!data$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"LB"
data$condition<-mapvalues(substr(data$Var2,1,1), 
                          from=c("1","2","3","4","5","6"), 
                          to=c("All","-ST","-LC","-LLm2","-LLm1","-LB"))
data$clean_ids_names<-mapvalues(data$clean_ids, 
                                from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                to=c("LC","LLm2","LLm1","LB","ST"))
#colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')
colors_t0s<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837')
size_all<-25
data$condition<-factor(data$condition,levels = c( "All", "-LLm1","-LC","-ST","-LLm2","-LB"))
data$clean_ids_names<-factor(data$clean_ids_names,levels = c("LC","LLm1","LLm2","LB","ST"))
at3<-ggplot(data,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="Diacetyl reductase (ACTD)")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)+ theme(legend.position = "none")
at3
#ggsave(".../Diacetyl_reductase.pdf",at, units="in", width=10, height=6, dpi=600)



########################### 
# unormalized check Sorbitol dehydrogenase
Sorbitol<-"OG4434"
OG_geneid_clean[OG_geneid_clean$V1%in%Sorbitol,]
OG1286_ids<-as.character(OG_geneid_clean[OG_geneid_clean$V1%in%Sorbitol,])
OG1286_ids<-OG1286_ids[-which(OG1286_ids=="")]
OG1286_ids<-gsub(" ","", OG1286_ids)

OG1286_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%OG1286_ids[-1],]

library(reshape2)
library(ggplot2)
library(plyr)
################## 
data<-melt(OG1286_data)
data$clean_ids<-gsub("-.*","", data$Var1)

data$clean_ids[!data$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"LB"
data$condition<-mapvalues(substr(data$Var2,1,1), 
                          from=c("1","2","3","4","5","6"), 
                          to=c("All","-ST","-LC","-LLm2","-LLm1","-LB"))
data$clean_ids_names<-mapvalues(data$clean_ids, 
                                from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                to=c("LC","LLm2","LLm1","LB","ST"))
#colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')
colors_t0s<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837')
size_all<-25
data$condition<-factor(data$condition,levels = c( "All", "-LLm1","-LC","-ST","-LLm2","-LB"))
data$clean_ids_names<-factor(data$clean_ids_names,levels = c("LC","LLm1","LLm2","LB","ST"))
at4<-ggplot(data,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="2,3-butanediol dehydrogenase")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)+ theme(legend.position = "none")
at4
#ggsave(".../Sorbitol_dehydrogenase.pdf",at, units="in", width=10, height=6, dpi=600)


# NLANBGKA_00960
# Ketol-acid reductoisomerase (NADP(+))
# OG789
selection_og<-c("OG789")
data_ids<-unlist(lapply(selection_og, function(x) as.character(OG_geneid_clean[OG_geneid_clean$V1%in%x,])))
data_ids<-data_ids[-which(data_ids=="")]
data_ids<-gsub(" ","", data_ids)

OG_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%data_ids[-grep("^OG.*",data_ids)],]
ST_data1<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00960",]
data_og<-rbind(OG_data,ST_data1)
rownames(data_og)[match(c("ST_data1"),rownames(data_og))]<-c("NLANBGKA_00960")
################## 
data<-melt(data_og)
data$clean_ids<-gsub("-.*|_00960","", data$Var1)

data$clean_ids[!data$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"LB"
data$condition<-mapvalues(substr(data$Var2,1,1), 
                          from=c("1","2","3","4","5","6"), 
                          to=c("All","-ST","-LC","-LLm2","-LLm1","-LB"))
data$clean_ids_names<-mapvalues(data$clean_ids, 
                                from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                to=c("LC","LLm2","LLm1","LB","ST"))
#colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')
colors_t0s<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837')
size_all<-25
data$condition<-factor(data$condition,levels = c( "All", "-LLm1","-LC","-ST","-LLm2","-LB"))
data$clean_ids_names<-factor(data$clean_ids_names,levels = c("LC","LLm1","LLm2","LB","ST"))

at5<-ggplot(data,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="Ketol-acid reductoisomerase (KARA1)")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)+ theme(legend.position = "none")
at5
#ggsave(".../Ketol-acid_reductoisomerase.pdf",at, units="in", width=10, height=6, dpi=600)
#2-hydroxy-2-methyl-3-oxobutanoate


#Citrate synthase
# OG4497 OG4339 NLANBGKA_00090
selection_og<-c("OG4497","OG4339")
data_ids<-unlist(lapply(selection_og, function(x) as.character(OG_geneid_clean[OG_geneid_clean$V1%in%x,])))
data_ids<-data_ids[-which(data_ids=="")]
data_ids<-gsub(" ","", data_ids)

OG_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%data_ids[-1],]
ST_data1<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00090",]
data_og<-rbind(OG_data,ST_data1)
rownames(data_og)[match(c("ST_data1"),rownames(data_og))]<-c("NLANBGKA_00090")
################## 
data<-melt(data_og)
data$clean_ids<-gsub("-.*|_00090","", data$Var1)

data$clean_ids[!data$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"LB"
data$condition<-mapvalues(substr(data$Var2,1,1), 
                          from=c("1","2","3","4","5","6"), 
                          to=c("All","-ST","-LC","-LLm2","-LLm1","-LB"))
data$clean_ids_names<-mapvalues(data$clean_ids, 
                                from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                to=c("LC","LLm2","LLm1","LB","ST"))
#colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')
colors_t0s<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837')
size_all<-25
data$condition<-factor(data$condition,levels = c( "All", "-LLm1","-LC","-ST","-LLm2","-LB"))
data$clean_ids_names<-factor(data$clean_ids_names,levels = c("LC","LLm1","LLm2","LB","ST"))

at6<-ggplot(data,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="Citrate synthase (CS)")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)+ theme(legend.position = "none")
at6
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/diacetyl_acetoin/Citrate_synthase.pdf",at, units="in", width=10, height=6, dpi=600)


#,Aconitate hydratase A
# OG1 OG4350 OG6501 NLANBGKA_00091
selection_og<-c("OG4497","OG4339")
data_ids<-unlist(lapply(selection_og, function(x) as.character(OG_geneid_clean[OG_geneid_clean$V1%in%x,])))
data_ids<-data_ids[-which(data_ids=="")]
data_ids<-gsub(" ","", data_ids)

OG_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%data_ids[-grep("^OG.*",data_ids)],]
ST_data1<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00091",]
data_og<-rbind(OG_data,ST_data1)
rownames(data_og)[match(c("ST_data1"),rownames(data_og))]<-c("NLANBGKA_00091")
################## 
data<-melt(data_og)
data$clean_ids<-gsub("-.*|_00091","", data$Var1)

data$clean_ids[!data$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"LB"
data$condition<-mapvalues(substr(data$Var2,1,1), 
                          from=c("1","2","3","4","5","6"), 
                          to=c("All","-ST","-LC","-LLm2","-LLm1","-LB"))
data$clean_ids_names<-mapvalues(data$clean_ids, 
                                from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                to=c("LC","LLm2","LLm1","LB","ST"))
#colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')
colors_t0s<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837')
size_all<-25
data$condition<-factor(data$condition,levels = c( "All", "-LLm1","-LC","-ST","-LLm2","-LB"))
data$clean_ids_names<-factor(data$clean_ids_names,levels = c("LC","LLm1","LLm2","LB","ST"))

at7<-ggplot(data,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="Aconitate hydratase (ACONT)")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)+ theme(legend.position = "none")
at7
#ggsave(".../Aconitate_hydratase.pdf",at, units="in", width=10, height=6, dpi=600)



#Isocitrate dehydrogenase
# OG946 NLANBGKA_00089
selection_og<-c("OG946")
data_ids<-unlist(lapply(selection_og, function(x) as.character(OG_geneid_clean[OG_geneid_clean$V1%in%x,])))
data_ids<-data_ids[-which(data_ids=="")]
data_ids<-gsub(" ","", data_ids)

OG_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%data_ids[-grep("^OG.*",data_ids)],]
ST_data1<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00089",]
data_og<-rbind(OG_data,ST_data1)
rownames(data_og)[match(c("ST_data1"),rownames(data_og))]<-c("NLANBGKA_00089")
################## 
data<-melt(data_og)
data$clean_ids<-gsub("-.*|_00089","", data$Var1)

data$clean_ids[!data$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"LB"
data$condition<-mapvalues(substr(data$Var2,1,1), 
                          from=c("1","2","3","4","5","6"), 
                          to=c("All","-ST","-LC","-LLm2","-LLm1","-LB"))
data$clean_ids_names<-mapvalues(data$clean_ids, 
                                from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                to=c("LC","LLm2","LLm1","LB","ST"))
#colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')
colors_t0s<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837')
size_all<-25
data$condition<-factor(data$condition,levels = c( "All", "-LLm1","-LC","-ST","-LLm2","-LB"))
data$clean_ids_names<-factor(data$clean_ids_names,levels = c("LC","LLm1","LLm2","LB","ST"))
at8<-ggplot(data,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="Isocitrate dehydrogenase (ICDHyr)")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)+ theme(legend.position = "none")
at8
#ggsave(".../Isocitrate_dehydrogenase.pdf",at, units="in", width=10, height=6, dpi=600)



#Pyruvate carboxylase
# OG1378 
selection_og<-c("OG1378")
data_ids<-unlist(lapply(selection_og, function(x) as.character(OG_geneid_clean[OG_geneid_clean$V1%in%x,])))
data_ids<-data_ids[-which(data_ids=="")]
data_ids<-gsub(" ","", data_ids)

OG_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%data_ids[-grep("^OG.*",data_ids)],]
################## 
data<-melt(OG_data)
data$clean_ids<-gsub("-.*","", data$Var1)

data$clean_ids[!data$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"LB"
data$condition<-mapvalues(substr(data$Var2,1,1), 
                          from=c("1","2","3","4","5","6"), 
                          to=c("All","-ST","-LC","-LLm2","-LLm1","-LB"))
data$clean_ids_names<-mapvalues(data$clean_ids, 
                                from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                to=c("LC","LLm2","LLm1","LB","ST"))
#colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')
colors_t0s<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837')
size_all<-25
data$condition<-factor(data$condition,levels = c( "All", "-LLm1","-LC","-ST","-LLm2","-LB"))
data$clean_ids_names<-factor(data$clean_ids_names,levels = c("LC","LLm1","LLm2","LB","ST"))
at9<-ggplot(data,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="Pyruvate carboxylase (PC)")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)+ theme(legend.position = "none")
at9
#ggsave(".../Pyruvate_carboxylase.pdf",at, units="in", width=10, height=6, dpi=600)


#Citrate lysese
# OG6715  OG6716 OG6717
selection_og<-c("OG6715","OG6716","OG6717")
data_ids<-unlist(lapply(selection_og, function(x) as.character(OG_geneid_clean[OG_geneid_clean$V1%in%x,])))
data_ids<-data_ids[-which(data_ids=="")]
data_ids<-gsub(" ","", data_ids)

OG_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%data_ids[-grep("^OG.*",data_ids)],]
####################


#Citrate-sodium symporter
selection_og<-c("OG1323")
data_ids<-unlist(lapply(selection_og, function(x) as.character(OG_geneid_clean[OG_geneid_clean$V1%in%x,])))
data_ids<-data_ids[-which(data_ids=="")]
data_ids<-gsub(" ","", data_ids)

OG_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%data_ids[-grep("^OG.*",data_ids)],]
################## 
data<-melt(OG_data)
data$clean_ids<-gsub("-.*","", data$Var1)

data$clean_ids[!data$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"LB"
data$condition<-mapvalues(substr(data$Var2,1,1), 
                          from=c("1","2","3","4","5","6"), 
                          to=c("All","-ST","-LC","-LLm2","-LLm1","-LB"))
data$clean_ids_names<-mapvalues(data$clean_ids, 
                                from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                to=c("LC","LLm2","LLm1","LB","ST"))
#colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')
colors_t0s<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837')
size_all<-25
data$condition<-factor(data$condition,levels = c( "All", "-LLm1","-LC","-ST","-LLm2","-LB"))
data$clean_ids_names<-factor(data$clean_ids_names,levels = c("LC","LLm1","LLm2","LB","ST"))
at11<-ggplot(data,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="Citrate-sodium symporter")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)+ theme(legend.position = "none")
at11
#ggsave(".../Citrate-sodium_symporter.pdf",at, units="in", width=10, height=6, dpi=600)


############### load manual KO list to connect with metabolomics
library(KEGGREST)
Citrate_M<-c("M00009","M00010","M00012","M00173","M00740")
mod_kegg<-lapply(Citrate_M, function(x) keggGet(x))

h<-1
plots_KOs<-list()
table<-data.frame()
size1<-15
for (i in (1:length(mod_kegg)))
{
  name<-mod_kegg[[i]][[1]]$NAME
  temp_KO_list<-unlist(strsplit(names(mod_kegg[[i]][[1]]$ORTHOLOGY),","))
  for (j in (1:length(temp_KO_list))){
    ids_ko<-emapper.annotations$V1[grep(temp_KO_list[j],emapper.annotations$V12)]
    if (!isEmpty(ids_ko)){
      if(length(ids_ko)==1){
        print(paste("found only once: ",temp_KO_list[j],ids_ko,sep = " "))
      }else{
        KO_RNA_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%ids_ko,]
        #reformat
        data_KO_RNA<-melt(KO_RNA_data)
        data_KO_RNA$clean_ids<-gsub("-.*|_.*","", data_KO_RNA$Var1)
        if (!all((data_KO_RNA$clean_ids== "NLANBGKA") ==F)){
          data_KO_RNA$clean_ids[!data_KO_RNA$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"SICO"
        }else{
          data_KO_RNA$clean_ids[!data_KO_RNA$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675"))]<-"SICO"
        }
        data_KO_RNA$condition<-mapvalues(substr(data_KO_RNA$Var2,1,1), 
                                         from=c("1","2","3","4","5","6"), 
                                         to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
        data_KO_RNA$clean_ids_names<-mapvalues(data_KO_RNA$clean_ids, 
                                               from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                               to=c("cremoris(5614)","lactis(6086)","lactis(10675)","SICO","ST"))
        
        data_KO_RNA$condition<-factor(data_KO_RNA$condition,levels = c( "All", "-lactis(10675)","-cremoris(5614)","-ST","-lactis(6086)","-SICO"))
        plots_KOs[[h]]<- ggplot(data_KO_RNA,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
          labs(x="", y= "mRNA TPM",fill="strains",title=paste(name,temp_KO_list[j],sep = " "))+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
          theme(text = element_text(size=size1),
                plot.title = element_text(size=size1),
                axis.title.x = element_text(size=size1),
                axis.title.y = element_text( size=size1), 
                legend.title = element_text(size=size1),
                axis.text.x = element_text(angle = 45,hjust=1),
                # axis.text.x=element_blank(),
                # axis.ticks.x=element_blank(),
                legend.text = element_text(size=size1))+facet_wrap(~condition)
        
        print(plots_KOs[[h]])
        data_KO_RNA$KO<-rep(temp_KO_list[j],nrow(data_KO_RNA))
        data_KO_RNA$pathway<-rep(name,nrow(data_KO_RNA))
        table<-rbind(table,data_KO_RNA)
        h<-h+1
      }
     
    }else{
      
    }
  }
}

prow <- plot_grid(at11 + theme(legend.position="none"),
                  at2 + theme(legend.position="none"),
                  at9 + theme(legend.position="none"),
                  at1 + theme(legend.position="none"),
                  at11 + theme(legend.position="none"),
                  at3 + theme(legend.position="none") ,
                  labels=c("a", "b","c","d","e","f"), ncol = 2, nrow = 3, label_fontface = "bold",label_size=20,align = "hv")

prow
#ggsave(".../Diacetyl_Acetoin_details.pdf",prow, units="in", width=14, height=18, dpi=600)



