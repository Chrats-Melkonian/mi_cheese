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
################################


################################# st NLANBGKA_00526
# unormalized check  Glutamine synthetase (glnA)
Glutamine<-c("OG1313")
OG_geneid_clean[OG_geneid_clean$V1%in%Glutamine,]
OG1313_ids<-as.character(OG_geneid_clean[OG_geneid_clean$V1%in%Glutamine,])
OG1313_ids<-OG1313_ids[-which(OG1313_ids=="")]
OG1313_ids<-gsub(" ","", OG1313_ids)

OG1313_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%OG1313_ids[-1],]
ST_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00526",]
data_Glutamine<-rbind(OG1313_data,ST_data)
rownames(data_Glutamine)[nrow(data_Glutamine)]<-rownames(htseq_count_minusCHCC3053_TPM)[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00526"]

data_Glutamine<-melt(data_Glutamine)
data_Glutamine$clean_ids<-gsub("-.*|_00526","", data_Glutamine$Var1)
data_Glutamine$clean_ids[!data_Glutamine$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"SICO"
data_Glutamine$condition<-mapvalues(substr(data_Glutamine$Var2,1,1), 
                                    from=c("1","2","3","4","5","6"), 
                                    to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
data_Glutamine$clean_ids_names<-mapvalues(data_Glutamine$clean_ids, 
                                          from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                          to=c("cremoris(5614)","lactis(6086)","lactis(10675)","SICO","ST"))

colors_t0s<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837')
size_all<-20

data_Glutamine$condition<-factor(data_Glutamine$condition,levels = c( "All", "-lactis(10675)","-cremoris(5614)","-ST","-lactis(6086)","-SICO"))
gs<-ggplot(data_Glutamine,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="Glutamine synthetase (glnA)")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)
gs
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/Glutamine_synthetase.pdf",gs, units="in", width=10, height=6, dpi=600)

############## alternative
data_Glutamine2<-data_Glutamine[data_Glutamine$condition%in%c("All","-ST"),]
data_Glutamine2$CN<-factor(paste(data_Glutamine2$condition,data_Glutamine2$clean_ids_names,sep = " "),levels = 
                             c("All cremoris(5614)", "-ST cremoris(5614)",
                                "All lactis(10675)", "-ST lactis(10675)" ,"All lactis(6086)", "-ST lactis(6086)" , 
                                "All SICO", "-ST SICO", "All ST"   ,  "-ST ST" ))
method2 <- "wilcox.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("All cremoris(5614)", "-ST cremoris(5614)"), c("All lactis(10675)", "-ST lactis(10675)"),c("All lactis(6086)", "-ST lactis(6086)")
                       ,c("All SICO", "-ST SICO"),c("All ST"   ,  "-ST ST")) # comparisons for post-hoc tests
#dat$groups<-rownames(temp_data)
g_g <- ggboxplot(data_Glutamine2,
               x = "CN", y = "value",
               color = "clean_ids_names",
               legend = "none",
               palette = colors_t0s,
               add = "jitter", size = 1) + theme(
                 axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y="mRNA TPM", x="strains",title="Glutamine synthetase (glnA)")+theme(text = element_text(size=size_all),
                                                                          plot.title = element_text(size=size_all),
                                                                          axis.title.x = element_text(size=size_all),
                                                                          axis.title.y = element_text( size=size_all), 
                                                                          legend.title = element_text(size=size_all),
                                                                          legend.text = element_text(size=size_all))+
  scale_x_discrete(labels= rep(c("All","-ST"),5))
g_g
gs_stat<-print(
  g_g + 
   stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.signif",size=5) # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
) 
gs_stat
################################# st NLANBGKA_00527
# unormalized check  HTH-type transcriptional regulator GlnR

Glutamine_TR<-c("OG1318")
OG_geneid_clean[OG_geneid_clean$V1%in%Glutamine_TR,]
OG1318_ids<-as.character(OG_geneid_clean[OG_geneid_clean$V1%in%Glutamine_TR,])
OG1318_ids<-OG1318_ids[-which(OG1318_ids=="")]
OG1318_ids<-gsub(" ","", OG1318_ids)

OG1318_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%OG1318_ids[-1],]
ST_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00527",]
data_GlnR<-rbind(OG1318_data,ST_data)
rownames(data_GlnR)[nrow(data_GlnR)]<-rownames(htseq_count_minusCHCC3053_TPM)[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00527"]

data_GlnR<-melt(data_GlnR)
data_GlnR$clean_ids<-gsub("-.*|_00527","", data_GlnR$Var1)
data_GlnR$clean_ids[!data_GlnR$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"SICO"
data_GlnR$condition<-mapvalues(substr(data_GlnR$Var2,1,1), 
                               from=c("1","2","3","4","5","6"), 
                               to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
data_GlnR$clean_ids_names<-mapvalues(data_GlnR$clean_ids, 
                                     from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                     to=c("cremoris(5614)","lactis(6086)","lactis(10675)","SICO","ST"))



size_all<-20

data_GlnR$condition<-factor(data_GlnR$condition,levels = c( "All", "-lactis(10675)","-cremoris(5614)","-ST","-lactis(6086)","-SICO"))
g_tr<-ggplot(data_GlnR,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="Glutamine related transcriptional regulator GlnR")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)
g_tr
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/Glutamine_tr.pdf",g_tr, units="in", width=10, height=6, dpi=600)


############## alternative
data_GlnR2<-data_GlnR[data_GlnR$condition%in%c("All","-ST"),]
data_GlnR2$CN<-factor(paste(data_GlnR2$condition,data_GlnR2$clean_ids_names,sep = " "),levels = 
                             c("All cremoris(5614)", "-ST cremoris(5614)",
                               "All lactis(10675)", "-ST lactis(10675)" ,"All lactis(6086)", "-ST lactis(6086)" , 
                               "All SICO", "-ST SICO", "All ST"   ,  "-ST ST" ))
method2 <- "wilcox.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("All cremoris(5614)", "-ST cremoris(5614)"), c("All lactis(10675)", "-ST lactis(10675)"),c("All lactis(6086)", "-ST lactis(6086)")
                       ,c("All SICO", "-ST SICO"),c("All ST"   ,  "-ST ST")) # comparisons for post-hoc tests
#dat$groups<-rownames(temp_data)
g_tr2 <- ggboxplot(data_GlnR2,
               x = "CN", y = "value",
               color = "clean_ids_names",
               legend = "none",
               palette = colors_t0s,
               add = "jitter", size = 1) + theme(
                 axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y="mRNA TPM", x="strains",title="Glutamine related TR (GlnR)")+theme(text = element_text(size=size_all),
                                                                            plot.title = element_text(size=size_all),
                                                                            axis.title.x = element_text(size=size_all),
                                                                            axis.title.y = element_text( size=size_all), 
                                                                            legend.title = element_text(size=size_all),
                                                                            legend.text = element_text(size=size_all))+
  scale_x_discrete(labels= rep(c("All","-ST"),5))
g_tr2


#########################
legend <- get_legend(g_tr)
prow <- plot_grid(gs  + theme(legend.position="none"),
                  g_tr+ theme(legend.position="none"),
                  legend,
                  labels=c("a", "b",""), ncol = 1, nrow = 3, rel_heights = c(10,10,1),label_fontface = "bold",label_size=20,align = "hv")

prow
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/Glutamine.pdf",prow, units="in", width=10, height=14, dpi=600)

# bcaa_plot bcaaA_plot
# NR_plot
# g_g
# g_tr2
# ammo

library(patchwork)
patchwork_figure3<- ((bcaaA_plot|
               bcaa_plot)
               /(ammo|
                   NR_plot)/(g_g|g_tr2))

patchwork_figure3f<-patchwork_figure3+plot_annotation(tag_levels = 'a')
patchwork_figure3f
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/FINAL_4main/figure3d.pdf",patchwork_figure3f, units="in", width=9, height=14, dpi=600)
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/FINAL_4main/figure1.png", units="in", width=24, height=12, dpi=600)
