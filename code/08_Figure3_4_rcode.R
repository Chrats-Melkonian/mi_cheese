setwd("/media/chrats/Elements/CHEESE_foodtranscriptomics/Projects_FoodTranscriptomics/Cheese/Tier0/Experiment")
library(readxl)
library(plyr)
library(viridis)
library(ggplot2)
acid<-read_excel(path = "190516-004 ACID.xlsx",sheet = "Sheet2")

colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')
acid<-acid[,c(5,7,8)]
colnames(acid)<-gsub(" ","_",acid[1,])
acid<-acid[-1,]
acid_f<-acid[-grep("Milk",acid$Unique_Sample_ID),]
acid_f$samples<-as.integer(unlist(lapply(strsplit(acid_f$Unique_Sample_ID,"-"),function(x) x[1])))
acid_f$rep<-unlist(lapply(strsplit(acid_f$Unique_Sample_ID,"-"),function(x) x[2]))
acid_f$samples_names <- mapvalues(acid_f$samples,
                                  from = c(1,2,3,4,5,6),
                                  to = c("All","-ST","-cremoris","-lactis(6086)","-lactis(10675)","-SICO"))

acid_f$Lactic_acid<-as.numeric(acid_f$Lactic_acid)

library(ggpubr)

dat <- acid_f
# Edit from here
x <- which(names(dat) == "samples_names") # name of grouping variable
y <- which(names(dat) == "Lactic_acid") # names of variables to test
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list( c("All", "-ST"),c("All", "-cremoris"), c("All", "-lactis(6086)"), c("All", "-lactis(10675)"), c("All", "-SICO")) # comparisons for post-hoc tests
# Edit until here
# Edit at your own risk

size_all<-15
p <- ggboxplot(dat,
               x = colnames(dat)[x], y = colnames(dat)[y],
               color = colnames(dat)[x],
               legend = "none",
               palette = colors_t0,
               add = "jitter")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title="",x="Condition",y="Lactic acid (mg/mL)")+ theme(text = element_text(size=size_all),plot.title = element_text(size=size_all),
                                                               axis.title.x = element_text(size=size_all),
                                                               axis.title.y = element_text( size=size_all),
                                                               legend.title = element_text(size=size_all),
                                                               legend.text = element_text(size=size_all))+ guides(colour = guide_legend(override.aes = list(size=4)),shape = guide_legend(override.aes = list(size=4)))

print(
  p #+ stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
  #                     method = method1, label.y = max(dat[, i], na.rm = TRUE)
  #)
  + stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
)


##########################################################
carb<-read_excel(path = "190520-001 CARB.xlsx",sheet = "Sheet2")
carb<-carb[,c(5,7,9)]
colnames(carb)<-gsub(" ","_",colnames(carb))
#carb<-carb[-1,]
carb_f<-carb[-grep("Milk",carb$Unique_Sample_ID),]
carb_f$samples<-as.integer(unlist(lapply(strsplit(carb_f$Unique_Sample_ID,"-"),function(x) x[1])))
carb_f$rep<-unlist(lapply(strsplit(carb_f$Unique_Sample_ID,"-"),function(x) x[2]))
carb_f$samples_names <- mapvalues(carb_f$samples,
                                  from = c(1,2,3,4,5,6),
                                  to = c("All","-ST","-cremoris","-lactis(6086)","-lactis(10675)","-SICO"))
carb_f$Galactose<-as.numeric(carb_f$Galactose)
carb_f[is.na(carb_f)]<-0

#carb_f$Galactose[is.na(carb_f$Galactose)]<-0
size_all<-15

my_comparisons1 <- list( c("All", "-ST"),c("All", "-cremoris"), c("All", "-lactis(6086)"), c("All", "-lactis(10675)"), c("All", "-SICO")) # comparisons for post-hoc tests
g <- ggboxplot(carb_f,
               x = "samples_names", y = "Galactose",
               color = "samples_names",
               legend = "none",
               palette = colors_t0,
               add = "jitter"
) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title="",x="",y="Galactose (mg/mL)")+ theme(text = element_text(size=size_all),plot.title = element_text(size=size_all),
                                                             axis.title.x = element_text(size=size_all),
                                                             axis.title.y = element_text( size=size_all),
                                                             legend.title = element_text(size=size_all),
                                                             legend.text = element_text(size=size_all))+ guides(colour = guide_legend(override.aes = list(size=4)),shape = guide_legend(override.aes = list(size=4)))

print(
  g #+ stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
  #                     method = method1, label.y = max(dat[, i], na.rm = TRUE)
  #)
  + stat_compare_means(comparisons = my_comparisons1, method = method2, label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
)
################Lactose
carb_f$Lactose<-as.numeric(carb_f$Lactose)
carb_f[is.na(carb_f)]<-0
size_all<-15

my_comparisons1 <- list( c("All", "-ST"),c("All", "-cremoris"), c("All", "-lactis(6086)"), c("All", "-lactis(10675)"), c("All", "-SICO")) # comparisons for post-hoc tests
l <- ggboxplot(carb_f,
               x = "samples_names", y = "Lactose",
               color = "samples_names",
               legend = "top",
               palette = colors_t0,
               add = "jitter"
) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title="",x="",y="Lactose (mg/mL)",color="Condition")+ theme(text = element_text(size=size_all),plot.title = element_text(size=size_all),
                                                             axis.title.x = element_text(size=size_all),
                                                             axis.title.y = element_text( size=size_all),
                                                             legend.title = element_text(size=size_all),
                                                             legend.text = element_text(size=size_all))
print(
  l #+ stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
  #                     method = method1, label.y = max(dat[, i], na.rm = TRUE)
  #)
  + stat_compare_means(comparisons = my_comparisons1, method = method2, label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
)
################################
#t1
#################################
setwd("/media/chrats/Elements/CHEESE_foodtranscriptomics/Projects_FoodTranscriptomics/Cheese/Tier1/Biochemical")
data <- read.csv("Cheese Tier 1 biochemical rawdata_v3_14022020.csv")
data_cnames<-colnames(data)
data_cnamesnew<-data[1,]
data<-data[-1,]
colnames(data)<-data_cnamesnew

data_meta<-data[,c(1:4)]

data_03<-data[,16:55] ### metabolites mg/g
col_names<-colnames(data[,16:55])
index<-!apply(data_03, 2, function(x) all(x==0))
data_03<-data_03[,index]
data_03<-t(apply(as.matrix(data_03), 1, as.numeric))
colnames(data_03)<-col_names[index]
data_03[,1:27]<-data_03[,1:27]/1000
##########
Conf_level<-0.95
size_all<-20
h<-1
plot<-list()
for (j in colnames(data_03)){
  data_temp<-data_03[,colnames(data_03)==j]
  data_temp_meta<-cbind(data_temp,data_meta)
  data_temp_meta$batch<-  ifelse(grepl("-3-",data_temp_meta$Sample)==F,"A","B")
  data_temp_meta$`Ripening months`<-factor(data_temp_meta$`Ripening months`, levels = c(unique(data_temp_meta$`Ripening months`)))
  data_temp_meta$`Culture `<-gsub("A3020","All",data_temp_meta$`Culture `)
  df_summary_all<-c()
  for(i in unique(data_temp_meta$Culture)){
    temp<-data_temp_meta[data_temp_meta$Culture==i,]
    df_summary <- data.frame(Time=unique(temp$`Ripening months`), n=tapply(temp$data_temp, temp$`Ripening months`, length),Culture=rep(i,length(unique(temp$`Ripening months`))), mean=tapply(temp$data_temp, temp$`Ripening months`, mean))
    df_summary$sd <- tapply(temp$data_temp, temp$`Ripening months`, sd)
    df_summary$sem <- df_summary$sd/sqrt(df_summary$n-1)
    df_summary$CI_lower <- df_summary$mean + qt((1-Conf_level)/2, df=df_summary$n-1)*df_summary$sem
    df_summary$CI_upper <- df_summary$mean - qt((1-Conf_level)/2, df=df_summary$n-1)*df_summary$sem
    
    df_summary_all<-rbind(df_summary_all,df_summary)
    
  }
  
  plot[[h]] <- ggplot(df_summary_all, aes(x=Time, y=mean)) +
    
    #### plot individual measurements ####
  geom_point(data=data_temp_meta, aes(x=`Ripening months`, y=data_temp,shape=batch,color=`Culture `),size=4) +
    
    #### plot average response over time ####
  geom_line(data=df_summary_all, aes(x=Time, y=mean,group=Culture,color=Culture,linetype=Culture), size=2, alpha=0.8)+
    scale_linetype_manual(name="Condition",values=c(4,3,2,1))+
    
    #### plot error (95%CI) of the response over time ####
  geom_ribbon(data=df_summary_all, aes(ymin=CI_lower, ymax=CI_upper,group=Culture,fill=Culture), alpha=0.2)+
    labs(x="",y=paste(j,"(mg/g)",sep = " "),color="Condition",shape="Batch")+
    theme(text = element_text(size=size_all),
           plot.title = element_text(size=size_all),
           axis.title.x = element_text(size=size_all),
           axis.title.y = element_text( size=size_all), 
           legend.title = element_text(size=size_all),
           legend.text = element_text(size=size_all),
          legend.position="top")+scale_color_manual(values = c('#fc8d62','#8da0cb','#e78ac3','#66c2a5'))+scale_fill_manual(values = c('#fc8d62','#8da0cb','#e78ac3','#66c2a5'))+
    guides(fill=F,linetype=guide_legend(keywidth = 4, keyheight = 1))+scale_shape_manual(values = c(19,8))
  h<-h+1
}

library(cowplot)
pp<-print(
  p 
  + stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 3.8,size=7) )
gg<-print(
  g 
  + stat_compare_means( method = method2, label = "p.signif",ref.group = "All",label.y = 0.5,size=7)+theme(axis.text.x = element_blank()) )

ll<-print(
  l 
  + stat_compare_means( method = method2, label = "p.signif",ref.group = "All",label.y = 36,size=7)+theme(axis.text.x = element_blank()) )
############
lc<-plot[[30]]+labs(x="Time(months)")

gal<-plot[[32]]+theme(axis.text.x = element_blank())

lac<-plot[[34]]+theme(axis.text.x = element_blank())
prow <- plot_grid(plot[[32]] + theme(legend.position="none"),
                  gg + theme(legend.position="none"),
                  plot[[34]] + theme(legend.position="none"),
                  ll+ theme(legend.position="none"),
                  lc + theme(legend.position="none"),
                  
                  pp+ theme(legend.position="none"),
                  labels=c("a", "b","c","d","e","f"), ncol =2 , nrow = 3, label_fontface = "bold",label_size=20,align = "hv")
prow
legend1 <- get_legend(plot[[32]] )
legend2 <- get_legend(ll )
prow2<-plot_grid(legend1,legend2, ncol = 1, nrow = 2)
prow3<-plot_grid(prow,prow2, rel_widths=c(12,1),rel_heights=c(12,1), ncol = 1, nrow = 2)
prow3
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/figure3.pdf",prow3, units="in", width=11, height=14, dpi=600)
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/figure3.png",prow3, units="in", width=11, height=14, dpi=600)
plot_grid(prow2)
###########3 connect t1-t0 significant 

############################
#next plot

################# MS data PCA
ms_data<-read_excel(path = "/home/chrats/Desktop/Git/cheese-ft/Data/T0/Biochemical/190524-006 HA VOC MS1 results.xlsx",sheet = "Sheet2")
ms_data_r<-ms_data[,c(7,9:ncol(ms_data))]


empty_mean<-apply(ms_data_r[20:22,-1], 2, function(x) mean(x, na.rm=TRUE))
empty_mean[3]<-0
ms_data_r[,-1]-empty_mean

ms_data_r[,-1]<- sweep(ms_data_r[,-1], 2, empty_mean, "-")
ms_data_r[,-1][ms_data_r[,-1]<0] <- 0

plot(apply(ms_data_r[-c(20:22),-1], 2,function(x) sd(x,na.rm = T)))

rm_names<-names(which(apply(ms_data_r[-c(20:22),-1], 2,function(x) sd(x,na.rm = T))<10))
ms_data_f<-ms_data_r[,!colnames(ms_data_r)%in%rm_names]
ms_data_f<-ms_data_f[-c(20:22),]
ms_data_f[is.na(ms_data_f)]<-0

library(ggbiplot)
ms_data_f.pca <- prcomp(ms_data_f[,-1], center = TRUE,scale. = TRUE)

groups_pca<-ms_data_f$`Unique Sample ID`
groups_pca<-unlist(lapply(strsplit(groups_pca,"-"),function(x) x[1]))
carb_f$rep<-unlist(lapply(strsplit(carb_f$Unique_Sample_ID,"-"),function(x) x[2]))
groups_pca <- mapvalues(groups_pca,
                        from = c("1","2","3","4","5","6","milk"),
                        to = c("All","minus_ST","minus_5614_c","minus_6086_l","minus_10675_l","minus_SICO","Milk"))


color<-c("grey70",viridis(6))
ggbiplot(ms_data_f.pca,ellipse=TRUE, groups=groups_pca,var.axes=F)+ geom_point(size=3,aes(color=groups_pca))+ xlim(-2.5, 2)+ theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_manual(name="Conditions", values= color)+ theme(text = element_text(size=size_all),plot.title = element_text(size=size_all),
                                                               axis.title.x = element_text(size=size_all),
                                                               axis.title.y = element_text( size=size_all),
                                                               legend.title = element_text(size=size_all),
                                                               legend.text = element_text(size=size_all))+ guides(colour = guide_legend(override.aes = list(size=4)),shape = guide_legend(override.aes = list(size=4)))


ggbiplot(ms_data_f.pca,ellipse=TRUE,choices = c(3,4), groups=groups_pca)+ geom_point(size=3,aes(color=groups_pca))+ xlim(-2.5, 2)+ theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_manual(name="Conditions", values= color)+ theme(text = element_text(size=size_all),plot.title = element_text(size=size_all),
                                                               axis.title.x = element_text(size=size_all),
                                                               axis.title.y = element_text( size=size_all),
                                                               legend.title = element_text(size=size_all),
                                                               legend.text = element_text(size=size_all))+ guides(colour = guide_legend(override.aes = list(size=4)),shape = guide_legend(override.aes = list(size=4)))




#####################################
################# carb, acids, MS
data_all_meta<-cbind(ms_data_f,carb[-3,c("Galactose","Lactose")],acid[-3,2:3])
data_all_meta.r<-data_all_meta[,-1]
data_all_meta.r[data_all_meta.r=="<LOD"] <-0
data<-apply(as.matrix(data_all_meta.r), 2, as.numeric)
data_all_meta.pca <- prcomp(data, center = TRUE,scale. = TRUE)

groups_pca<-data_all_meta$`Unique Sample ID`
groups_pca<-unlist(lapply(strsplit(groups_pca,"-"),function(x) x[1]))
# data_all_meta$rep<-unlist(lapply(strsplit(data_all_meta$Unique_Sample_ID,"-"),function(x) x[2]))
groups_pca <- mapvalues(groups_pca,
                        from = c("1","2","3","4","5","6","milk"),
                        to = c("All","-ST","-cremoris","-lactis(6086)","-lactis(10675)","-SICO","Milk"))


color<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837','#a50026',"black")
#color<-c('#e41a1c',"grey70",'#4daf4a','#984ea3','#ff7f00','#a65628','#377eb8')

top5<-names(sort(data_all_meta.pca$scale,decreasing = T)[1:5])
#data_all_meta.pca$scale<-data_all_meta.pca$scale[names(data_all_meta.pca$scale)%in%top5]
data_all_meta.pca$rotation<-data_all_meta.pca$rotation[rownames(data_all_meta.pca$rotation)%in%top5,]

pca_metab<-ggbiplot(data_all_meta.pca,ellipse=T, groups=groups_pca,var.axes=T,varname.adjust = 1.1,varname.size = 5)+ geom_point(size=4,aes(color=groups_pca))+ theme_minimal()+
  theme(legend.position = "right")+xlim(-2,3)+ylim(-2,3)+labs(x="PC1 (33.0% exp.var.)",y="PC2 (24.2% exp.var.)")+
  scale_colour_manual(name="Conditions", values= color)+ theme(text = element_text(size=size_all),plot.title = element_text(size=size_all),
                                                               axis.title.x = element_text(size=size_all),
                                                               axis.title.y = element_text( size=size_all),
                                                               legend.title = element_text(size=size_all),
                                                               legend.text = element_text(size=size_all))+ guides(colour = guide_legend(override.aes = list(size=4)),shape = guide_legend(override.aes = list(size=4)))



#carb_f$Galactose[is.na(carb_f$Galactose)]<-0
size_all<-20

ms_data_f_nomilk<-ms_data_f[-c(18,19),]
groups_pca1<-ms_data_f_nomilk$`Unique Sample ID`
groups_pca1<-unlist(lapply(strsplit(groups_pca1,"-"),function(x) x[1]))
ms_data_f_nomilk$rep<-unlist(lapply(strsplit(groups_pca1,"-"),function(x) x[2]))
ms_data_f_nomilk$groups <- mapvalues(groups_pca1,
                                     from = c("1","2","3","4","5","6"),
                                     to =  c("All","-ST","-cremoris","-lactis(6086)","-lactis(10675)","-SICO"))

dat <- ms_data_f_nomilk
colnames(dat)[2:22]<-gsub("/","",colnames(dat[2:22]))
colnames(dat)[2:22]<-gsub(" ","_",colnames(dat[2:22]))
# Edit from here
x <- which(names(dat) == "groups") # name of grouping variable
y <- colnames(dat)[2:22]
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons1 <- list( c("All", "-ST"),c("All", "-cremoris"), c("All", "-lactis(6086)"), c("All", "-lactis(10675)"), c("All", "-SICO")) # comparisons for post-hoc tests
# Edit until here
# Edit at your own risk

milk_data<-ms_data_f[c(18,19),2:22]
milk_data_mean<-apply(milk_data, 2, mean)

j<-1

plot_metMS<-list()
for (i in y) {
  plot_metMS[[j]] <- ggboxplot(dat,
                 x = "groups", y = i,
                 color = "groups",
                 legend = "none",
                 palette = colors_t0,
                 add = "jitter"
  ) + theme(text = element_text(size=size_all),plot.title = element_text(size=size_all),
            axis.title.x = element_text(size=size_all),
            axis.title.y = element_text( size=size_all),
            legend.title = element_text(size=size_all),
            legend.text = element_text(size=size_all),
    axis.text.x = element_text(angle = 45, hjust = 1))+ geom_hline(yintercept=milk_data_mean[j])
  j<-j+1
}

dia<-plot_metMS[[4]]+stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 4000,size=7)+theme(axis.title.x=element_blank())
P23<-plot_metMS[[8]]+stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 1500,size=7)+theme(axis.title.x=element_blank(),
                                                                                                                              axis.text.x=element_blank())
ace<-plot_metMS[[10]]+stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 4000,size=7)+theme(axis.title.x=element_blank())
legend2 <- get_legend(pca_metab)
prow4 <- plot_grid(pca_metab + theme(legend.position="none"),
                   # P23+ theme(legend.position="none"),
                  dia + theme(legend.position="none"),
                  ace + theme(legend.position="none"),
                  labels=c("a", "b","c"), ncol = 3, nrow = 1, label_fontface = "bold",label_size=30,align = "hv")


prow4
prow5<-plot_grid(prow4,legend2,ncol = 2, nrow = 1,rel_widths = c(4,1),rel_heights = c(4,1))
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/figure3.5.2.pdf",prow5, units="in", width=25, height=8, dpi=600)

#### all compounds relavent to cremoris
ethyl_h<-plot_metMS[[18]]+theme(axis.title.x=element_blank(),
                                axis.text.x=element_blank(),
                                axis.ticks.x=element_blank()) +stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 300,size=7)+theme(axis.title.x=element_blank())
 # cremoris produce 
ethylo_a<-plot_metMS[[6]]+theme(axis.title.x=element_blank(),
                                axis.text.x=element_blank(),
                                axis.ticks.x=element_blank()) +stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 400,size=7)+theme(axis.title.x=element_blank())


heptanone<-plot_metMS[[14]] +stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 40,size=7)+theme(axis.title.x=element_blank())
 # lower, inc. ST
methyl_thiolanone<-plot_metMS[[19]] +stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 50,size=7)+theme(axis.title.x=element_blank())
 #ST +cremoris
prow7 <- plot_grid(ethyl_h + theme(legend.position="none"),
                   ethylo_a + theme(legend.position="none"),
                   heptanone + theme(legend.position="none"),
                   methyl_thiolanone+theme(legend.position="none") ,
                   labels=c("a", "b","c","d"), ncol = 2, nrow = 2, label_fontface = "bold",label_size=30,align = "hv")
prow7<-plot_grid(prow7,legend2,ncol = 2, nrow = 1,rel_widths = c(4,1),rel_heights = c(4,1))
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/figure6_cremoris_st.pdf",prow7, units="in", width=11, height=11, dpi=600)


plot_metMS[[15]] #  - cremoris, high
plot_metMS[[13]] #  - cremoris, high
#plot_metMS[[10]] #  - cremoris, high main plot
plot_metMS[[9]] #  - cremoris, high
plot_metMS[[8]] 

heptanal<-plot_metMS[[15]] +stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 35,size=7)+theme(axis.title.x=element_blank())
hexanal<-plot_metMS[[13]]+stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 1000,size=7)+theme(axis.title.x=element_blank())
ethyl_furan<-plot_metMS[[9]]+stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 35,size=7)+theme(axis.title.x=element_blank())
pentanedione<-plot_metMS[[8]] +stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 1000,size=7)+theme(axis.title.x=element_blank())
prow8 <- plot_grid(heptanal + theme(legend.position="none")+theme(axis.title.x=element_blank(),
                                                                  axis.text.x=element_blank(),
                                                                  axis.ticks.x=element_blank()),
                   hexanal + theme(legend.position="none")+theme(axis.title.x=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank()),
                   ethyl_furan + theme(legend.position="none"),
                   pentanedione+theme(legend.position="none") ,
                   labels=c("a", "b","c","d"), ncol = 2, nrow = 2, label_fontface = "bold",label_size=30,align = "hv")
prow8<-plot_grid(prow8,legend2,ncol = 2, nrow = 1,rel_widths = c(4,1),rel_heights = c(4,1))
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/figure6_cremoris_only.pdf",prow8, units="in", width=11, height=11, dpi=600)


#plot_metMS[[4]]  main plot
#### all compounds relavent to SICO 10675 
diaethyl_dis<-plot_metMS[[11]]+stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 80,size=7)+theme(axis.title.x=element_blank())
pentanone<-plot_metMS[[7]] +theme(axis.title.x=element_blank(),
                                  axis.text.x=element_blank(),
                                  axis.ticks.x=element_blank())+stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 190,size=7)+theme(axis.title.x=element_blank())
butanone<-plot_metMS[[5]]+theme(axis.title.x=element_blank(),
                                axis.text.x=element_blank(),
                                axis.ticks.x=element_blank()) +stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 150,size=7)+theme(axis.title.x=element_blank())
acetone<-plot_metMS[[2]] +stat_compare_means( method = method2, label = "p.signif", ref.group = "All",label.y = 150,size=7)+theme(axis.title.x=element_blank())

prow6 <- plot_grid(pentanone + theme(legend.position="none"),
                   butanone + theme(legend.position="none"),
                   acetone + theme(legend.position="none"),
                   diaethyl_dis+theme(legend.position="none") ,
                   labels=c("a", "b","c","d"), ncol = 2, nrow = 2, label_fontface = "bold",label_size=30,align = "hv")
prow6<-plot_grid(prow6,legend2,ncol = 2, nrow = 1,rel_widths = c(4,1),rel_heights = c(4,1))
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/figure6_S_10675_SICO.pdf",prow6, units="in", width=11, height=11, dpi=600)
### all ST


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


###########################
########################### ST NLANBGKA_00232 NrgA
# unormalized check Ammonium transporter  (nrgA/amtB)
Amo_transporter_og<-"OG1286"
OG_geneid_clean[OG_geneid_clean$V1%in%Amo_transporter_og,]
OG1286_ids<-as.character(OG_geneid_clean[OG_geneid_clean$V1%in%Amo_transporter_og,])
OG1286_ids<-OG1286_ids[-which(OG1286_ids=="")]
OG1286_ids<-gsub(" ","", OG1286_ids)

OG1286_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%OG1286_ids[-1],]
ST_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00232",]
data_ammonium<-rbind(OG1286_data,ST_data)
rownames(data_ammonium)[nrow(data_ammonium)]<-rownames(htseq_count_minusCHCC3053_TPM)[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00232"]
library(reshape2)
library(ggplot2)
library(plyr)
################## 
data<-melt(data_ammonium)
data$clean_ids<-gsub("-.*|_00232","", data$Var1)

data$clean_ids[!data$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"SICO"
data$condition<-mapvalues(substr(data$Var2,1,1), 
                          from=c("1","2","3","4","5","6"), 
                          to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
data$clean_ids_names<-mapvalues(data$clean_ids, 
                                from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                to=c("cremoris(5614)","lactis(6086)","lactis(10675)","SICO","ST"))
#colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')
colors_t0s<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837')
size_all<-20

data$condition<-factor(data$condition,levels = c( "All", "-lactis(10675)","-cremoris(5614)","-ST","-lactis(6086)","-SICO"))
at<-ggplot(data,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="Ammonium transporter (nrgA/amtB)")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)
at
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/Ammonium_transporter.pdf",at, units="in", width=10, height=6, dpi=600)


############## alternative
ammonium_data2<-data[data$condition%in%c("All","-ST"),]
ammonium_data2$CN<-factor(paste(ammonium_data2$condition,ammonium_data2$clean_ids_names,sep = " "),levels = 
                        c("All cremoris(5614)", "-ST cremoris(5614)",
                          "All lactis(10675)", "-ST lactis(10675)" ,"All lactis(6086)", "-ST lactis(6086)" , 
                          "All SICO", "-ST SICO", "All ST"   ,  "-ST ST" ))
method2 <- "wilcox.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("All cremoris(5614)", "-ST cremoris(5614)"), c("All lactis(10675)", "-ST lactis(10675)"),c("All lactis(6086)", "-ST lactis(6086)")
                       ,c("All SICO", "-ST SICO"),c("All ST"   ,  "-ST ST")) # comparisons for post-hoc tests
#dat$groups<-rownames(temp_data)
ammo <- ggboxplot(ammonium_data2,
                   x = "CN", y = "value",
                   color = "clean_ids_names",
                   legend = "none",
                   palette = colors_t0s,
                   add = "jitter", size = 1) + theme(
                     axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y="mRNA TPM", x="strains",title="Ammonium T (nrgA/amtB)")+theme(text = element_text(size=size_all),
                                                                                                 plot.title = element_text(size=size_all),
                                                                                                 axis.title.x = element_text(size=size_all),
                                                                                                 axis.title.y = element_text( size=size_all), 
                                                                                                 legend.title = element_text(size=size_all),
                                                                                                 legend.text = element_text(size=size_all))+
  scale_x_discrete(labels= rep(c("All","-ST"),5))
ammo

################################### LL
# Branched-chain amino acid transport system 2 carrier protein: OG4556 OG4597 OG5376

################################### ST
# Branched-chain amino acid transport system 2 carrier protein: NLANBGKA_00561


BCAA<-c("OG4556")
OG_geneid_clean[OG_geneid_clean$V1%in%BCAA,]
BCAA_ids<-as.character(OG_geneid_clean[OG_geneid_clean$V1%in%BCAA,])
BCAA_ids<-BCAA_ids[-which(BCAA_ids=="")]
BCAA_ids<-gsub(" ","", BCAA_ids)
BCAA_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%BCAA_ids[-1],]

BCAA2<-c("OG4597")
OG_geneid_clean[OG_geneid_clean$V1%in%BCAA2,]
BCAA_ids2<-as.character(OG_geneid_clean[OG_geneid_clean$V1%in%BCAA2,])
BCAA_ids2<-BCAA_ids2[-which(BCAA_ids2=="")]
BCAA_ids2<-gsub(" ","", BCAA_ids2)
BCAA_data2<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%BCAA_ids2[-1],]

# BCAA3<-c("OG5376")
# OG_geneid_clean[OG_geneid_clean$V1%in%BCAA3,]
# BCAA_ids3<-as.character(OG_geneid_clean[OG_geneid_clean$V1%in%BCAA3,])
# BCAA_ids3<-BCAA_ids3[-which(BCAA_ids3=="")]
# BCAA_ids3<-gsub(" ","", BCAA_ids3)
# BCAA_data3<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%BCAA_ids3[-1],]

ST_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00561",]
data_BCAA<-rbind(BCAA_data,BCAA_data2,ST_data)
rownames(data_BCAA)[nrow(data_BCAA)]<-rownames(htseq_count_minusCHCC3053_TPM)[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00561"]

data_BCAA<-melt(data_BCAA)
data_BCAA$clean_ids<-gsub("-.*|_00561","", data_BCAA$Var1)
data_BCAA$clean_ids[!data_BCAA$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"SICO"
data_BCAA$condition<-mapvalues(substr(data_BCAA$Var2,1,1), 
                                    from=c("1","2","3","4","5","6"), 
                                    to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
data_BCAA$clean_ids_names<-mapvalues(data_BCAA$clean_ids, 
                                          from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                          to=c("cremoris(5614)","lactis(6086)","lactis(10675)","SICO","ST"))

colors_t0s<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837')
size_all<-20
# ggplot(data_OG1313,aes(x=condition,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
#   labs(x="", y= "mRNA count",fill="strains", title="Glutamine synthetase (glnA)")+#geom_jitter(color="black", size=0.4, alpha=0.9,aes(fill=clean_ids_names))+geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
#   theme(text = element_text(size=size_all),
#         plot.title = element_text(size=size_all),
#         axis.title.x = element_text(size=size_all),
#         axis.title.y = element_text( size=size_all), 
#         legend.title = element_text(size=size_all),
#         axis.text.x = element_text(angle = 45,hjust=1),
#         legend.text = element_text(size=size_all))
data_BCAA$condition<-factor(data_BCAA$condition,levels = c( "All", "-lactis(10675)","-cremoris(5614)","-ST","-lactis(6086)","-SICO"))
bcaa<-ggplot(data_BCAA,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="Branched-chain amino acid transport system 2 carrier protein")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)
bcaa


############## alternative
bcaa_data2<-data_BCAA[data_BCAA$condition%in%c("All","-ST"),]
bcaa_data2$CN<-factor(paste(bcaa_data2$condition,bcaa_data2$clean_ids_names,sep = " "),levels = 
                            c("All cremoris(5614)", "-ST cremoris(5614)",
                              "All lactis(10675)", "-ST lactis(10675)" ,"All lactis(6086)", "-ST lactis(6086)" , 
                              "All SICO", "-ST SICO", "All ST"   ,  "-ST ST" ))
method2 <- "wilcox.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("All cremoris(5614)", "-ST cremoris(5614)"), c("All lactis(10675)", "-ST lactis(10675)"),c("All lactis(6086)", "-ST lactis(6086)")
                       ,c("All SICO", "-ST SICO"),c("All ST"   ,  "-ST ST")) # comparisons for post-hoc tests
#dat$groups<-rownames(temp_data)
bcaa_plot <- ggboxplot(bcaa_data2,
                  x = "CN", y = "value",
                  color = "clean_ids_names",
                  legend = "none",
                  palette = colors_t0s,
                  add = "jitter", size = 1) + theme(
                    axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y="mRNA TPM", x="strains",title="BCAA transport system 2 carrier protein")+theme(text = element_text(size=size_all),
                                                                                                 plot.title = element_text(size=size_all),
                                                                                                 axis.title.x = element_text(size=size_all),
                                                                                                 axis.title.y = element_text( size=size_all), 
                                                                                                 legend.title = element_text(size=size_all),
                                                                                                 legend.text = element_text(size=size_all))+
  scale_x_discrete(labels= rep(c("All","-ST"),5))
bcaa_plot

##########################3


# Branched-chain-amino-acid aminotransferase OG1572
# Branched-chain-amino-acid aminotransferase NLANBGKA_00875
BCAAa<-c("OG1572")
OG_geneid_clean[OG_geneid_clean$V1%in%BCAAa,]
BCAAa_ids<-as.character(OG_geneid_clean[OG_geneid_clean$V1%in%BCAAa,])
BCAAa_ids<-BCAAa_ids[-which(BCAAa_ids=="")]
BCAAa_ids<-gsub(" ","", BCAAa_ids)
BCAAa_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%BCAAa_ids[-1],]


ST_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00875",]
data_BCAAa<-rbind(BCAAa_data,ST_data)
rownames(data_BCAAa)[nrow(data_BCAAa)]<-rownames(htseq_count_minusCHCC3053_TPM)[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00875"]

data_BCAAa<-melt(data_BCAAa)
data_BCAAa$clean_ids<-gsub("-.*|_00875","", data_BCAAa$Var1)
data_BCAAa$clean_ids[!data_BCAAa$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"SICO"
data_BCAAa$condition<-mapvalues(substr(data_BCAAa$Var2,1,1), 
                               from=c("1","2","3","4","5","6"), 
                               to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
data_BCAAa$clean_ids_names<-mapvalues(data_BCAAa$clean_ids, 
                                     from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                     to=c("cremoris(5614)","lactis(6086)","lactis(10675)","SICO","ST"))

colors_t0s<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837')
size_all<-20

data_BCAAa$condition<-factor(data_BCAAa$condition,levels = c( "All", "-lactis(10675)","-cremoris(5614)","-ST","-lactis(6086)","-SICO"))
bcaaA<-ggplot(data_BCAAa,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="Branched-chain-amino-acid aminotransferase")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)
bcaaA

############## alternative
bcaaA_data2<-data_BCAAa[data_BCAAa$condition%in%c("All","-ST"),]
bcaaA_data2$CN<-factor(paste(bcaaA_data2$condition,bcaaA_data2$clean_ids_names,sep = " "),levels = 
                        c("All cremoris(5614)", "-ST cremoris(5614)",
                          "All lactis(10675)", "-ST lactis(10675)" ,"All lactis(6086)", "-ST lactis(6086)" , 
                          "All SICO", "-ST SICO", "All ST"   ,  "-ST ST" ))
method2 <- "wilcox.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("All cremoris(5614)", "-ST cremoris(5614)"), c("All lactis(10675)", "-ST lactis(10675)"),c("All lactis(6086)", "-ST lactis(6086)")
                       ,c("All SICO", "-ST SICO"),c("All ST"   ,  "-ST ST")) # comparisons for post-hoc tests
#dat$groups<-rownames(temp_data)
bcaaA_plot <- ggboxplot(bcaaA_data2,
                       x = "CN", y = "value",
                       color = "clean_ids_names",
                       legend = "none",
                       palette = colors_t0s,
                       add = "jitter", size = 1) + theme(
                         axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y="mRNA TPM", x="strains",title="BCAA aminotransferase")+theme(text = element_text(size=size_all),
                                                                                                             plot.title = element_text(size=size_all),
                                                                                                             axis.title.x = element_text(size=size_all),
                                                                                                             axis.title.y = element_text( size=size_all), 
                                                                                                             legend.title = element_text(size=size_all),
                                                                                                             legend.text = element_text(size=size_all))+
  scale_x_discrete(labels= rep(c("All","-ST"),5))
bcaaA_plot
##########################################
legend1 <- get_legend(bcaaA+theme(legend.position="bottom",legend.title = element_text(size=15),
                                  axis.text.x = element_text(angle = 45,hjust=1),
                                  legend.text = element_text(size=15)))
prow <- plot_grid(bcaaA+ theme(legend.position="none")+theme(axis.title.x=element_blank(),
                                                             axis.text.x=element_blank(),
                                                             axis.ticks.x=element_blank()),
                  bcaa + theme(legend.position="none")+theme(axis.title.x=element_blank(),
                                                             axis.text.x=element_blank(),
                                                             axis.ticks.x=element_blank()),
                  at  + theme(legend.position="none"),
                  labels=c("a", "b","c"), ncol = 1, nrow = 3, label_fontface = "bold",label_size=25,align = "hv")

prow2<-plot_grid(legend1,prow, rel_heights  = c( .2,3), ncol = 1, nrow = 2)
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/BCAA_AT.pdf",prow2, units="in", width=8, height=16, dpi=600)



###################
# Nitrogen regulatory protein P-II OG1285
# Nitrogen regulatory protein P-II NLANBGKA_00233
BCAAa<-c("OG1285")
OG_geneid_clean[OG_geneid_clean$V1%in%BCAAa,]
BCAAa_ids<-as.character(OG_geneid_clean[OG_geneid_clean$V1%in%BCAAa,])
BCAAa_ids<-BCAAa_ids[-which(BCAAa_ids=="")]
BCAAa_ids<-gsub(" ","", BCAAa_ids)
BCAAa_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%BCAAa_ids[-1],]


ST_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00233",]
data_BCAAa<-rbind(BCAAa_data,ST_data)
rownames(data_BCAAa)[nrow(data_BCAAa)]<-rownames(htseq_count_minusCHCC3053_TPM)[rownames(htseq_count_minusCHCC3053_TPM)%in%"NLANBGKA_00233"]

data_BCAAa<-melt(data_BCAAa)
data_BCAAa$clean_ids<-gsub("-.*|_00233","", data_BCAAa$Var1)
data_BCAAa$clean_ids[!data_BCAAa$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675","NLANBGKA"))]<-"SICO"
data_BCAAa$condition<-mapvalues(substr(data_BCAAa$Var2,1,1), 
                                from=c("1","2","3","4","5","6"), 
                                to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
data_BCAAa$clean_ids_names<-mapvalues(data_BCAAa$clean_ids, 
                                      from=c("CHCC5614","CHCC6086","CHCC10675","SICO","NLANBGKA"), 
                                      to=c("cremoris(5614)","lactis(6086)","lactis(10675)","SICO","ST"))

colors_t0s<-c('#fdae61','#313695','#4575b4','#abd9e9','#1b7837')
size_all<-20

data_BCAAa$condition<-factor(data_BCAAa$condition,levels = c( "All", "-lactis(10675)","-cremoris(5614)","-ST","-lactis(6086)","-SICO"))
nrp2<-ggplot(data_BCAAa,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="Nitrogen regulatory protein P-II")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)
nrp2
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/nrp2.pdf",prow2, units="in", width=8, height=16, dpi=600)
###########33 alternative
NR<-data_BCAAa[data_BCAAa$condition%in%c("All","-ST"),]
NR$CN<-factor(paste(NR$condition,NR$clean_ids_names,sep = " "),levels = 
                        c("All cremoris(5614)", "-ST cremoris(5614)",
                          "All lactis(10675)", "-ST lactis(10675)" ,"All lactis(6086)", "-ST lactis(6086)" , 
                          "All SICO", "-ST SICO", "All ST"   ,  "-ST ST" ))
method2 <- "wilcox.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("All cremoris(5614)", "-ST cremoris(5614)"), c("All lactis(10675)", "-ST lactis(10675)"),c("All lactis(6086)", "-ST lactis(6086)")
                       ,c("All SICO", "-ST SICO"),c("All ST"   ,  "-ST ST")) # comparisons for post-hoc tests
#dat$groups<-rownames(temp_data)
NR_plot <- ggboxplot(NR,
                       x = "CN", y = "value",
                       color = "clean_ids_names",
                       legend = "none",
                       palette = colors_t0s,
                       add = "jitter", size = 1) + theme(
                         axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y="mRNA TPM", x="strains",title="Nitrogen regulatory protein P-II")+theme(text = element_text(size=size_all),
                                                                                                             plot.title = element_text(size=size_all),
                                                                                                             axis.title.x = element_text(size=size_all),
                                                                                                             axis.title.y = element_text( size=size_all), 
                                                                                                             legend.title = element_text(size=size_all),
                                                                                                             legend.text = element_text(size=size_all))+
  scale_x_discrete(labels= rep(c("All","-ST"),5))
NR_plot

##########################################
legend1 <- get_legend(bcaaA+theme(legend.position="bottom",legend.title = element_text(size=15),
                                  axis.text.x = element_text(angle = 45,hjust=1),
                                  legend.text = element_text(size=15)))
prow <- plot_grid(bcaaA+ theme(legend.position="none")+theme(axis.title.x=element_blank(),
                                                             axis.text.x=element_blank(),
                                                             axis.ticks.x=element_blank()),
                  bcaa + theme(legend.position="none")+theme(axis.title.x=element_blank(),
                                                             axis.text.x=element_blank(),
                                                             axis.ticks.x=element_blank()),
                  at  + theme(legend.position="none"),
                  labels=c("a", "b","c"), ncol = 1, nrow = 3, label_fontface = "bold",label_size=25,align = "hv")

prow2<-plot_grid(legend1,prow, rel_heights  = c( .2,3), ncol = 1, nrow = 2)
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/BCAA_AT.pdf",prow2, units="in", width=8, height=16, dpi=600)


#################################
# unormalized check  3-phenylpropionate-dihydrodiol/cinnamic acid-dihydrodiol dehydrogenase
colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')
phenylpropionate<-c("OG2576")
OG_geneid_clean[OG_geneid_clean$V1%in%phenylpropionate,]
OG2576_ids<-as.character(OG_geneid_clean[OG_geneid_clean$V1%in%phenylpropionate,])
OG2576_ids<-OG2576_ids[-which(OG2576_ids=="")]
OG2576_ids<-gsub(" ","", OG2576_ids)

OG2576_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%OG2576_ids[-1],]

data_OG2576<-melt(OG2576_data)
data_OG2576$clean_ids<-gsub("-.*","", data_OG2576$Var1)
data_OG2576$clean_ids[!data_OG2576$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675"))]<-"SICO"
data_OG2576$condition<-mapvalues(substr(data_OG2576$Var2,1,1), 
                                 from=c("1","2","3","4","5","6"), 
                                 to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
data_OG2576$clean_ids_names<-mapvalues(data_OG2576$clean_ids, 
                                       from=c("CHCC5614","CHCC6086","CHCC10675","SICO"), 
                                       to=c("cremoris(5614)","lactis(6086)","lactis(10675)","SICO"))

colors_t0s2<-c('#4575b4','#abd9e9')
size_all<-20
data_OG2576$condition<-factor(data_OG2576$condition,levels = c( "All", "-lactis(10675)","-cremoris(5614)","-ST","-lactis(6086)","-SICO"))
phe<-ggplot(data_OG2576,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
  labs(x="", y= "mRNA TPM",fill="strains",title="3-phenylpropionate-dihydrodiol/cinnamic acid-dihydrodiol dehydrogenase")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=15),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        axis.text.x = element_text(angle = 45,hjust=1),
        legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s2)
phe
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/phenylpropionate.pdf",phe, units="in", width=10, height=6, dpi=600)

############### load manual KO list to connect with metabolomics
table_metabolites_connect_KEGG <- read.csv("/media/chrats/Elements/CHEESE_foodtranscriptomics/Projects_FoodTranscriptomics/Cheese/Tier0/connect_metabolomics_metatranscriptomics/table_metabolites_connect_KEGG.csv")
emapper.annotations <- read.delim("~/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/3_count_tables/emapper.annotations.tsv", header=FALSE, comment.char="#")
emapper.annotations$V1<-unlist(lapply(strsplit(emapper.annotations$V1,"ID=|;"), function(x) x[2]))
emapper.annotations$V12<-gsub("ko:","",emapper.annotations$V12)
table_with_KOs<-table_metabolites_connect_KEGG[!is.na(table_metabolites_connect_KEGG$Kos),]

unique(unlist(strsplit(table_with_KOs$Kos,",")))
unique(unlist(strsplit(table_with_KOs$Kos,",")))%in%emapper.annotations$V12

h<-1
plots_KOs<-list()
size1<-15
for (i in 1:nrow(table_with_KOs)){
  temp<-table_with_KOs[i,]
  ALL_KOs<-unique(unlist(strsplit(temp$Kos,",")))
  for (j in (1:length(ALL_KOs))){
    ids_ko<-emapper.annotations$V1[grep(ALL_KOs[j],emapper.annotations$V12)]
    if (!isEmpty(ids_ko)){
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
        labs(x="", y= "mRNA TPM",fill="strains",title=ALL_KOs[j])+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
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
      h<-h+1
    }else{
      
    }
    
  }
}

prow <- plot_grid(plots_KOs[[1]]  + theme(legend.position="none"),
                  plots_KOs[[2]] + theme(legend.position="none"),
                  plots_KOs[[3]] + theme(legend.position="none"),
                  plots_KOs[[4]] + theme(legend.position="none") ,
                  plots_KOs[[5]]  + theme(legend.position="none"),
                  plots_KOs[[6]] + theme(legend.position="none") ,
                  labels=c("a", "b","c","d","e","f"), ncol = 2, nrow = 3, label_fontface = "bold",label_size=20,align = "hv")

prow
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/SUPP/Diacetyl_Acetoin.pdf",prow, units="in", width=10, height=14, dpi=600)

legend1 <- get_legend(plot[[1]])
prow2<-plot_grid(legend1,prow, rel_heights  = c( .2,3), ncol = 1, nrow = 2)


ALL_KOs<-unique(unlist(strsplit(table_metabolites_connect_KEGG$Kos,",")))
htseq_count_minusCHCC3053_TMP

library(reshape2)
library(ggplot2)
library(plyr)
#Citrate-sodium symporter OG1323

  Citrate_og<-"OG1323"
  OG_geneid_clean[OG_geneid_clean$V1%in%Citrate_og,]
  Citrate_og_ids<-as.character(OG_geneid_clean[OG_geneid_clean$V1%in%Citrate_og,])
  Citrate_og_ids<-Citrate_og_ids[-which(Citrate_og_ids=="")]
  Citrate_og_ids<-gsub(" ","", Citrate_og_ids)
  
  Citrate_data<-htseq_count_minusCHCC3053_TPM[rownames(htseq_count_minusCHCC3053_TPM)%in%Citrate_og_ids[-1],]

  ################## 
  data<-melt(Citrate_data)
  data$clean_ids<-gsub("-.*","", data$Var1)
  
  data$clean_ids[!data$clean_ids%in%(c("CHCC5614","CHCC6086","CHCC10675"))]<-"SICO"
  data$condition<-mapvalues(substr(data$Var2,1,1), 
                            from=c("1","2","3","4","5","6"), 
                            to=c("All","-ST","-cremoris(5614)","-lactis(6086)","-lactis(10675)","-SICO"))
  data$clean_ids_names<-mapvalues(data$clean_ids, 
                                  from=c("CHCC5614","CHCC6086","CHCC10675","SICO"), 
                                  to=c("cremoris(5614)","lactis(6086)","lactis(10675)","SICO"))
  #colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')
  colors_t0s<-c('#313695','#4575b4','#abd9e9','#1b7837')
  size_all<-20
  data$condition<-factor(data$condition,levels = c( "All", "-lactis(10675)","-cremoris(5614)","-ST","-lactis(6086)","-SICO"))
  at<-ggplot(data,aes(x=clean_ids_names,y=value,fill=clean_ids_names))+geom_boxplot(position=position_dodge(1))+ theme_bw()+
    labs(x="", y= "mRNA TPM",fill="strains",title="Citrate-sodium symporter")+geom_jitter(color="black", size=0.6, alpha=0.9,aes(fill=clean_ids_names))+#geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=2)+
    theme(text = element_text(size=size_all),
          plot.title = element_text(size=size_all),
          axis.title.x = element_text(size=size_all),
          axis.title.y = element_text( size=size_all), 
          legend.title = element_text(size=size_all),
          axis.text.x = element_text(angle = 45,hjust=1),
          # axis.text.x=element_blank(),
          # axis.ticks.x=element_blank(),
          legend.text = element_text(size=size_all))+facet_wrap(~condition)+scale_fill_manual(values = colors_t0s)
  at
