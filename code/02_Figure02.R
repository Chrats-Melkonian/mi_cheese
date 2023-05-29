################ Author: Chrats Melkonian
#first set directory to github repository location
################
setwd("~/Desktop/Git/")
data <- read.csv("./mi_cheese/data/T1/Biochemical/Cheese Tier 1 biochemical rawdata_v3_14022020.csv")
data_cnames<-colnames(data)
data_cnamesnew<-data[1,]
data<-data[-1,]
colnames(data)<-data_cnamesnew

data_sub<-as.data.frame(data[,c(16:55,59:ncol(data))])
data_meta<-data[,c(1:4)]
temp_col<-colnames(data_sub)
data_sub[data_sub=="<LOD"]<-0
data_sub<-t(apply(as.matrix(data_sub), 1, as.numeric))
colnames(data_sub)<-temp_col
data_sub[,31:40]<-data_sub[,31:40]*1000

index<-!apply(data_sub, 2, function(x) all(x==0)) # remove all 0
data_sub<-data_sub[,index]

qant<-apply(data_sub, 2, function(x) quantile(x,na.rm=T))
data_sub<-data_sub[,-which(qant[4,]==0)] #remove compound with quantile 75% = 0
#add mean values to NA for UMAP
data_sub_umap<-data_sub
for (i in 1:nrow(data_sub_umap)){
  x<-data_sub_umap[,i]
  ifelse(length(table(is.na(x)))==2, data_sub_umap[is.na(x),i]<-mean(x, na.rm = T) ,x)
}

library(uwot)
embedding<-umap(data_sub_umap,n_neighbors=23)
df_out <- as.data.frame(embedding)### nn 23
df_out<-cbind(df_out,data_meta)

colors_same<-c('black','#005824','#238b45','#41ae76','#66c2a4')
library(ggplot2)
library(ggrepel)
library(scales)
#
colnames(df_out)[4]<-"Culture"
colnames(df_out)[5]<-"Ripening_time"
colnames(df_out)[6]<-"Ripening_months"
theme_set(theme_bw())

df_out$Ripening_months<-factor(df_out$Ripening_months,levels = c(0.5,3,6,9,12))
size_all<-20

umap1<-ggplot(df_out,aes(x=V1,y=V2,color=Ripening_months))+geom_point(aes( shape=Culture),size=3,stroke=1)+scale_color_manual(name="Time", values= colors_same)+
  stat_ellipse(aes(x=V1, y=V2,lty=Ripening_months),type = "norm",size=1)+scale_linetype_manual(name="Time",values=c(5,4,3,2,1))+
  guides(colour = guide_legend(override.aes = list(size=5)),shape = guide_legend(override.aes = list(size=5)))+ labs(shape="Condition",x="UMAP1",y="UMAP2")+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        legend.text = element_text(size=size_all))
umap1

df_out_zoom<-df_out[!df_out$Ripening_months%in%"0.5",]
colors_same<-c('#005824','#238b45','#41ae76','#66c2a4')
umap2<-ggplot(df_out_zoom,aes(x=V1,y=V2,color=Ripening_months))+geom_point(aes( shape=Culture),size=3,stroke=1)+scale_color_manual(name="Time", values= colors_same)+
  stat_ellipse(aes(x=V1, y=V2,lty=Ripening_months),type = "norm",size=1)+scale_linetype_manual(name="Time",values=c(4,3,2,1))+
  guides(colour = guide_legend(override.aes = list(size=5)),shape = guide_legend(override.aes = list(size=5)))+ labs(shape="Condition",x="UMAP1",y="UMAP2")+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        legend.text = element_text(size=size_all))

umap2

##########################

embedding2<-umap(t(data_sub_umap))
df_out2 <- as.data.frame(embedding2)### nn 23
df_out2$names<-colnames(data_sub_umap)
Class<-c(
  rep("Acids",30),
  rep("Carbohydrates",31),
  rep("Peptides",248))
df_out2$Class<-as.factor(Class)

umap3<-ggplot(df_out2,aes(x=V1,y=V2,color=Class))+geom_point(aes( shape=Class),size=3,stroke=1)+scale_color_manual(name="Time", values= c('#7fc97f','#beaed4','#fdc086'))+
 # stat_ellipse(aes(x=V1, y=V2,lty=Class, group=Class),type = "norm",size=1.5)+scale_linetype_manual(name="Time",values=c(6,7,8))+
  guides(colour = guide_legend(override.aes = list(size=5)),shape = guide_legend(override.aes = list(size=5)))+ labs(shape="",x="UMAP1",y="UMAP2")+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        legend.text = element_text(size=size_all))

umap3
library(apcluster)
apres <- apcluster(negDistMat(r=2), embedding2, details=TRUE,q=0.01)
plot(apres,embedding2)
apres@clusters[[1]]
temp_scluster<-data_sub_umap[,apres@clusters[[1]]]

temp_scluster<-cbind(df_out$Culture,temp_scluster)

################


pca_res <- prcomp(data_sub_umap, scale. = TRUE,center = T)
library(ggfortify)
autoplot(pca_res, data=df_out,colour="Ripening_months",frame = TRUE,frame.type = 'norm',shape="Culture")

data_sub_umap_scale<-apply(data_sub_umap, 2, function(x) rescale(x, to = c(0, 1)))
pca_res2 <- prcomp(t(data_sub_umap), scale. = TRUE,center = T)
Class_PCA<-c(
  rep("Acids",30),
  rep("Carbohydrates",31),
  rep("Peptides",248))
df<-data.frame(1:309)
df$Class_PCA<-Class_PCA
autoplot(pca_res2,df, colour="Class_PCA",label = TRUE, label.size = 3)

################
################ relative change between times figure 2 f
################

data_sub
datac<-cbind(data_meta,data_sub)
colnames(datac)[2]<-"Culture"
colnames(datac)[3]<-"Ripening_time"
colnames(datac)[4]<-"Ripening_months"
library(dplyr)

datac_avg_time<-datac %>%                                        # Specify data frame
  dplyr::group_by(across(all_of(c("Ripening_months")))) %>%                         # Specify group indicator             # Specify function
  dplyr::summarise(across(everything(), list(mean)))
datac_avg_time<-datac_avg_time[,-c(2,3,4)]
datac_avg_time<-datac_avg_time[c(1,3,4,5,2),]
datac_avg_time<-datac_avg_time[,!apply(datac_avg_time, 2, function(x) all(is.na(x)))]

library(scales)
library(reshape2)
datac_avg_time_nm<-apply(as.matrix(datac_avg_time), 2, as.numeric)
# datac_avg_time_nm_rescale<-apply(datac_avg_time_nm, 2, function(x) rescale(x, to = c(1, 1)))
datac_avg_time_nm[,-1]<-datac_avg_time_nm[,-1]+1
rm_05_3<-apply(datac_avg_time_nm[,-1], 2, function(x) (x[2]-x[1])/x[1])
rm_3_6<-apply(datac_avg_time_nm[,-1], 2, function(x) (x[3]-x[2])/x[2])
rm_6_9<-apply(datac_avg_time_nm[,-1], 2, function(x) (x[4]-x[3])/x[3])
rm_9_12<-apply(datac_avg_time_nm[,-1], 2, function(x) (x[5]-x[4])/x[4])
data_rc_time<-rbind(rm_05_3,rm_3_6,rm_6_9,rm_9_12)
rownames(data_rc_time)<-c("05_3","3_6","6_9","9_12")
data_rc_time_df<-melt(data_rc_time)

Class<-c(
rep(rep("Acids",30),4),
rep(rep("Sugars",3),4),
rep(rep("Flavor-related compounds",28),4),
rep(rep("Peptides",248),4))
data_rc_time_df$Class<-as.factor(Class)
data_rc_time_df_noS<-data_rc_time_df[!data_rc_time_df$Class=="Sugars",]

data_rc_time_df_noS$value[is.infinite(data_rc_time_df_noS$value)]<-NA

################ 
color_met<-c("#d6604d", "#878787","#4393c3" )
data_rc_time_df_noS$value_log<-log10(data_rc_time_df_noS$value)
data_rc_time_df_noS$value_log[is.infinite(data_rc_time_df_noS$value_log)]<-NA


rc_m<-ggplot(data_rc_time_df_noS, aes(Var1, value_log,  fill = Class, shape = Class)) +geom_boxplot(outlier.shape=NA,show.legend = FALSE)+
  geom_point(position=position_jitterdodge(), alpha = 0.8, size = 2) +
  scale_fill_manual(values=color_met) +
  scale_shape_manual(values= c(23, 24, 25)) +
  theme_bw()+labs(shape="Class",x="Time intervals",y="Log10(relative change)")+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        legend.text = element_text(size=size_all))+ guides(fill = guide_legend(override.aes = list(size=5)))
rc_m

subset(data_rc_time_df_noS,Class=="Acids")
#05_3 	Acetic acid_1, 3-6 Tyramine_1, g-Amino Butyric acid_1 	g-Amino Butyric acid (GABA)_1 Putrescine_1 Cadaverine_1 Butyric acid_1
subset(data_rc_time_df_noS,Class=="Flavor-related compounds")
#	2-Heptanone_1, 2,6-dimethylpyrazine_1 Hexanal_1 , 	2-Nonanone_1 , 	2-Pentanone_1, 2-Heptanone_1 >1 
subset(data_rc_time_df_noS,Class=="Peptides")


################ Supplementary

data_meta<-data[,c(1:4)]

data_sub
################ 
Conf_level<-0.95
size_all<-15
h<-1
plot<-list()
for (j in colnames(data_sub)){
  data_temp<-data_sub[,colnames(data_sub)==j]
  data_temp_meta<-cbind(data_temp,data_meta)
  data_temp_meta$batch<-  ifelse(grepl("-3-",data_temp_meta$Sample)==F,"A","B")
  data_temp_meta$`Ripening months`<-factor(data_temp_meta$`Ripening months`, levels = c(unique(data_temp_meta$`Ripening months`)))
  data_temp_meta$`Culture`<-gsub("A3020","All",data_temp_meta$`Culture`)
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
  geom_point(data=data_temp_meta, aes(x=`Ripening months`, y=data_temp,shape=batch,color=`Culture`),size=4) +
    
    #### plot average response over time ####
  geom_line(data=df_summary_all, aes(x=Time, y=mean,group=Culture,color=Culture,linetype=Culture), size=2, alpha=0.8)+
    scale_linetype_manual(name="Condition",values=c(4,3,2,1))+
    
    #### plot error (95%CI) of the response over time ####
  geom_ribbon(data=df_summary_all, aes(ymin=CI_lower, ymax=CI_upper,group=Culture,fill=Culture), alpha=0.2)+
    labs(x="",y=paste(j,"(mg/Kg)",sep = " "),color="Condition",shape="Batch")+
    theme(text = element_text(size=size_all),
          plot.title = element_text(size=size_all),
          axis.title.x = element_text(size=size_all),
          axis.title.y = element_text( size=size_all), 
          legend.title = element_text(size=15),
          legend.text = element_text(size=15),
          legend.position="top")+scale_color_manual(values = c('#fc8d62','#8da0cb','#e78ac3','#66c2a5'))+scale_fill_manual(values = c('#fc8d62','#8da0cb','#e78ac3','#66c2a5'))+
    guides(fill=F,linetype=guide_legend(keywidth = 4, keyheight = 1))+scale_shape_manual(values = c(19,8))
  h<-h+1
}
names(plot)<-colnames(data_sub)
#05_3 	Acetic acid_1, 3-6 Tyramine_1, g-Amino Butyric acid_1 	g-Amino Butyric acid (GABA)_1 Putrescine_1 Cadaverine_1 Butyric acid_1

library(cowplot)
legend1 <- get_legend(plot$`Acetic acid`)

prow_acids <- plot_grid(plot$`Acetic acid` + theme(legend.position="none"),
                  plot$Tyramine + theme(legend.position="none"), 
                  plot$`g-Amino Butyric acid` + theme(legend.position="none"), 
                  plot$Putrescine+theme(legend.position="none"), 
                  plot$Cadaverine+theme(legend.position="none"), 
                  labels=c("a", "b","c","d","e",""), ncol = 2, nrow = 3, label_fontface = "bold",label_size=25,align = "hv")
prow_acids2<- plot_grid(prow_acids,legend1,ncol = 1, nrow = 2,rel_heights=c(2,0.2),rel_widths =c(1,0.5))
prow_acids2
# ggsave(".../SUPP/ACIDS_RC_supp.pdf",prow_acids2, units="in", width=10, height=15, dpi=600)


#Carbohydrates
#	2-Heptanone_1, 2,6-dimethylpyrazine_1 Hexanal_1 , 	2-Nonanone_1 , 	2-Pentanone_1, 2-Heptanone_1 >1 


prow_carb <- plot_grid(plot$`2-Heptanone` + theme(legend.position="none"),
                        plot$`2,6-dimethylpyrazine` + theme(legend.position="none"), 
                        plot$Hexanal + theme(legend.position="none"), 
                        plot$`2-Nonanone`+theme(legend.position="none"), 
                        plot$`2-Pentanone`+theme(legend.position="none"), 
                        labels=c("a", "b","c","d","e",""), ncol = 2, nrow = 3, label_fontface = "bold",label_size=25,align = "hv")
prow_carb2<- plot_grid(prow_carb,legend1,ncol = 1, nrow = 2,rel_heights=c(2,0.2),rel_widths =c(1,0.5))
prow_carb2
#ggsave("Git/cheese-ft/PLOTS/SUPP/CARB_RC_supp.pdf",prow_carb2, units="in", width=10, height=15, dpi=600)


################ feature selection figure 2g

#############################################
index<-data_meta$`Ripening months`%in%c("0.5","3")
data_sub2<-data_sub[!index,]
data_meta_sub2<-data_meta[!index,]


y<-as.factor(data_meta_sub2$`Culture`)

library(Boruta)

data_sub_boruta<-data_sub2
for (i in 1:nrow(data_sub_boruta)){
  x<-data_sub_boruta[,i]
  ifelse(length(table(is.na(x)))==2, data_sub_boruta[is.na(x),i]<-mean(x, na.rm = T) ,x)
}

colnames(data_sub_boruta)<-paste(1:309,colnames(data_sub_boruta),sep = "_")
boruta.train <- Boruta(x=data_sub_boruta, y=y,doTrace = 2,maxRuns=10000)

plot(boruta.train, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta.train$ImpHistory),function(i)
  boruta.train$ImpHistory[is.finite(boruta.train$ImpHistory[,i]),i])
names(lz) <- colnames(boruta.train$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta.train$ImpHistory), cex.axis = 0.7)

final.boruta <- TentativeRoughFix(boruta.train)
boruta.df <- attStats(final.boruta)
boruta.df.ord<-boruta.df[order(boruta.df$meanImp,decreasing = T),]
boruta.df.ord.c<-boruta.df.ord[boruta.df.ord$decision=="Confirmed",]
boruta.df.ord.c.h<-boruta.df.ord.c[boruta.df.ord.c$meanImp>mean(boruta.df.ord.c$meanImp),]
rownames(boruta.df.ord.c.h)<-gsub("^X","",rownames(boruta.df.ord.c.h))



colnames(data_sub2)<-paste(1:309,colnames(data_sub2),sep = "_")
select<-match(unlist(lapply(strsplit(rownames(boruta.df.ord.c.h),"_"), function(x) x[1])),unlist(lapply(strsplit(colnames(data_sub2),"_"), function(x) x[1])))
data_sub2_select<-data_sub2[,select]
Class_fs<-Class_PCA[select]


Condition<-data_meta_sub2$`Culture`
Ripening_months<-data_meta_sub2$`Ripening months`
data_sub2_select_rs<-apply(data_sub2_select, 2, rescale)
data_sub3_select<-cbind(Condition,Ripening_months,as.data.frame(data_sub2_select_rs))
data_sub3_select_df<-melt(data_sub3_select)


library(ggridges)

temp<-factor(Class_fs)
data_bar<-as.data.frame(temp)
data_bar$temp<-factor(data_bar$temp)
data_bar$size<-rep(1,length(temp))
data_bar$map<-rep("A",length(temp))
ggplot(data_bar,aes(y=size,x=map,fill=temp,order=temp))+  geom_bar( stat="identity")
a<-plyr::mapvalues(Class_fs, 
          from=unique(Class_fs), 
          to=color_met)
data_sub3_select_df$variable<-factor(data_sub3_select_df$variable,levels = rev(colnames(data_sub2_select)))

data_sub3_select_df$Condition<-gsub("A3020","All",data_sub3_select_df$Condition)
data_sub3_select_df$Condition<-gsub("SICO","LM",data_sub3_select_df$Condition)

sel_ridges<-c("31_Galactose","33_Lactose","30_Lactic acid","205_EEEKNRLNF_α-S2_170-178","170_VNELS(+79.97)KD_α-S1_52-58","106_ELSKDIGS(+79.97)ES(+79.97)TE_α-S1_54-65")
sel_ridges1<-c("205_EEEKNRLNF_α-S2_170-178","170_VNELS(+79.97)KD_α-S1_52-58","106_ELSKDIGS(+79.97)ES(+79.97)TE_α-S1_54-65","30_Lactic acid","33_Lactose","31_Galactose")

data_sub3_select_df_sel<-data_sub3_select_df[data_sub3_select_df$variable%in%sel_ridges,]
data_sub3_select_df_sel$variable<-factor(data_sub3_select_df_sel$variable,levels = sel_ridges1)
df_bc<-ggplot(data_sub3_select_df_sel, aes(x=value, y=variable))+ 
  geom_density_ridges(aes(fill=Condition,color=Condition),scale = 1,alpha = .3,jittered_points = TRUE, point_size =2, size = 0.1,
                      rel_min_height = .05, point_alpha = 1)+
  scale_color_manual(values = c('#fc8d62','#8da0cb','#e78ac3','#66c2a5'))+scale_fill_manual(values = c('#fc8d62','#8da0cb','#e78ac3','#66c2a5'))+
  theme_ridges(center = TRUE)+coord_cartesian(clip = "off") +scale_y_discrete(expand = c(0, 0)) +labs(x="Scaled concentrations",y="Discriminative features")+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        legend.text = element_text(size=size_all),axis.text.y=element_blank())

#selected peptides ELSKDIGSESTE , VNELSKD, EEEKNRLNF

#sub select original peptide data
data_pep<-cbind(data[,c(2:4)],apply(data[,87:334], 2, function(x) rescale(as.integer(x))))
data_pep_melt<-melt(data_pep, id=c("Culture ", "Ripening time","Ripening months"))
temp_size<-unlist(lapply(strsplit(as.character(data_pep_melt$variable),'_'), function(x) x[length(x)]))
temp_seq<-unlist(lapply(strsplit(as.character(data_pep_melt$variable),'_'), function(x) x[1]))
data_pep_melt$seq<-gsub("[(].*[)]","",temp_seq)
data_pep_melt$pep_size<-nchar(data_pep_melt$seq)
data_pep_melt$value<-as.numeric(data_pep_melt$value)
ggplot(data_pep_melt,aes(value,pep_size))+geom_point()+facet_wrap(~`Ripening time`)
#### one-way ANOVA


library(rstatix)
data_sub3_select_df_sel_list<-split(data_sub3_select_df_sel,factor(data_sub3_select_df_sel$variable))
lapply(data_sub3_select_df_sel_list, function(x) x %>% anova_test(value ~ Condition*Ripening_months))
mean(c(33.556, 60.853, 63.572 ,66.294 ,137.525 ,886.441))
sd(c(33.556, 60.853, 63.572 ,66.294 ,137.525 ,886.441))
mean(c(0.835,0.841,0.847,0.737,0.987,0.920))
sd(c(0.835,0.841,0.847,0.737,0.987,0.920))
#F(3, 36) = 208+-334 , p<0.001, eta2[g] = 0.86+-0.08

####
ggplot(data_sub3_select_df_sel, aes(x=value, y=variable))+ 
  geom_density_ridges(aes(fill=Condition,color=Condition),scale = 1,alpha = .3,jittered_points = TRUE, point_size =2, size = 0.1,
                      rel_min_height = .05, point_alpha = 1)+
  scale_color_manual(values = c('#fc8d62','#8da0cb','#e78ac3','#66c2a5'))+scale_fill_manual(values = c('#fc8d62','#8da0cb','#e78ac3','#66c2a5'))+
  theme_ridges(center = TRUE)+coord_cartesian(clip = "off") +scale_y_discrete(expand = c(0, 0)) +labs(x="Scaled concentrations",y="Discriminative features")+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        legend.text = element_text(size=size_all))

### needs to be modified based on the order 
df_bc_final<-df_bc + annotate("text", x = 1.5, y = 6, label = "bold(Galactose)",color="black",parse=T,size=3)+
  annotate("text", x = 1.5, y = 6-1, label = "bold(Lactose)",color="black",parse=T,size=3)+
  annotate("text", x = 1.5, y = 6-2, label = "bold(Lactic_acid)",color="#d6604d",parse=T,fontface='bold',size=3)+
  annotate(geom = "point", x = 1.5, y = 6-3,shape=25,color="#4393c3" ,fill="#4393c3",size=3 )+
  annotate(geom = "point", x = 1.5, y = 6-4,shape=25,color="#4393c3" ,fill="#4393c3",size=3)+
  annotate(geom = "point", x = 1.5, y = 6-5,shape=25,color="#4393c3" ,fill="#4393c3",size=3)

df_bc_final

color_met
library(cowplot)
legend1 <- get_legend(umap1)
legend2 <- get_legend(rc_m)
legend3 <- get_legend(df_bc_final)
prow <- plot_grid(umap1 + theme(legend.position="none"),
                  umap2 + theme(legend.position="none"), 
                  rc_m_final + theme(legend.position="none"), 
                  df_bc_final+theme(legend.position="none"), 
                  labels=c("a", "b","c","d"), ncol = 2, nrow = 2, label_fontface = "bold",label_size=25,align = "hv")
#ggsave("figure1.pdf", units="in", width=14, height=9, dpi=600)
legend_row2<-plot_grid(legend1,legend2,nrow =2,ncol = 1,align = "hv")

final<-plot_grid(prow,legend,ncol=2,rel_widths = c(2,0.5),align = "hv")
# ggsave("Git/cheese-ft/PLOTS/figure2.1.pdf",final, units="in", width=14, height=11, dpi=600)
# ggsave("Git/cheese-ft/PLOTS/figure2.1.png",final, units="in", width=14, height=11, dpi=600)
plot(prow)
#alternative supplementary

unique(data_sub3_select_df$variable)
legend1 <- get_legend(plot$`Acetic acid`)
prow_pep <- plot_grid(plot$`EEEKNRLNF_α-S2_170-178` + theme(legend.position="none"),
                       plot$`S(+79.97)S(+79.97)S(+79.97)EESITRINK_β_32-43` + theme(legend.position="none"), 
                       plot$`ELSKDIGS(+79.97)ES(+79.97)TE_α-S1_54-65` + theme(legend.position="none"), 
                       plot$`KHPIKHQGLPQ_α-S1_18-28`+theme(legend.position="none"), 
                       plot$`TPVVVPP_β_95-101`+theme(legend.position="none"), 
                       plot$`IVPNS(+79.97)AEERL_α-S1_126-135`+theme(legend.position="none"),
                       labels=c("a", "b","c","d","e","f"), ncol = 2, nrow = 3, label_fontface = "bold",label_size=25,align = "hv")
prow_pep<- plot_grid(prow_pep,legend1,ncol = 1, nrow = 2,rel_heights=c(2,0.2),rel_widths =c(1,0.5))
prow_pep
# ggsave("Git/cheese-ft/PLOTS/SUPP/PEP_cond_supp.pdf",prow_pep, units="in", width=10, height=15, dpi=600)


library(patchwork)
patchwork<- (m + theme(legend.position="none")|
               p_12+theme(legend.position="none")|
               t + theme(legend.position="none")|
               n + theme(legend.position="none"))/(
               umap1 + theme(legend.position="none")|
               rc_m + theme(legend.position="none")|
               df_bc_final + theme(legend.position="none")|
                 plot[[32]] + theme(legend.position="none"))
patchworkf<-patchwork+plot_annotation(tag_levels = 'a')
patchworkf
# ggsave("Git/cheese-ft/PLOTS/FINAL_4main/figure1_05082022.pdf",patchworkf, units="in", width=24, height=12, dpi=600)
# ggsave("Git/cheese-ft/PLOTS/FINAL_4main/figure1_05082022.png", units="in", width=24, height=12, dpi=600)

legend_row1<-plot_grid(legend,nrow = 1,ncol = 1,align = "hv")
legend_row1
# ggsave("Git/cheese-ft/PLOTS/FINAL_4main/figure1_legenedr1.pdf", units="in", width=24, height=12, dpi=600)
legend_row2
# ggsave("Git/cheese-ft/PLOTS/FINAL_4main/figure1_legenedr2.pdf", units="in", width=24, height=12, dpi=600)
