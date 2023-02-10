################ Author: Chrats Melkonian
#first set directory to github repository location
################
# Remark: through the analysis or data the follow abbreviation can be used in exchange 10675: LLm1, 6086: LLm2, 5614: LC, SICO: LB.
################
setwd("~/Desktop/Git/")
data <- read.csv("./mi_cheese/data/T1/Biochemical/Cheese Tier 1 biochemical rawdata_v3_14022020.csv")
data_cnames<-colnames(data)
data_cnamesnew<-data[1,]
data<-data[-1,]
colnames(data)<-data_cnamesnew
data_meta<-data[,c(1:4)]

library(viridis)
color<-viridis(4)
library(cowplot)
library(ggpubr)
library(reshape2)

data_colonies<-data[,56:58] ### CFU/G
data_colonies<-t(apply(as.matrix(data_colonies), 1, function(x) as.numeric(x)+1))
data_colonies_log10<-log10(data_colonies)
colnames(data_colonies_log10)<-colnames(data[,56:58])
data_colonies_log10m<-melt(data_colonies_log10)

df_out_cfu<-cbind(data_colonies_log10m,rbind(data_meta,data_meta,data_meta))
colnames(df_out_cfu)[5]<-"Culture"
colnames(df_out_cfu)[6]<-"Ripening_time"
colnames(df_out_cfu)[7]<-"Ripening_months"
size_all<-20
df_out_cfu$Ripening_months<-factor(df_out_cfu$Ripening_months,levels = c(0.5,3,6,9,12))
ggplot(df_out_cfu,aes(x=Ripening_months,y=value,color=Var2))+geom_point()+geom_jitter(width = 0.05,size=2)+facet_grid(vars(Culture))+theme_bw()+
  labs(y="log10(CFU/g)", x="Time (months)",color="Taxonomy group")+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        legend.text = element_text(size=size_all))+scale_color_manual(values = c('#fc8d62','#66c2a5','#8da0cb'))

################
#2 weeks
################
newdata_m_05 <- subset(df_out_cfu, Ripening_months =="0.5"& Var2=="Mesophilic cocci" ) 

dat <- as.data.frame(newdata_m_05)
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("All", "All_HP"), c("All", "All-ST"),c("All", "All-LB")) # comparisons for post-hoc tests
milk_data_mean<-mean(dat$value)

dat$batch<-  ifelse(grepl("-3-",dat$Sample)==F,"A","B")
dat$Culture<-factor(dat$Culture,levels = c( "All",     "All_HP", "All-ST", "All-LB"))
p <- ggboxplot(dat,
               x = "Culture", y = "value",
               color = "Culture",
               legend = "none",
               palette = c('#fc8d62','#8da0cb','#66c2a5','#e78ac3'),
               add = "jitter", size = 1) + theme(
                 axis.text.x = element_text(angle = 45, hjust = 1))+ geom_hline(yintercept=milk_data_mean)+
  labs(y="log10(CFU/g)", x="",title="Mesophilic cocci - month 0.5")+theme(text = element_text(size=size_all),
                                                                          plot.title = element_text(size=size_all),
                                                                          axis.title.x = element_text(size=size_all),
                                                                          axis.title.y = element_text( size=size_all), 
                                                                          legend.title = element_text(size=size_all),
                                                                          legend.text = element_text(size=size_all))

p_05<-print(
  p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
                         method = method1, label.y = min(dat$value, na.rm = TRUE),size=4
  )
  + stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format",size=3) # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
) 


################
# newdata_m_9 <- subset(df_out_cfu, Ripening_months =="9"& Var2=="Mesophilic cocci" ) 
# 
# dat <- as.data.frame(newdata_m_9)
# method1 <- "anova" # one of "anova" or "kruskal.test"
# method2 <- "t.test" # one of "wilcox.test" or "t.test"
# my_comparisons <- list(c("A3020", "A3020_HP"), c("A3020", "A3020-ST"), c("A3020", "A3020-SICO")) # comparisons for post-hoc tests
# milk_data_mean<-mean(dat$value)

# p <- ggboxplot(dat,
#                x = "groups", y = "value",
#                color = "groups",
#                legend = "none",
#                palette = "npg",
#                add = "jitter") + theme(
#                  axis.text.x = element_text(angle = 45, hjust = 1))+ geom_hline(yintercept=milk_data_mean)
# print(
#   p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
#                        method = method1, label.y = max(dat$value, na.rm = TRUE)
#   )
#   + stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
# ) +labs(y="log2(CFU/g)", x="Culture condition",title="Mesophilic cocci - month 9")


################
newdata_m_12 <- subset(df_out_cfu, Ripening_months =="12"& Var2=="Mesophilic cocci" ) 

dat <- as.data.frame(newdata_m_12)
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <-  list(c("All", "All_HP"), c("All", "All-ST"),c("All", "All-LB")) # comparisons for post-hoc tests
milk_data_mean<-mean(dat$value)

p <- ggboxplot(dat,
               x = "Culture", y = "value",
               color = "Culture",
               legend = "none",
               palette =c('#fc8d62','#8da0cb','#66c2a5','#e78ac3'),
               add = "jitter",size = 1) + theme(
                 axis.text.x = element_text(angle = 45, hjust = 1))+ geom_hline(yintercept=milk_data_mean)+
  labs(y="", x="",title="Mesophilic cocci - month 12")+theme(text = element_text(size=size_all),
                                                             plot.title = element_text(size=size_all),
                                                             axis.title.x = element_text(size=size_all),
                                                             axis.title.y = element_text( size=size_all), 
                                                             legend.title = element_text(size=size_all),
                                                             legend.text = element_text(size=size_all))
p_12<- print(
  p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
                         method = method1,label.y = min(dat$value, na.rm = TRUE),size=5)
  + stat_compare_means(aes(label = ..p.signif..),
                       method = "t.test", ref.group = "All",size=7)
    #stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.signif",size=7) # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
) 


##################################3
newdata_nslab <- subset(df_out_cfu, Var2=="NSLAB" ) 
Conf_level <-  0.95
newdata_nslab$batch<-  ifelse(grepl("-3-",newdata_nslab$Sample)==F,"A","B")
df_summary_all<-c()
for(i in unique(newdata_nslab$Culture)){
  temp<-newdata_nslab[newdata_nslab$Culture==i,]
  df_summary <- data.frame( n=tapply(temp$value, temp$Ripening_months, length),Culture=rep(i,length(unique(temp$Ripening_months))), mean=tapply(temp$value, temp$Ripening_months, function(x) mean(x,na.rm = T)))
  df_summary$sd <- tapply(temp$value, temp$Ripening_months, function(x) sd(x,na.rm = T))
  df_summary$sem <- df_summary$sd/sqrt(df_summary$n-1)
  df_summary$CI_lower <- df_summary$mean + qt((1-Conf_level)/2, df=df_summary$n-1)*df_summary$sem
  df_summary$CI_upper <- df_summary$mean - qt((1-Conf_level)/2, df=df_summary$n-1)*df_summary$sem
  df_summary$Time=rownames(df_summary)
  df_summary_all<-rbind(df_summary_all,df_summary)
  
}


#### Command to prepare the plot ####
n<- ggplot(df_summary_all, aes(x=Time, y=mean)) + theme_bw()+ylim(-0.7,10)+
  
  #### plot individual measurements ####
geom_point(data=newdata_nslab, aes(x=Ripening_months, y=value,color=Culture,shape=batch),size=4) +
  
  #### plot average response over time ####
geom_line(data=df_summary_all, aes(x=Time, y=mean,group=Culture,color=Culture,linetype=Culture), size=2, alpha=0.8)+
  scale_linetype_manual(name="Condition",values=c(4,3,2,1))+
  
  #### plot error (95%CI) of the response over time ####
geom_ribbon(data=df_summary_all, aes(ymin=CI_lower, ymax=CI_upper,group=Culture,fill=Culture), alpha=0.2)+
  labs(y="", x="Time (months)",color="Conditions",shape="Batch",title = "NSLAB")+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        legend.text = element_text(size=size_all))+scale_color_manual(values = c('#fc8d62','#8da0cb','#e78ac3','#66c2a5'))+scale_fill_manual(values = c('#fc8d62','#8da0cb','#e78ac3','#66c2a5'))+
  guides(fill=F,linetype=guide_legend(keywidth = 4, keyheight = 1))+scale_shape_manual(values = c(19,8))
n
##############################

newdata_m <- subset(df_out_cfu, Var2=="Thermophilic cocci" ) 
Conf_level <-  0.95
newdata_m$batch<-  ifelse(grepl("-3-",newdata_m$Sample)==F,"A","B")
df_summary_all<-c()
for(i in unique(newdata_m$Culture)){
  temp<-newdata_m[newdata_m$Culture==i,]
  df_summary <- data.frame( n=tapply(temp$value, temp$Ripening_months, length),Culture=rep(i,length(unique(temp$Ripening_months))), mean=tapply(temp$value, temp$Ripening_months, mean))
  df_summary$sd <- tapply(temp$value, temp$Ripening_months, sd)
  df_summary$sem <- df_summary$sd/sqrt(df_summary$n-1)
  df_summary$CI_lower <- df_summary$mean + qt((1-Conf_level)/2, df=df_summary$n-1)*df_summary$sem
  df_summary$CI_upper <- df_summary$mean - qt((1-Conf_level)/2, df=df_summary$n-1)*df_summary$sem
  df_summary$Time=rownames(df_summary)
  df_summary_all<-rbind(df_summary_all,df_summary)
  
}


#### Command to prepare the plot ####
t<- ggplot(df_summary_all, aes(x=Time, y=mean)) + theme_bw()+ylim(-0.7,10)+
  
  #### plot individual measurements ####
geom_point(data=newdata_m, aes(x=Ripening_months, y=value,color=Culture,shape=batch),size=4) +
  
  #### plot average response over time ####
geom_line(data=df_summary_all, aes(x=Time, y=mean,group=Culture,color=Culture,linetype=Culture), size=2, alpha=0.8)+
  scale_linetype_manual(name="Condition",values=c(4,3,2,1))+
  
  #### plot error (95%CI) of the response over time ####
geom_ribbon(data=df_summary_all, aes(ymin=CI_lower, ymax=CI_upper,group=Culture,fill=Culture), alpha=0.2)+
  labs(y="", x="Time (months)",color="Condition",shape="Batch",title = "Thermophilic cocci")+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        legend.text = element_text(size=size_all))+scale_color_manual(values = c('#fc8d62','#8da0cb','#e78ac3','#66c2a5'))+scale_fill_manual(values = c('#fc8d62','#8da0cb','#e78ac3','#66c2a5'))+
  guides(fill=F,linetype=guide_legend(keywidth = 4, keyheight = 1))+scale_shape_manual(values = c(19,8))
t
############################
newdata_nslab <- subset(df_out_cfu, Var2=="Mesophilic cocci" ) 
Conf_level <-  0.95
newdata_nslab$batch<-  ifelse(grepl("-3-",newdata_nslab$Sample)==F,"A","B")
df_summary_all<-c()
for(i in unique(newdata_nslab$Culture)){
  temp<-newdata_nslab[newdata_nslab$Culture==i,]
  df_summary <- data.frame( n=tapply(temp$value, temp$Ripening_months, length),Culture=rep(i,length(unique(temp$Ripening_months))), mean=tapply(temp$value, temp$Ripening_months, function(x) mean(x,na.rm = T)))
  df_summary$sd <- tapply(temp$value, temp$Ripening_months, function(x) sd(x,na.rm = T))
  df_summary$sem <- df_summary$sd/sqrt(df_summary$n-1)
  df_summary$CI_lower <- df_summary$mean + qt((1-Conf_level)/2, df=df_summary$n-1)*df_summary$sem
  df_summary$CI_upper <- df_summary$mean - qt((1-Conf_level)/2, df=df_summary$n-1)*df_summary$sem
  df_summary$Time=rownames(df_summary)
  df_summary_all<-rbind(df_summary_all,df_summary)
  
}


#### Command to prepare the plot ####
m<- ggplot(df_summary_all, aes(x=Time, y=mean)) + theme_bw()+ylim(-0.7,10)+
  
  #### plot individual measurements ####
geom_point(data=newdata_nslab, aes(x=Ripening_months, y=value,color=Culture,shape=batch),size=4) +
  
  #### plot average response over time ####
geom_line(data=df_summary_all, aes(x=Time, y=mean,group=Culture,color=Culture,linetype=Culture), size=2, alpha=0.8)+
  scale_linetype_manual(name="Condition",values=c(4,3,2,1))+
  
  #### plot error (95%CI) of the response over time ####
geom_ribbon(data=df_summary_all, aes(ymin=CI_lower, ymax=CI_upper,group=Culture,fill=Culture), alpha=0.2)+
  labs(y="log10(CFU/g)", x="Time (months)",color="Condition",shape="Batch",title = "Mesophilic cocci")+
  theme(text = element_text(size=size_all),
        plot.title = element_text(size=size_all),
        axis.title.x = element_text(size=size_all),
        axis.title.y = element_text( size=size_all), 
        legend.title = element_text(size=size_all),
        legend.text = element_text(size=size_all))+scale_color_manual(values = c('#fc8d62','#8da0cb','#e78ac3','#66c2a5'))+scale_fill_manual(values = c('#fc8d62','#8da0cb','#e78ac3','#66c2a5'))+
  guides(fill=F,linetype=guide_legend(keywidth = 4, keyheight = 1))+scale_shape_manual(values = c(19,8))
m
########################
legend <- get_legend(m)
prow <- plot_grid(m + theme(legend.position="none"),
                  t + theme(legend.position="none"), 
                  n + theme(legend.position="none"), 
                  p_05+theme(legend.position="none"), 
                  p_12+theme(legend.position="none"), 
                  legend,
                  labels=c("a", "b","c","d","e",""), ncol = 3, nrow = 2, label_fontface = "bold",label_size=25,align = "hv")
# ggsave(".../PLOTS/figure1.pdf", units="in", width=14, height=9, dpi=600)
# ggsave(".../PLOTS/figure1.png", units="in", width=14, height=9, dpi=600)
plot(prow)


prow <- plot_grid(m + theme(legend.position="none"),
                  p_12+theme(legend.position="none"),
                  t + theme(legend.position="none"), 
                  n + theme(legend.position="none"), 
                  legend,
                  labels=c("a", "b","c","d",""), ncol = 5, nrow = 1, label_fontface = "bold",label_size=25,align = "hv")

plot(prow)
