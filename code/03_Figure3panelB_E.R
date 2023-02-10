# volcano plots for figure 3 b-e
setwd("~/Desktop/Git/")
minus.LCvsAll.complete <- read.delim("./mi_cheese/data/T0/metatranscriptomics/4_Diff_expression/Diff_Expression_ST_3/tables/minus LCvsAll.complete.txt")
minus.LLm2vsAll.complete <- read.delim("./mi_cheese/data/T0/metatranscriptomics/4_Diff_expression/Diff_Expression_ST_3/tables/minus LLm2vsAll.complete.txt")
minus.LLm1vsAll.complete <- read.delim("./mi_cheese/data/T0/metatranscriptomics/4_Diff_expression/Diff_Expression_ST_3/tables/minus LLm1vsAll.complete.txt")
minus.LBvsAll.complete <- read.delim("./mi_cheese/data/T0/metatranscriptomics/4_Diff_expression/Diff_Expression_ST_3/tables/minus LBvsAll.complete.txt")

library(ggplot2)
library(ggrepel)
library(grid)

size_text<-20

data<-minus.LLm2vsAll.complete
data$pvalue.minus.log10<- -log10(data$pvalue)
sig_temp<-rep(0,nrow(data))
sig_temp[which(data$padj<0.01)]<-1
data$sig<-factor(sig_temp)

#colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')

deLLm2<-ggplot(data,aes(x=log2FoldChange, y=pvalue.minus.log10, color=sig)) +scale_color_manual("",values=c('grey','#4575b4'))+
  geom_point(aes(size=sig)) +scale_size_manual(values=c(1,3))+
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  labs(title="-lactis(6086)vsAll")+
  theme(text = element_text(size=size_text),
        plot.title = element_text(size = rel(1), hjust = 0.5),
        axis.title.x = element_text(size=size_text ),
        axis.title.y = element_text( size=size_text),legend.title = element_text(size=size_text),
        legend.text = element_text(size_text))+guides(size = FALSE,color=F)+ 
  geom_text(label=paste("Total D.E. features: ", length(which(data$padj<0.01)),sep = ""), x=-3,  y=6, hjust=0, color='#4575b4', size=7)
deLLm2

######################
# ggsave("DE.pdf", units="in", width=7, height=5, dpi=300)
#######################

data<-minus.LCvsAll.complete
data$pvalue.minus.log10<- -log10(data$pvalue)
sig_temp<-rep(0,nrow(data))
sig_temp[which(data$padj<0.01)]<-1
data$sig<-factor(sig_temp)

deLC<-ggplot(data,aes(x=log2FoldChange, y=pvalue.minus.log10, color=sig)) +scale_color_manual("",values=c('grey','#fdae61'))+
  geom_point(aes(size=sig),shape=17) +scale_size_manual(values=c(1,3))+
  #geom_text_repel(data=subset(res_data,sig==1),aes(label=names),parse  = TRUE,box.padding = 0.2)+
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  labs(title="-cremoris(5614)vsAll")+
  theme(text = element_text(size=size_text),
        plot.title = element_text(size = rel(1), hjust = 0.5),
        axis.title.x = element_text(size=size_text ),
        axis.title.y = element_text( size=size_text),legend.title = element_text(size=size_text),
        legend.text = element_text(size_text))+guides(size = FALSE,color=F)+ 
  geom_text(label=paste("Total D.E. features: ", length(which(data$padj<0.01)),sep = ""), x=-3,  y=15, hjust=0, color='#fdae61', size=7)
deLC


############################

data<-minus.LLm1vsAll.complete
data$pvalue.minus.log10<- -log10(data$pvalue)
sig_temp<-rep(0,nrow(data))
sig_temp[which(data$padj<0.01)]<-1
data$sig<-factor(sig_temp)

deLLm1<-ggplot(data,aes(x=log2FoldChange, y=pvalue.minus.log10, color=sig)) +scale_color_manual("",values=c('grey','#313695'))+
  geom_point(aes(size=sig)) +scale_size_manual(values=c(1,3))+
  #geom_text_repel(data=subset(res_data,sig==1),aes(label=names),parse  = TRUE,box.padding = 0.2)+
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  labs(title="-lactis(10675)vsAll")+
  theme(text = element_text(size=size_text),
        plot.title = element_text(size = rel(1), hjust = 0.5),
        axis.title.x = element_text(size=size_text ),
        axis.title.y = element_text( size=size_text),legend.title = element_text(size=size_text),
        legend.text = element_text(size_text))+guides(size = FALSE,color=F)+ 
  geom_text(label=paste("Total D.E. features: ", length(which(data$padj<0.01)),sep = ""), x=-4,  y=70, hjust=0, color='#313695', size=7)
deLLm1
######################

data<-minus.LBvsAll.complete
data$pvalue.minus.log10<- -log10(data$pvalue)
sig_temp<-rep(0,nrow(data))
sig_temp[which(data$padj<0.01)]<-1
data$sig<-factor(sig_temp)


deLB<-ggplot(data,aes(x=log2FoldChange, y=pvalue.minus.log10, color=sig)) +scale_color_manual("",values=c('grey','#abd9e9'))+
  geom_point(aes(size=sig),shape=15) +scale_size_manual(values=c(1,3))+
  #geom_text_repel(data=subset(res_data,sig==1),aes(label=names),parse  = TRUE,box.padding = 0.2)+
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  labs(title="-SICOvsAll")+
  theme(text = element_text(size=size_text),
        plot.title = element_text(size = rel(1), hjust = 0.5),
        axis.title.x = element_text(size=size_text ),
        axis.title.y = element_text( size=size_text),legend.title = element_text(size=size_text),
        legend.text = element_text(size_text))+guides(size = FALSE,color=F)+ 
  geom_text(label=paste("Total D.E. features: ", length(which(data$padj<0.01)),sep = ""), x=-3,  y=18, hjust=0, color='#abd9e9', size=7)
deLB

prow <- plot_grid(deLLm1  + theme(legend.position="none"),
                  deLLm2 + theme(legend.position="none"),
                  deLC + theme(legend.position="none"),
                  deLB + theme(legend.position="none") ,
                  labels=c("b", "c","d","e"), ncol = 2, nrow =2, label_fontface = "bold",label_size=25,align = "hv")

prow
#ggsave("/PLOTS/figure3be.pdf",prow, units="in", width=9, height=9, dpi=600)
