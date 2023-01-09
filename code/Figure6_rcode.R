# volcano plots for figure 3

setwd("/home/chrats/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/4_Diff_expression/Diff_Expression_ST_3/tables")
minus.5614vsAll.complete <- read.delim("~/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/4_Diff_expression/Diff_Expression_ST_3/tables/minus 5614vsAll.complete.txt")
minus.6086vsAll.complete <- read.delim("~/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/4_Diff_expression/Diff_Expression_ST_3/tables/minus 6086vsAll.complete.txt")
minus.10675vsAll.complete <- read.delim("~/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/4_Diff_expression/Diff_Expression_ST_3/tables/minus 10675vsAll.complete.txt")
minus.SICOvsAll.complete <- read.delim("~/Desktop/Git/cheese-ft/Data/T0/metatranscriptomics/4_Diff_expression/Diff_Expression_ST_3/tables/minus SICOvsAll.complete.txt")



library(ggplot2)
library(ggrepel)
library(grid)
res_data$sig<-colr
size_text<-20

data<-minus.6086vsAll.complete
data$pvalue.minus.log10<- -log10(data$pvalue)
sig_temp<-rep(0,nrow(data))
sig_temp[which(data$padj<0.01)]<-1
data$sig<-factor(sig_temp)

#colors_t0<-c('#a50026','#1b7837','#fdae61','#4575b4','#313695','#abd9e9')

de6086<-ggplot(data,aes(x=log2FoldChange, y=pvalue.minus.log10, color=sig)) +scale_color_manual("",values=c('grey','#4575b4'))+
  geom_point(aes(size=sig)) +scale_size_manual(values=c(1,3))+
  #geom_text_repel(data=subset(res_data,sig==1),aes(label=names),parse  = TRUE,box.padding = 0.2)+
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  labs(title="-lactis(6086)vsAll")+
  theme(text = element_text(size=size_text),
        plot.title = element_text(size = rel(1), hjust = 0.5),
        axis.title.x = element_text(size=size_text ),
        axis.title.y = element_text( size=size_text),legend.title = element_text(size=size_text),
        legend.text = element_text(size_text))+guides(size = FALSE,color=F)+ 
  geom_text(label=paste("Total D.E. features: ", length(which(data$padj<0.01)),sep = ""), x=-3,  y=6, hjust=0, color='#4575b4', size=7)
de6086

######################
# ggsave("DE.pdf", units="in", width=7, height=5, dpi=300)
#######################

data<-minus.5614vsAll.complete
data$pvalue.minus.log10<- -log10(data$pvalue)
sig_temp<-rep(0,nrow(data))
sig_temp[which(data$padj<0.01)]<-1
data$sig<-factor(sig_temp)

de5614<-ggplot(data,aes(x=log2FoldChange, y=pvalue.minus.log10, color=sig)) +scale_color_manual("",values=c('grey','#fdae61'))+
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
de5614


############################

data<-minus.10675vsAll.complete
data$pvalue.minus.log10<- -log10(data$pvalue)
sig_temp<-rep(0,nrow(data))
sig_temp[which(data$padj<0.01)]<-1
data$sig<-factor(sig_temp)

de10675<-ggplot(data,aes(x=log2FoldChange, y=pvalue.minus.log10, color=sig)) +scale_color_manual("",values=c('grey','#313695'))+
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
de10675
######################

data<-minus.SICOvsAll.complete
data$pvalue.minus.log10<- -log10(data$pvalue)
sig_temp<-rep(0,nrow(data))
sig_temp[which(data$padj<0.01)]<-1
data$sig<-factor(sig_temp)


deSICO<-ggplot(data,aes(x=log2FoldChange, y=pvalue.minus.log10, color=sig)) +scale_color_manual("",values=c('grey','#abd9e9'))+
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
deSICO

prow <- plot_grid(de10675  + theme(legend.position="none"),
                  de6086 + theme(legend.position="none"),
                  de5614 + theme(legend.position="none"),
                  deSICO + theme(legend.position="none") ,
                  labels=c("b", "c","d","e"), ncol = 2, nrow =2, label_fontface = "bold",label_size=25,align = "hv")

prow
ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/figure6A.pdf",prow, units="in", width=9, height=9, dpi=600)
