# Cheese paper figure 4 panel c showing smetana interactions
setwd("~/Desktop/Git/")
# load libs
library(tidyverse)
library(ggalluvial)

# load data, join metabolite names and species taxonomy with extracellular fluxes
milk_exp = read.delim("./mi_cheese/data/Models/data/interactions.tsv") %>% mutate(compound=gsub("M_","",compound)) %>% mutate(compound=gsub("_e","",compound))
smetana_species = read.delim("./mi_cheese/data/Models/data/smetana_taxonomy.tsv",header = TRUE)
metnames=read.delim("./mi_cheese/data/Models/data/bigg_classes.tsv",header=TRUE) %>% dplyr::rename(.,compound=reaction)

# plot data
ggplot(milk_exp %>% filter(smetana>0.25) %>% left_join(.,smetana_species%>%dplyr::select(donor,donor_species),by="donor") %>% left_join(.,smetana_species%>%dplyr::select(receiver,receiver_species),by="receiver") %>% left_join(.,metnames,by="compound") %>% filter(name!="Ammonia") %>% group_by(receiver_species,donor_species,name) %>% dplyr::summarize(ave_smet=mean(smetana),count=n()) %>% mutate(count=as.factor(count)) ,
       aes(axis1 = donor_species, axis2 = name, axis3 = receiver_species,
           y = ave_smet)) +
  scale_x_discrete(limits = c("Donor", "Metabolite", "Reciever")) +
  xlab("Interaction") +
  geom_alluvium(aes(fill = count)) +
  geom_stratum(width=0.4) +
  theme_minimal() + geom_text(stat = "stratum", aes(label = after_stat(stratum)),min.y=0.2)+theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "bottom") + scale_fill_discrete(breaks=c("1","2","3")) + labs(fill = "Counts across simulation conditions")

# ggsave("figures/fig4panelC.pdf",height = 6, width =12)