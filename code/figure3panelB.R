# Cheese paper figure 3 panel B

# set wd
setwd("/home/chrats/Desktop/Git/cheese-ft/Data/Models")

# load libs
library(tidyverse)

# load data, join metabolite names and species taxonomy with extracellular fluxes
metnames=read.delim("data/bigg_classes.tsv",header=TRUE)
species = read.delim("data/species_names.tsv",header = TRUE)
ex_flux = read.delim("data/individual_fluxes.tsv",header=TRUE) %>% 
  left_join(species,ex_flux,by="species") %>% 
  mutate(direction=ifelse(flux>0,"export","import")) %>% 
  left_join(.,metnames,by="reaction") %>%
  mutate(name=ifelse(name=="NA",reaction,name)) %>% filter(name!="(±)-Tryptophan",
                                                           name!="1-deoxy-1-(N6-lysino)-D-fructose",
                                                           name!="N-acetyl-seryl-aspartate",
                                                           name!="L-lysinium(1+)",
                                                           name!="L-argininium(1+)")

ex_flux$name[which(ex_flux$name=="gamma-Aminobutyric acid")] <- "GABA"
# plot data, play around with geom_dotplot dotsize parameter if dots are too small after re-sizing
fig4aold_nowb<-ex_flux %>% 
  filter(name %in% c("Formic acid","Glycolic acid","D-Galactose","Ammonia", "D-Alanine", "Glycine", "L-Serine", "L-Valine", "Citric acid", "Folic acid", "DL-Glutamate", "Diacetyl", "L-Glutamine", "Lactose" , "Succinic acid", "GABA", "(S)-Propane-1,2-diol","Fumaric acid")) %>% 
  ggplot(aes(x=as.factor(name),fill=direction)) + scale_fill_manual(values= c("#1b9e77","#d95f02"))+
  geom_dotplot(method="dotdensity",dotsize=.9,stackgroups = TRUE) + coord_flip()+ 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "bottom") + 
  facet_wrap(~taxonomy,scales="free_x",nrow=1) +ylab("Frequency") + xlab("Exchanged metabolite") + labs(fill = "Reaction directionality") + 
  theme(legend.title=element_text(size=20),legend.text=element_text(size=20),axis.text = element_text(size=20),
        strip.text.x = element_text(size = 20),axis.title = element_text(size=20),axis.text.x = element_blank())

ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/FINAL_4main/fig3panelB.pdf",height = 6, width =10)
