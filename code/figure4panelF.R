# Cheese paper figure 4 panel f

# set wd
setwd("/home/chrats/Desktop/Git/cheese-ft/Data/Models")

# load libs
library(tidyverse)

# load data, join metabolite names and species taxonomy with extracellular fluxes
bigg_rxns=read.delim("data/bigg_models_reactions.txt")
species = read.delim("data/species_names.tsv",header = TRUE)
fba=read.delim("data/summary_fba.tsv",header=FALSE) %>% dplyr::rename(.,model=V1,gf_media=V2,sim_media=V3,rxn=V4,flux=V5)

# identify all reactions of interest
all_rxns=c("R_VALTA","R_VALt2r","R_EX_val__L_e","R_GLUDC","R_GLUDy","R_GLUSy","R_GLUSx","R_GLUABUTt7","R_GLNS","R_ICDHyr","R_ACONT","R_CS","R_PC","R_CITt7","R_BTDD_RR","R_ACTD","R_ACTD2","R_ACTDa","R_ACLDC","R_ACLS","R_KARA1i","R_KARA1")
colors_t0s<-c('#fdae61','#313695','#4575b4','#1b7837')

# plot data
data<-fba %>% 
  filter(rxn %in% all_rxns) %>% 
  left_join(.,species%>%dplyr::rename(model=species)) %>% 
  drop_na() %>% filter(gf_media!="no_gf") %>%
  mutate(rxn=gsub("R_","",rxn),rxn=gsub("_e","",rxn)) %>%
  mutate(sim_media=gsub("milk_aer","Milk (aerobic)",sim_media),sim_media=gsub("milk_ana","Milk (anaerobic)",sim_media),sim_media=gsub("mm_aer","Minimal milk (aerobic)",sim_media),sim_media=gsub("mm_ana","Minimal milk (anaerobic)",sim_media)) %>%
  mutate(gf_media=gsub("milk_aer","Milk (aerobic)",gf_media),gf_media=gsub("milk_ana","Milk (anaerobic)",gf_media),gf_media=gsub("mm_aer","Minimal milk (aerobic)",gf_media),gf_media=gsub("mm_ana","Minimal milk (anaerobic)",gf_media))
data$taxonomy<-gsub("CHCC","",data$taxonomy)


data_sel<-data[data$rxn%in%c("ICDHyr","ACONT","CS","PC" , "ACLS" ,"ACLDC","KARA1" ,"KARA1i","ACTD","ACTDa" ,"BTDD_RR" ),]
data_sel$rxn<-factor(data_sel$rxn,levels = c("ICDHyr","ACONT","CS","PC" , "ACLS" ,"ACLDC","KARA1" ,"KARA1i","ACTD","ACTDa" ,"BTDD_RR" ))
ggplot(data_sel) + 
  geom_boxplot(aes(x=taxonomy,y=abs(flux),fill=taxonomy),alpha=0.6,outlier.shape = NA) + 
  geom_jitter(aes(x=taxonomy,y=abs(flux),color=taxonomy,shape=taxonomy),alpha=0.5,size=5) + 
  scale_y_log10() + facet_wrap(~rxn,nrow = 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("|Flux| (mmol/gDCW*hr)") + xlab("Species") +
  # labs(color="Simulation media",shape="Gap-filling media") + 
  labs(fill="Taxonomy") + 
  theme(legend.title=element_text(size=20),legend.text=element_text(size=20),axis.text = element_text(size=20),
        strip.text.x = element_text(size = 20),axis.title = element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values = colors_t0s)+scale_color_manual(values = colors_t0s)

ggsave("/home/chrats/Desktop/Git/cheese-ft/PLOTS/FINAL_4main/fig4f.pdf",height = 5,width = 25)

