#Collapse table using Qiime2
Check this: https://github.com/Anegin24/Microbiome-analysis/blob/main/Grouped-feature-table.docx
#Cal frequencies using Qiime2
qiime taxa collapse \
--i-table table.qza \
--o-collapsed-table collapse.table.qza \
--p-level 6 \
--i-taxonomy taxonomy.qza
qiime feature-table relative-frequency \
--i-table collapse.table.qza \
--o-relative-frequency-table collapse.frequency.table.qza \
--output-dir collapse.frequency/
  qiime tools export \
--input-path collapse.frequency.table.qza \
--output-path collapse.frequency/
biom convert 
--input-fp collapsegroup.frequency/feature-table.biom 
--output-fp collapsegroup.frequency.table.txt 
--header-key "taxonomy" 
--to-tsv
#Library
library(qiime2R)
library(tidyverse)
library(ggplot2)
library(ggtext)
library(glue)
library(RColorBrewer)
library(yaml)
library(Readr)
library(readxl)
#Add ASVs column (ASV1,ASV2,...) -> save as freq.tsv
taxonomy <- read_tsv("freqgroup.tsv") %>%
  separate(Taxon,
  into=c("kingdom", "phylum", "class", "order", "family", "genus"),sep=";") %>%
  mutate(genus = str_replace(string=genus,
                             pattern="(.*)", 
                             replacement="*\\1*"),
         taxon = glue("{genus}<br>({ASVs})"))
taxonomy<-taxonomy%>%select(-kingdom,-phylum,-class,-order,-family,-genus)
good<-taxonomy%>%filter(ASVs=="ASV71"|ASVs=="ASV200"|ASVs=="ASV115")
goodplot<-good%>%pivot_longer(c("Antibiotic","Control","Heal","Bacimix","Fermented-Product"))%>%mutate(Abundance=value*100)
ggplot(goodplot,aes(x=name,y=Abundance,fill=taxon)) + 
  geom_col() + 
  labs(x="Experimental Group",y="Relative Abundance (%)", title = "Percentages of Good Bacteria in Experimental groups") + 
  scale_y_continuous(expand=c(0,0),limits = c(0,30), breaks = seq(0, 30, by=5)) + 
  theme_classic() + 
  theme (plot.title = element_textbox_simple(size=20,margin=margin(0,0,0,25))) +
  scale_fill_manual(name=NULL, values = c(brewer.pal(4,"Dark2"),"gray"))
ggsave("goodplot.pdf", height=4, width=8, device="pdf")
#Do the same with Harmful bacteria