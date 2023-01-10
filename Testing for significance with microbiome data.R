#Cal frequencies using Qiime2
qiime taxa collapse \
--i-table table.qza \
--o-collapsed-table collapse.table.qza \
--p-level 6 \
--i-taxonomy taxonomy.qza
qiime feature-table relative-frequency \
--i-table collapse.table.qza \
--o-relative-frequency-table frequency.table.qza \
--output-dir collapse.frequency/
  qiime tools export \
--input-path frequency.table.qza \
--output-path frequency/
  biom convert \
--input-fp frequency/feature-table.biom \
--output-fp freqsam.txt \
--header-key "taxonomy" \
--to-tsv
#Library
library(qiime2R)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggtext)
library(glue)
library(RColorBrewer)
library(yaml)
library(readr)
library(readxl)
#transpose data in freqsam.tsv, add sampleid, rename genus to ASV1,ASV2,...ASVn, save as freq.tsv
freq <- read_tsv("significant differential/freq.tsv") %>%
rename_all(tolower) %>%
select(sampleid, starts_with("asvs")) %>%
pivot_longer(-sampleid, names_to="asvs", values_to="count")
#freqsam.tsv: add Taxon header, add column ASVs, save as taxonomy.tsv
taxonomy <- read_tsv("significant differential//taxonomy.tsv") %>%
separate(Taxon,
        into=c("kingdom", "phylum", "class", "order", "family", "genus"),sep=";") %>%
        mutate(genus = str_replace(string=genus,
        pattern="(.*)", 
        replacement="*\\1*"),
        taxon = glue("{genus} ({ASVs})"))
#Generate composite dataframe:
composite <- inner_join(freq,taxonomy,by="ASVs") %>% 
mutate (abundance = count*100) %>% 
select (-count) %>% 
inner_join(.,metadata, by = c("sampleid"="SampleID"))
#Statistical test: Kruskall wallis or Wilcoxon
sig<-composite %>% 
nest(data = -taxon) %>% 
mutate(test = map(.x=data, ~kruskal.test(abundance~Type,data=.x) %>% tidy)) %>% 
unnest(test) %>% 
filter (p.value < 0.05) %>% 
select (taxon,p.value)
#merge sig with composite and plot:
composite %>%
inner_join(sig,by="taxon") %>% 
ggplot(aes(x=abundance,y=taxon,fill=Type)) + 
geom_col(position = position_dodge())

