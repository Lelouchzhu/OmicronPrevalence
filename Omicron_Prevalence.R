#devtools::install_github("outbreak-info/R-outbreak-info")

library(outbreakinfo)
library(ggbreak)
library(tidyverse)
library(ggpubr)
outbreakinfo::authenticateUser()
#  Provide GISAID credentials using authenticateUser()


locations = c("Singapore","United States", "South Africa")


Omicron_sub<-c("BA.1","BA.2","BA.2.12.1","BA.4","BA.5")

Omicron_labels = Omicron_sub
BA1_lineages = lookupSublineages("BA.1", returnQueryString = TRUE)
BA2_lineages = lookupSublineages("BA.2", returnQueryString = TRUE)
BA2121_lineages = lookupSublineages("BA.2.12.1", returnQueryString = TRUE)
BA4_lineages = lookupSublineages("BA.4", returnQueryString = TRUE)
BA5_lineages = lookupSublineages("BA.5", returnQueryString = TRUE)

names(Omicron_labels) = c(BA1_lineages, BA2_lineages, BA2121_lineages,BA4_lineages, BA5_lineages)

Omicron_sub_by_state = purrr::map_df(locations, function(loc) getPrevalence(pangolin_lineage =  c(BA1_lineages, BA2_lineages,BA2121_lineages, BA4_lineages, BA5_lineages), location = loc))

Omicron_sub_world = getPrevalence(pangolin_lineage =  c(BA1_lineages, BA2_lineages,BA2121_lineages, BA4_lineages, BA5_lineages) )

plotPrevalenceOverTime(rbind(Omicron_sub_by_state,Omicron_sub_world) %>% filter(date > "2021-11-01"), labelDictionary = Omicron_labels)+
  facet_grid(rows=vars(location))