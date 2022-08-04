library(outbreakinfo)
library(ggbreak)
library(tidyverse)
library(ggpubr)
#outbreakinfo::authenticateUser() #Provide GISAID credentials using authenticateUser()

BA275_india = getPrevalence(pangolin_lineage = "BA.2.75", location = "India")
BA275_us = getPrevalence(pangolin_lineage = "BA.2.75", location = "United States")
BA275_uk = getPrevalence(pangolin_lineage = "BA.2.75", location = "United Kingdom")
BA275_jp = getPrevalence(pangolin_lineage = "BA.2.75", location = "Japan")
BA275_ca = getPrevalence(pangolin_lineage = "BA.2.75", location = "Canada")
BA275_au = getPrevalence(pangolin_lineage = "BA.2.75", location = "Australia")

BA275 = dplyr::bind_rows(BA275_india, BA275_us, BA275_uk,BA275_jp,BA275_ca,BA275_au)
plotPrevalenceOverTime(BA275, colorVar = "location", title="BA.2.75 prevalence over time")

ggsave(file="Extend Fig 1_July22_1.pdf",device="pdf",width=10,height=8)
