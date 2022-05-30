#  Omicron Prevalence

#devtools::install_github("outbreak-info/R-outbreak-info")

library(outbreakinfo)
library(ggbreak)
library(tidyverse)
library(ggpubr)
outbreakinfo::authenticateUser() #Provide GISAID credentials using authenticateUser()


locations = c("Israel","Singapore","United States", "South Africa")


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

ggsave(file="Extend Fig 1.pdf",device="pdf",width=10,height=8)

#  Antigenic cartography

library(tidyverse)
library(openxlsx)
library(readr)
library(broom)
library(plotly)
library(corrplot)
library(ggpubr)
library(Racmacs)

##  For naive sera
names_df_pvnt<-read.xlsx("pVNT 250522 Updated.xlsx",sheet = "P1")%>%
  select(Serum.panel,ID)
unique(names_df_pvnt$Serum.panel)

pvnt_df<-read.xlsx("pVNT 250522 Updated.xlsx",sheet = "P1")%>%
  select(-Serum.panel)%>%
  #filter(ID!="90001075")%>%
  pivot_longer(-ID,names_to = "Virus",values_to = "NT50") %>% 
  mutate(NT50=round(NT50))%>%
  pivot_wider(names_from=ID, values_from=NT50) %>%
  mutate(Virus=gsub("\\.", " ", Virus))%>%
  mutate(Virus=gsub("BA ", "BA.", Virus))
pvnt_df2<-pvnt_df%>%select(-Virus)
rownames(pvnt_df2)<-pvnt_df$Virus
pvnt_df2<-as.matrix(pvnt_df2)

rownames(pvnt_df2)

map_pvnt <- acmap(
  titer_table = pvnt_df2,
  ag_names = rownames(pvnt_df2),
  sr_names = colnames(pvnt_df2),
)

map_pvnt <- optimizeMap(
  map                     = map_pvnt,
  number_of_dimensions    = 2,
  number_of_optimizations = 1000,
  minimum_column_basis    = "none",
)
view(map_pvnt)

map_pvnt

map_pvnt2<-relaxMap(map_pvnt)
view(map_pvnt2)

plotly_map_table_distance(map_pvnt2)
map_pvnt2 <- bootstrapMap(
  map                      = map_pvnt,
  method = "resample",
  bootstrap_repeats        = 1000,
  optimizations_per_repeat = 100,
  ag_noise_sd              = 0.7,
  titer_noise_sd           = 0.7
)

boostrap_ag_coords_list <- mapBootstrap_agCoords(map_pvnt2)
boostrap_sr_coords_list <- mapBootstrap_srCoords(map_pvnt2)

agGroups(map_pvnt2)<-c(2,3,3,3,1)
srGroups(map_pvnt2)<-names_df_pvnt$Serum.panel
plotly_map_table_distance(map_pvnt2)

view(map_pvnt2)

pdf(file = "Fig 2e.pdf", width = 18, height = 14)

plot(
  map_pvnt2,
  optimization_number = 1,
  xlim = c(-4,5),  #Adjust the xlim and ylim to plot the figure properly
  ylim = c(-4,3),  #Adjust the xlim and ylim to plot the figure properly
  plot_ags = TRUE,
  plot_sr = TRUE,
  plot_labels = "antigens",
  plot_blobs = TRUE,
  show_procrustes = TRUE,
  show_error_lines = FALSE,
  plot_stress = FALSE,
  indicate_outliers = "arrowheads",
  grid.col = "grey90",
  grid.margin.col = "grey50",
  outlier.arrow.col = grid.col,
  fill.alpha = 0.8,
  outline.alpha = 0.8,
  label.offset = 0,
  padding = 1,
  cex = 1,
)

dev.off()

##  For vaccinated and breakthrough sera

names_df_vaccine_pvnt<-read.xlsx("pvnt 250522 Updated.xlsx",sheet = "P2")%>%
  select(Serum.panel,ID)
unique(names_df_vaccine_pvnt$Serum.panel)

vaccine_pvnt_df<-read.xlsx("pvnt 250522 Updated.xlsx",sheet = "P2")%>%
  select(-Serum.panel)%>%
  #filter(ID!="90001075")%>%
  pivot_longer(-ID,names_to = "Virus",values_to = "NT50") %>% 
  mutate(NT50=round(NT50))%>%
  pivot_wider(names_from=ID, values_from=NT50) %>%
  mutate(Virus=gsub("\\.", " ", Virus))%>%
  mutate(Virus=gsub("BA ", "BA.", Virus))
vaccine_pvnt_df2<-vaccine_pvnt_df%>%select(-Virus)
rownames(vaccine_pvnt_df2)<-vaccine_pvnt_df$Virus
vaccine_pvnt_df2<-as.matrix(vaccine_pvnt_df2)

rownames(vaccine_pvnt_df2)

map_vaccine_pvnt <- acmap(
  titer_table = vaccine_pvnt_df2,
  ag_names = rownames(vaccine_pvnt_df2),
  sr_names = colnames(vaccine_pvnt_df2),
)

map_vaccine_pvnt <- optimizeMap(
  map                     = map_vaccine_pvnt,
  number_of_dimensions    = 2,
  number_of_optimizations = 1000,
  minimum_column_basis    = "none",
)
view(map_vaccine_pvnt)

map_vaccine_pvnt

map_vaccine_pvnt2<-relaxMap(map_vaccine_pvnt)
view(map_vaccine_pvnt2)

plotly_map_table_distance(map_vaccine_pvnt2)
map_vaccine_pvnt2 <- bootstrapMap(
  map                      = map_vaccine_pvnt,
  method = "resample",
  bootstrap_repeats        = 1000,
  optimizations_per_repeat = 100,
  ag_noise_sd              = 0.7,
  titer_noise_sd           = 0.7
)

boostrap_ag_coords_list <- mapBootstrap_agCoords(map_vaccine_pvnt2)
boostrap_sr_coords_list <- mapBootstrap_srCoords(map_vaccine_pvnt2)

agGroups(map_vaccine_pvnt2)<-c(2,3,3,3,3,3,1)
srGroups(map_vaccine_pvnt2)<-names_df_vaccine_pvnt$Serum.panel
plotly_map_table_distance(map_vaccine_pvnt2)

view(map_vaccine_pvnt2)
pdf(file = "Extend Fig 5.pdf", width = 20, height = 22)
plot(
  map_vaccine_pvnt2,
  optimization_number = 1,
  xlim = c(-5,5),
  ylim = c(-6,6),
  plot_ags = TRUE,
  plot_sr = TRUE,
  plot_labels = "antigens",
  plot_blobs = TRUE,
  show_procrustes = TRUE,
  show_error_lines = FALSE,
  plot_stress = FALSE,
  indicate_outliers = "arrowheads",
  grid.col = "grey90",
  grid.margin.col = "grey50",
  outlier.arrow.col = grid.col,
  fill.alpha = 0.8,
  outline.alpha = 0.8,
  label.offset = 0,
  padding = 1,
  cex = 1,
)
dev.off()
