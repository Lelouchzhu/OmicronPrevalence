library(outbreakinfo)
library(ggbreak)
library(tidyverse)
library(ggpubr)
#outbreakinfo::authenticateUser() #Provide GISAID credentials using authenticateUser()


locations = c("Israel","Singapore","United States", "South Africa")

Omicron_sub<-c("BA.1","BA.2","BA.5","BA.2.75","BQ.1.1","XBB","XBB.1.5")

Omicron_labels = Omicron_sub
BA1_lineages = lookupSublineages("BA.1", returnQueryString = TRUE)
BA2_lineages = lookupSublineages("BA.2", returnQueryString = TRUE)
BA5_lineages = lookupSublineages("BA.5", returnQueryString = TRUE)
BA275_lineages = lookupSublineages("BA.2.75", returnQueryString = TRUE)
BQ11_lineages = lookupSublineages("BQ.1.1", returnQueryString = TRUE)
XBB_lineages = lookupSublineages("XBB", returnQueryString = TRUE)
XBB15_lineages = lookupSublineages("XBB.1.5", returnQueryString = TRUE)

names(Omicron_labels) = c(BA1_lineages, BA2_lineages,BA5_lineages,BA275_lineages,BQ11_lineages,XBB_lineages,XBB15_lineages)

Omicron_sub_by_state = purrr::map_df(locations, function(loc) getPrevalence(pangolin_lineage =  c(BA1_lineages, BA2_lineages,BA5_lineages,BA275_lineages,BQ11_lineages,XBB_lineages,XBB15_lineages), location = loc))

Omicron_sub_world = getPrevalence(pangolin_lineage =  c(BA1_lineages, BA2_lineages,BA5_lineages,BA275_lineages,BQ11_lineages,XBB_lineages,XBB15_lineages) )

plotPrevalenceOverTime(rbind(Omicron_sub_by_state,Omicron_sub_world) %>% filter(date > "2021-12-01",date < "2023-03-28"), labelDictionary = Omicron_labels)+
  facet_grid(rows=vars(location))

ggsave(file="2023/Fig 1A_Mar28.pdf",device="pdf",width=10,height=8)



# Plot the mutations as a heatmap
mutations = getMutationsByLineage(pangolin_lineage=Omicron_sub, frequency=0.75, logInfo = FALSE)
plotMutationHeatmap(mutations, title = "S-gene mutations in lineages")
ggsave(file="2023/Fig1B_Mar28.pdf",device="pdf",width=12,height=4)

#  Stats on Fig2
library(rstatix)
library(openxlsx)

pvnt_df<-read.xlsx("2023/pVNT_serotype_202208.xlsx",sheet = "Sheet1")%>%
  pivot_longer(-(Serum.panel:ID),names_to = "Virus",values_to = "NT50")%>%
  mutate(logNT50=log10(NT50))%>%
  mutate(Virus=gsub("\\.", " ", Virus))%>%
  mutate(Virus=gsub("BA ", "BA.", Virus))
virus_filter3<-c("Ancestral","Delta","Beta","BA.1","BA.2","BA.5","SARS-CoV-1")

Group_filter3<-c("Ancestral","Delta","Beta","BA.1","BA.2","BA.5","SARS1")

pvnt_df$Virus <- factor(pvnt_df$Virus, levels=virus_filter3)
pvnt_df$Serum.panel <- factor(pvnt_df$Serum.panel, levels=Group_filter3)

stat.test <- pvnt_df %>%
  group_by(Serum.panel) %>%
  wilcox_test(logNT50 ~ Virus) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

write.xlsx(stat.test,file = "2023/Fig2_p_value.xlsx")

#  Stats on Fig3
pvnt_df2<-read.xlsx("2023/For antigenic cartograph_VAX_20230327.xlsx",sheet = "df")%>%
  pivot_longer(-(Serum.panel:ID),names_to = "Virus",values_to = "NT50")%>%
  mutate(logNT50=log10(NT50))%>%
  mutate(Virus=gsub("\\.", " ", Virus))%>%
  mutate(Virus=gsub("BA ", "BA.", Virus))%>%
  mutate(Virus=gsub("BQ ", "BQ.", Virus))

virus_filter3<-c("SARS-CoV-2","BA.1","BA.2",
                 "BA.2 L452R",
                 "BA.2 F486V",
                 "BA.5",
                 "BA.4.6.1",
                 "BA.2.75.2",
                 "XBB",
                 "BQ.1.1",
                 "SARS-CoV-1"
)

Group_filter3<-c("PP","PPP")

pvnt_df2$Virus <- factor(pvnt_df2$Virus, levels=virus_filter3)
pvnt_df2$Serum.panel <- factor(pvnt_df2$Serum.panel, levels=Group_filter3)

stat.test2 <- pvnt_df2 %>%
  group_by(Serum.panel) %>%
  pairwise_wilcox_test(logNT50 ~ Virus) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

write.xlsx(stat.test2,file = "2023/Fig3_p_value.xlsx")

#  Antigenic cartography

library(openxlsx)
library(readr)
library(broom)
library(Racmacs)

##  For naive sera
names_df_pvnt<-read.xlsx("2023/pVNT_serotype_202208.xlsx",sheet = "Sheet1")%>%
  select(Serum.panel,ID)
unique(names_df_pvnt$Serum.panel)

pvnt_df<-read.xlsx("2023/pVNT_serotype_202208.xlsx",sheet = "Sheet1")%>%
  select(-Serum.panel)%>%
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
  number_of_optimizations = 2000,
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
  optimizations_per_repeat = 400,
  ag_noise_sd              = 0.68,
  titer_noise_sd           = 0.55
)

boostrap_ag_coords_list <- mapBootstrap_agCoords(map_pvnt2)
boostrap_sr_coords_list <- mapBootstrap_srCoords(map_pvnt2)

agGroups(map_pvnt2)<-c(2,2,2,3,3,3,1)
srGroups(map_pvnt2)<-names_df_pvnt$Serum.panel
plotly_map_table_distance(map_pvnt2)

view(map_pvnt2)


# save.acmap(
#   map_pvnt2,
#   filename = "naive_serotype.ace",
#   compress = FALSE,
#   round_titers = FALSE
# )

# Load from previous generated acmap file
map_pvnt2 <- read.acmap("2023/naive_serotype.ace")
map_pvnt2_plot<-map_pvnt2

naive_ag_cords<-map_pvnt2_plot[["optimizations"]][[1]][["ag_base_coords"]]
rownames(naive_ag_cords)<-rownames(pvnt_df2)

naive_ag_cords<-as.data.frame(naive_ag_cords)%>%
  rename(X=V1,Y=V2)%>%
  mutate(Antigen=rownames(naive_ag_cords))

# Create a function to calculate the Euclidean distance between two points
euclidean_distance <- function(x1, y1, x2, y2) {
  return(sqrt((x1 - x2)^2 + (y1 - y2)^2))
}

# Calculate the pairwise distances between Antigens based on their X and Y coordinates
n <- nrow(naive_ag_cords)
distances <- matrix(0, n, n)

for (i in 1:n) {
  for (j in 1:n) {
    distances[i, j] <- euclidean_distance(naive_ag_cords$X[i], naive_ag_cords$Y[i], naive_ag_cords$X[j], naive_ag_cords$Y[j])
  }
}

# Convert the matrix to a data frame and label rows and columns
naive_distances_df <- as.data.frame(distances)
colnames(naive_distances_df)<-rownames(naive_ag_cords)
rownames(naive_distances_df)<-rownames(naive_ag_cords)
write.xlsx(naive_distances_df,"2023/Fig2I_naive_AU.xlsx",overwrite = T)


agFill(map_pvnt2_plot)<-c("blue","dark blue","light blue","red","dark red","pink","yellow")
test1<-c(rep("blue", 20),rep("light blue",10),rep("dark blue",10),rep("red",11),rep("dark red",10),rep("pink",5),rep("yellow",14))
srOutline(map_pvnt2_plot)<-test1
pdf(file = "2023/Fig 2H.pdf", width = 10, height = 9)

plot(
  map_pvnt2_plot,
  optimization_number = 1,
  xlim = c(-4,6),  #Adjust the xlim and ylim to plot the figure properly
  ylim = c(-4,5),  #Adjust the xlim and ylim to plot the figure properly
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

#---vaccinated samples
##  For PP samples
names_df_pvnt<-read.xlsx("2023/For antigenic cartograph_VAX_20230327.xlsx",sheet = "df")%>%
  select(Serum.panel,ID)%>%
  filter(Serum.panel=="PP")
unique(names_df_pvnt$Serum.panel)

pvnt_df<-read.xlsx("2023/For antigenic cartograph_VAX_20230327.xlsx",sheet = "df")%>%
  filter(Serum.panel=="PP")%>%
  select(-Serum.panel)%>%
  pivot_longer(-ID,names_to = "Virus",values_to = "NT50") %>% 
  mutate(NT50=round(NT50))%>%
  pivot_wider(names_from=ID, values_from=NT50) %>%
  mutate(Virus=gsub("\\.", " ", Virus))%>%
  mutate(Virus=gsub("BA ", "BA.", Virus))%>%
  mutate(Virus=case_when(Virus=="BA.4 6 1"~"BA.4.6.1",
                         Virus=="BA.2 75 2"~"BA.2.75.2",
                         Virus=="BQ 1 1"~"BQ.1.1",
                         T~Virus
  ))
pvnt_df2<-pvnt_df%>%select(-Virus)
rownames(pvnt_df2)<-pvnt_df$Virus
pvnt_df2<-as.matrix(pvnt_df2)

pvnt_df_test<-pvnt_df%>%pivot_longer(-Virus,names_to = "ID",values_to = "NT50") %>% 
  group_by(Virus)%>%
  summarise(sd=sd(log10(NT50)))

pvnt_df_test2<-pvnt_df%>%pivot_longer(-Virus,names_to = "ID",values_to = "NT50") %>% 
  group_by(ID)%>%
  summarise(sd=sd(log10(NT50)))


rownames(pvnt_df2)

map_pvnt <- acmap(
  titer_table = pvnt_df2,
  ag_names = rownames(pvnt_df2),
  sr_names = colnames(pvnt_df2),
)

map_pvnt <- optimizeMap(
  map                     = map_pvnt,
  number_of_dimensions    = 2,
  number_of_optimizations = 2000,
  minimum_column_basis    = "none",
)
view(map_pvnt)

map_pvnt

map_pvnt2<-relaxMap(map_pvnt)
view(map_pvnt2)

plotly_map_table_distance(map_pvnt2)
map_pvnt2 <- bootstrapMap(
  map                      = map_pvnt,
  method = "noisy",
  bootstrap_repeats        = 2000,
  optimizations_per_repeat = 400,
  ag_noise_sd              = 0.7,
  titer_noise_sd           = 0.7
)

boostrap_ag_coords_list <- mapBootstrap_agCoords(map_pvnt2)
boostrap_sr_coords_list <- mapBootstrap_srCoords(map_pvnt2)

agGroups(map_pvnt2)<-c(2,3,3,3,3,3,3,3,3,3,1)
srGroups(map_pvnt2)<-names_df_pvnt$Serum.panel
plotly_map_table_distance(map_pvnt2)

view(map_pvnt2)


# save.acmap(
#   map_pvnt2,
#   filename = "PP_serotype_Mar27_1.ace",
#   compress = FALSE,
#   round_titers = FALSE
# )

# Load from previous generated acmap file
map_pvnt2 <- read.acmap("2023/PP_serotype_Mar27_1.ace")
map_pvnt2_plot<-map_pvnt2

test1<-c(rep("light blue",20))

srOutline(map_pvnt2_plot)<-test1
pdf(file = "2023/Fig3C_PP_serotype_Mar27.pdf", width = 9, height = 11)

plot(
  map_pvnt2_plot,
  optimization_number = 1,
  xlim = c(-3,6),  #Adjust the xlim and ylim to plot the figure properly
  ylim = c(-8,3),  #Adjust the xlim and ylim to plot the figure properly
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

PP_ag_cords<-map_pvnt2_plot[["optimizations"]][[1]][["ag_base_coords"]]
rownames(PP_ag_cords)<-rownames(pvnt_df2)

PP_ag_cords<-as.data.frame(PP_ag_cords)%>%
  rename(X=V1,Y=V2)%>%
  mutate(Antigen=rownames(PP_ag_cords))
# Calculate the pairwise distances between Antigens based on their X and Y coordinates
n <- nrow(PP_ag_cords)
distances <- matrix(0, n, n)

for (i in 1:n) {
  for (j in 1:n) {
    distances[i, j] <- euclidean_distance(PP_ag_cords$X[i], PP_ag_cords$Y[i], PP_ag_cords$X[j], PP_ag_cords$Y[j])
  }
}
# Convert the matrix to a data frame and label rows and columns
PP_distances_df <- as.data.frame(distances)
colnames(PP_distances_df)<-rownames(PP_ag_cords)
rownames(PP_distances_df)<-rownames(PP_ag_cords)
write.xlsx(PP_distances_df,"2023/Fig4A_PP_AU.xlsx",overwrite = T)

##  For PPP samples
names_df_pvnt<-read.xlsx("2023/For antigenic cartograph_VAX_20230327.xlsx",sheet = "df")%>%
  select(Serum.panel,ID)%>%
  filter(Serum.panel=="PPP")
unique(names_df_pvnt$Serum.panel)

pvnt_df<-read.xlsx("2023/For antigenic cartograph_VAX_20230327.xlsx",sheet = "df")%>%
  filter(Serum.panel=="PPP")%>%
  select(-Serum.panel)%>%
  pivot_longer(-ID,names_to = "Virus",values_to = "NT50") %>% 
  mutate(NT50=round(NT50))%>%
  pivot_wider(names_from=ID, values_from=NT50) %>%
  mutate(Virus=gsub("\\.", " ", Virus))%>%
  mutate(Virus=gsub("BA ", "BA.", Virus))%>%
  mutate(Virus=case_when(Virus=="BA.4 6 1"~"BA.4.6.1",
                         Virus=="BA.2 75 2"~"BA.2.75.2",
                         Virus=="BQ 1 1"~"BQ.1.1",
                         T~Virus
  ))
pvnt_df2<-pvnt_df%>%select(-Virus)
rownames(pvnt_df2)<-pvnt_df$Virus
pvnt_df2<-as.matrix(pvnt_df2)

pvnt_df_test<-pvnt_df%>%pivot_longer(-Virus,names_to = "ID",values_to = "NT50") %>% 
  group_by(Virus)%>%
  summarise(sd=sd(log10(NT50)))

pvnt_df_test2<-pvnt_df%>%pivot_longer(-Virus,names_to = "ID",values_to = "NT50") %>% 
  group_by(ID)%>%
  summarise(sd=sd(log10(NT50)))


rownames(pvnt_df2)

map_pvnt <- acmap(
  titer_table = pvnt_df2,
  ag_names = rownames(pvnt_df2),
  sr_names = colnames(pvnt_df2),
)

map_pvnt <- optimizeMap(
  map                     = map_pvnt,
  number_of_dimensions    = 2,
  number_of_optimizations = 2000,
  minimum_column_basis    = "none",
)
view(map_pvnt)

map_pvnt

map_pvnt2<-relaxMap(map_pvnt)
view(map_pvnt2)

plotly_map_table_distance(map_pvnt2)
map_pvnt2 <- bootstrapMap(
  map                      = map_pvnt,
  method = "noisy",
  bootstrap_repeats        = 2000,
  optimizations_per_repeat = 400,
  ag_noise_sd              = 0.7,
  titer_noise_sd           = 0.7
)

boostrap_ag_coords_list <- mapBootstrap_agCoords(map_pvnt2)
boostrap_sr_coords_list <- mapBootstrap_srCoords(map_pvnt2)

agGroups(map_pvnt2)<-c(2,3,3,3,3,3,3,3,3,3,1)
srGroups(map_pvnt2)<-names_df_pvnt$Serum.panel
plotly_map_table_distance(map_pvnt2)

view(map_pvnt2)


# save.acmap(
#   map_pvnt2,
#   filename = "PPP_serotype_Mar27_1.ace",
#   compress = FALSE,
#   round_titers = FALSE
# )

map_pvnt2 <- read.acmap("2023/PPP_serotype_Mar27_1.ace")
map_pvnt2_plot<-map_pvnt2


test1<-c(rep("blue", 18))

srOutline(map_pvnt2_plot)<-test1
pdf(file = "2023/Fig3D_PPP_serotype_Mar27.pdf", width = 9, height = 11)

plot(
  map_pvnt2_plot,
  optimization_number = 1,
  xlim = c(-4,3),  #Adjust the xlim and ylim to plot the figure properly
  ylim = c(-5,4),  #Adjust the xlim and ylim to plot the figure properly
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
PPP_ag_cords<-map_pvnt2_plot[["optimizations"]][[1]][["ag_base_coords"]]
rownames(PPP_ag_cords)<-rownames(pvnt_df2)

PPP_ag_cords<-as.data.frame(PPP_ag_cords)%>%
  rename(X=V1,Y=V2)%>%
  mutate(Antigen=rownames(PPP_ag_cords))
# Calculate the pairwise distances between Antigens based on their X and Y coordinates
n <- nrow(PPP_ag_cords)
distances <- matrix(0, n, n)

for (i in 1:n) {
  for (j in 1:n) {
    distances[i, j] <- euclidean_distance(PPP_ag_cords$X[i], PPP_ag_cords$Y[i], PPP_ag_cords$X[j], PPP_ag_cords$Y[j])
  }
}
# Convert the matrix to a data frame and label rows and columns
PPP_distances_df <- as.data.frame(distances)
colnames(PPP_distances_df)<-rownames(PPP_ag_cords)
rownames(PPP_distances_df)<-rownames(PPP_ag_cords)
write.xlsx(PPP_distances_df,"2023/Fig4B_PPP_AU.xlsx")