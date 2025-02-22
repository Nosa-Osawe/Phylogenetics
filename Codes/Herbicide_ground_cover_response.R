library(pheatmap)
library(V.PhyloMaker)
library(ape)
library(tidyverse)
library(readxl)
library(vegan)
library(phytools)

plant <- read_excel("C:\\Users\\DELL\\Documents\\Git in R\\Phylogenetics\\Data\\Herbicide_Phylogenetic.xlsx", 
                    sheet = "Sheet1")
view(plant)
sum(is.na(plant$Q11))

plant$Q11 <- as.numeric(plant$Q11)

plant.list <- plant %>% 
  filter(Week==0, Belt=="B1") %>% 
  select(Scientific_name, Genus,Family ) %>% 
  rename(species = Scientific_name,
         genus = Genus, 
         family = Family) %>% 
  as.data.frame()


plant_tree <-  phylo.maker(plant.list, output.tree = TRUE, output.sp.list = TRUE)

plot(plant_tree$scenario.3, show.tip.label = TRUE)

# Compute the phylogenetic distance matrix
distance <- cophenetic(plant_tree$scenario.3)

# heatmap of phylogenetic distances
pheatmap(distance, 
         color = colorRampPalette(c("white", "red"))(50),
         main = "Phylogenetic Distance Matrix")

# Pre-spray ground cover, collapsed (averaged) by row [Q]--making 17 surveys

pre_spay <- plant %>% 
  filter(Spray=="Pre-spray") %>% 
  pivot_longer(
    cols = -c(1:9),
    names_to = "Q",
    values_to = "Rows"
  ) %>% 
  as.data.frame() %>%
  group_by(Q, Scientific_name) %>% 
  summarise(mean_ground = mean(Rows)) %>% 
  as.data.frame()

q1_cover <- pre_spay %>% 
  filter(Q=="Q1") %>% 
  select(mean_ground) %>% 
  as.data.frame()

species <-  pre_spay %>% 
  filter(Q=="Q1") %>% 
  select(Scientific_name) %>% 
  as.data.frame()
species<-species$Scientific_name
q1_cover<- as.vector(t(q1_cover))
names(q1_cover) <- species 

# Blomberg’s K for phylogenetic signal
q1.signal <- phylosig(plant_tree$scenario.3, q1_cover, method = "K")
print(q1.signal)





