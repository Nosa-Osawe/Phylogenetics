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


