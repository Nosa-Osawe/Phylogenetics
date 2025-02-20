install.packages("remotes")  # Needed to install V.PhyloMaker
remotes::install_github("jinyizju/V.PhyloMaker")
install.packages("ape")  # For visualizing the tree

library(pheatmap)
library(V.PhyloMaker)
library(ape)

taxa_list <- data.frame(
  species = c("Quercus_robur", "Quercus_alba", "Pinus_sylvestris", "Betula_pendula"),
  genus = c("Quercus", "Quercus", "Pinus", "Betula"),
  family = c("Fagaceae", "Fagaceae", "Pinaceae", "Betulaceae")
)

tree <- phylo.maker(taxa_list, output.tree = TRUE)

tree_calibrated <- phylo.maker(taxa_list, output.tree = TRUE,
                               output.sp.list = TRUE, scenario = "S3")
plot(tree$scenario.3, show.tip.label = TRUE)

plot(tree_calibrated$scenario.3, show.tip.label = TRUE)

# Compute the phylogenetic distance matrix
phylo_dist <- cophenetic(tree$scenario.3)

# View the distance matrix
print(phylo_dist)

# Example: Distance between Quercus_robur and Pinus_sylvestris
phylo_dist["Quercus_robur", "Pinus_sylvestris"]


# Create a heatmap of phylogenetic distances
pheatmap(phylo_dist, 
         color = colorRampPalette(c("white", "blue"))(50),
         main = "Phylogenetic Distance Matrix")


library(vegan)

# Perform Principal Coordinates Analysis (PCoA)
pcoa_res <- cmdscale(phylo_dist, k = 2)

# Plot PCoA
plot(pcoa_res, xlab = "PCoA 1", ylab = "PCoA 2", 
     main = "Phylogenetic PCoA", pch = 19, col = "darkgreen")
text(pcoa_res, labels = rownames(pcoa_res), pos = 3)


library(phytools)

# Example: Simulated trait data (e.g., leaf size, height)
trait_values <- c(10.2, 10.4, 4.0, 13.1)  
names(trait_values) <- taxa_list$species  # Assign names to match species

# Calculate Blomberg’s K for phylogenetic signal
phylo_signal <- phylosig(tree$scenario.3, trait_values, method = "K")
print(phylo_signal)
