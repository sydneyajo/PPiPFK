###### Library Prep and Reading the Data #############
######################################################

# INSTALL necessary packages (if not already installed)
install.packages("ape", dependencies = TRUE)
install.packages("phangorn", dependencies = TRUE)
install.packages("seqinr", dependencies = TRUE)

# LOAD the necessary packages
library(ape)
library(phangorn)
library(seqinr)

# LOAD the sample data
protein_data <- read.alignment("/Users/sydneyjohnston/Desktop/Botany_563/Project/PPiPFK/data/PFK-aligned_CW.fasta", format = "fasta")

# Check the structure of the data to ensure it's read correctly
str(protein_data)

######################################################
####### DISTANCE BASED PHYLOGENY #####################
######################################################

# Convert the protein data into a phyDat object for phangorn
protein_phyDat <- as.phyDat(protein_data)

# Calculate pairwise distance matrix using dist.hamming (for proteins)
D <- dist.hamming(protein_phyDat)

# Check the class and length of the distance object
class(D)
length(D)

# Convert the distance object to a matrix
D_matrix <- as.matrix(D)

# Plot the distance matrix
par(mar = c(4, 4, 2, 2))  # Adjust margins
image(x = 1:nrow(D_matrix), y = 1:ncol(D_matrix), D_matrix, col = rev(terrain.colors(100)), 
      xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(side = 2, at = 1:nrow(D_matrix), labels = rownames(D_matrix), las = 2, cex.axis = 0.6)
axis(side = 3, at = 1:ncol(D_matrix), labels = rownames(D_matrix), las = 3, cex.axis = 0.6)

# GET the NJ Tree
tre <- nj(D)

# Ladderize the tree for better visualization
tre <- ladderize(tre)

# PLOT the tree
plot(tre, cex = 0.6)
title("A simple NJ tree")

# Plot an unrooted NJ tree
plot(tre, show.tip = FALSE)
title("Unrooted NJ Tree")

# Root the tree and ladderize it
tre2 <- root(tre, out = 13)
tre2 <- ladderize(tre2)

-----------------------
  # Check for outgroup fittings
  
  # Root the tree with sequence 14 as the outgroup
  tre_outgroup_14 <- root(tre, out = 14)

# Ladderize the rooted tree
tre_outgroup_14_ladderized <- ladderize(tre_outgroup_14)

# Modify the tip labels to display only the last 5 characters
tre_outgroup_14_ladderized$tip.label <- substr(tre_outgroup_14_ladderized$tip.label, 
                                               nchar(tre_outgroup_14_ladderized$tip.label) - 4, 
                                               nchar(tre_outgroup_14_ladderized$tip.label))

# Plot the tree with reduced labels
plot(tre_outgroup_14_ladderized, show.tip = TRUE, edge.width = 2, cex = 0.6)
title("Ladderized Rooted NJ Tree with Sequence 14 as Outgroup")

# Add rounded branch lengths as labels (maximum 5 digits)
#rounded_edge_lengths <- signif(tre_outgroup_14_ladderized$edge.length, digits = 3)
#edgelabels(rounded_edge_lengths, cex = 0.5)
