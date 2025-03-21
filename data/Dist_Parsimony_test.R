######################################################
###### Library Prep and Reading the Data #############
######################################################


#INSTALL necessary packages
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)

#LOAD the necessary packages
library(ape)
library(ade4)
library(stats)
library(adegenet)
library(phangorn)

#LOAD the sample data
protein_data <- read.dna("/Users/sydneyjohnston/Desktop/Botany_563/Project/PPiPFK/data/PFK-aligned_CW.fasta", format = "fasta", as.character = TRUE)

# Check the structure of the data to ensure it's read correctly
str(protein_data)    # Show details of the object

# Convert the protein sequences into a matrix
protein_matrix <- do.call(rbind, lapply(protein_data, function(x) unlist(strsplit(x, ""))))

# View the first 5 sequences and first 10 positions
protein_matrix[1:5, 1:10]

#View the entire dataset in spreadsheet:
View(protein_matrix)

--------

#####################################################
####### DISTANCE BASED PHYLOGENY ###########
#####################################################

#COMPUTE genetic differences
#They choose a Tamura and Nei 1993 model which allows for different rates of transitions
#and transversions, heterogeneous base frequencies, and between-site variation of the 
#substitution rate (more on Models of Evolution).
D <- dist.dna(dna, model="TN93")
class(D)
length(D)

temp <- as.data.frame(as.matrix(D))
table.paint(temp, cleg=0, clabel.row = 5, clabel.col = 5)

temp <- t(as.matrix(D))
temp <- temp[,ncol(temp):1]

par(mar=c(1,5,5,1))
image(x=1:80, y=1:80, temp, col=rev(terrain.colors(100)), xaxt="n", yaxt="n", xlab = "", ylab = "")
      axis(side=2, at=1:80, lab=rownames(dna), las=2, cex.axis=.5)
      axis(side=3, at=1:80, lab=row.names(dna), las=3, cex.axis=.5)

#GET the NJ Tree
tre <- nj(D)

#Before plotting, we can use the ```ladderize``` function which reorganizes the internal 
#structure of the tree to get the ladderized effect when plotted
tre <- ladderize(tre)


#PLOT the tree
plot(tre, cex=.6)
title("A simple NJ tree")

plot(tre, show.tip=FALSE)
title("Unrooted NJ Tree")
myPal <- colorRampPalette(c("red", "yellow", "green", "blue"))
tiplabels(annot$year, bg=num2col(annot$year, col.pal = myPal), cex=.5)
temp <- pretty(1993:2008, 5)
legend("bottomleft", fill=num2col(temp,col.pal = myPal), leg=temp, ncol = 2)


plot(tre, type="unrooted", show.tip=FALSE)
title("Unrooted NJ tree")

head(annot)

tre2 <- root(tre, out=1)
tre2 <- ladderize(tre2)

plot(tre2, show.tip=FALSE, edge.width=2)
title("Rooted NJ tree")
tiplabels(tre$tip.label, bg=transp(num2col(annot$year, col.pal=myPal),.7), cex=.5,
           fg="transparent")

axisPhylo()
temp <- pretty(1993:2008, 5)
legend("topright", fill=transp(num2col(temp, col.pal=myPal),.7), leg=temp, ncol=2)




