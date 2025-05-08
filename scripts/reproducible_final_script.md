#!/bin/bash

# Reproducible workflow for PFK phylogeny
# Author: Sydney Bradley
# Environment: phylo_env (conda)

# Acquired a dataset from Kuil 2023 paper:
# "Pyrophosphate as allosteric regulator of ATP-phosphofructokinase in Clostridium thermocellum and other 
# bacteria with ATP- and PPi-phosphofructokinases" by Teun Kuil, Carolus  M.K. Nurminen , Antonius J.A. van Maris: 
# https://pubmed.ncbi.nlm.nih.gov/37380119/
# Downloaded genes from UniProt

# Step 1: Align sequences with ClustalW
clustalw -ALIGN -INFILE=data/PFK-samples.fasta -OUTFILE=data/PFK-aligned_CW.fasta -OUTPUT=FASTA

# Step 2: Align sequences with MAFFT (optional)
mafft --auto data/PFK-samples.fasta > data/PFK-aligned_MAFFT.fasta

# Run IQTree to find a suited model

iqtree -s PFK-aligned_CW.fasta


# Run RAxML as a comparison (optional)
./raxml-ng --msa PFK-aligned_CW.fasta \
           --data-type AA \
           --model PROTGTR \
           --prefix T3 \
           --threads 2 \
           --seed 2

#Assess the models and determine if the iqtree model was a better fit (based on near 0 total). 

# Roughly visualize and compare the IQTree outputs using FigTree

# Step 3: Remove problematic pfkB isozyme from E. coli (manual editing)
# Note: Save new edited alignment files as *_noEcoli.fasta

# Step 4: Re-Run IQ-TREE on ClustalW alignment
iqtree -s data/PFK-aligned_CW_noEcoli.fasta -m LG+F+G4 -bb 1000 -alrt 1000 -nt AUTO

# Step 5: Re-Run IQ-TREE on MAFFT alignment (for comparison)
iqtree -s data/PFK-aligned_MAFFT_noEcoli.fasta -m LG+F+G4 -bb 1000 -alrt 1000 -nt AUTO


# Step 6: Visualize tree
# Use FigTree locally to confirm then upload .treefile to iTOL (https://itol.embl.de)
# Reconfigure the tree to root from outgroup and be visually interpretable
# Export to Illustrator to adjust Tip Labels 


# For further detail, refer to notebook.log and/or README file as necessary. These contain download/install information, troubleshooting processes, and additional method attempts. 
