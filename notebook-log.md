This project will focus on analyzing the diversity of the PPi-PFK protein. 

### TITLE:

Comparative Genomics Uncovers the Evolutionary Trajectory and Broader Distribution of Pyrophosphate-Dependent Phosphofructokinase (PPi-PFK) as an Alternative Phosphorylation Substrate in Glycolysis

### DESCRIPTION:

This study employs comparative genomics and phylogenetics to examine the distribution, evolutionary dynamics, and functional role of pyrophosphate-dependent phosphofructokinase (PPi-PFK) across diverse organisms. By identifying homologs and previously unannotated sequences, we aim to expand recognition of PPi-PFK beyond known species. Our findings will elucidate its evolutionary trajectory, metabolic significance, and relationship to ATP-dependent phosphofructokinases in central carbon metabolism.

### DATA COLLECTION: 

Approximately 3500 sequences were collected from 2023 "Pyrophosphate as allosteric regulator of ATP-phosphofructokinase in Clostridium thermocellum and other bacteria with ATP- and PPi-phosphofructokinases" by Teun Kuil, Carolus  M.K. Nurminen , Antonius J.A. van Maris : https://pubmed.ncbi.nlm.nih.gov/37380119/

- Some of these sequences include ATP- and PPi-PFKs (with and without PPase), and ATP-only PFK. 
  
- I have also run BLASTp on my target organism, *C. thermocellum* using: 

### BASE SCRIPT

```
#!/bin/bash
source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh
unset $PYTHONPATH

#Activate the example BLAST environment from the Readme file (you need to follow the commands to create it first)
conda activate blast.env

blastp -query Ctherm_DSM1313.fasta -db nr -out CthermPPi_blastp_output.tsv -outfmt '6 qstart qend sstart send sacc ssciname stitle pident evalue qlen sseq' -evalue 0.00000001 -max_target_seqs 10000000 -remote

conda deactivate blast.env
```

#### Details as of Feburary 11th, 2025

I ran this script on C. thermocellum PPi-PFK (positive control/target sequence), *C. thermocellum* ATP-PFK (negative control), and Human Liver ATP-PFK (Eukaryotic sample). 

As of Feb 11th, 2025, I am taking note of more organisms that are referred to in literature to utilize PPi-PFK, and creating a spreadsheet for a pHMMR run to collect more homologs. 




