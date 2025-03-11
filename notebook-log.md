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

The UniProt KEGG Orthology numbers I have stored are: 
- K00850
- K16370
- K21071
- K24182

These are all links to the amino acid sequences

----------------------------------------------------------------------------

For the larger project, I will be utilizing a mass of NCBI representative genomes to get a broad idea of the distribution of PPi-PFK and similar homologs across genus. I want to compare organisms that harbor ATP- and PPi- dependent genes, and in some cases organisms have the genes for both versions of the protein. 

Details for the NCBI Reference Genome Download: 
- Eukaryote search filters: 
  Reference Genomes (Representative Genomes)
  Annotated Genomes
  Complete Assembly 
  All years
  Total: 121 Genomes
  8 GB 
  
- Archaea search filters: 
  Reference Genomes (Representative Genomes)
  Annotated Genomes
  Complete Assembly 
  All years
  Manual: Propionibacterium
  Total: 363 Genomes
  1 GB
  
- Bacteria search filters: 
  Reference Genomes (Representative Genomes)
  Annotated Genomes
  Complete Assembly 
  All years
  Total: 5,761 Genomes
  10 GB
  
  
  
Column selection: 
Assembly, RefSeq, GenBank, TaxID, Modifier, Size (Mb), Scientific name, Annotation, Level, Release date, WGS accession, Scaffolds count, Sequencing technology, Genes, Protein-coding, CheckM contamination (%), CheckM completedness (%), contig N50

Downloading the Package (Feb 18th, 2025): ALL file source, check boxes for genome sequences (FASTA), Protein (FASTA), and for euk and arch -- include: Sequence and Annotation (GBFF). 

  
Include: 
  Trichomonas vaginalis
  Thermoproteus
  Pyrococcus
  Thermococcales
  Propionibacterium
  

----------------------------------------------------------------------------

To reduce the number of organisms I analyze for the class project I will only choose the ~15 organisms specifically discussed in the paper listed in #### Data Collection ####. 

The Uniprot for those organisms were used to download protein sequences: 
Q8NR14 - 346 aa       ATP
P72830 - 361 aa       ATP
Q55988 - 384 aa       ATP
A0A081UD84 - 336 aa   ATP
E0RTD9 - 366 aa       ATP
P06999 - 309 aa       ATP
O51669 - 447 aa       ATP
B5YEZ4 - 319 aa       ATP
A0A0H2Z4Y4 - 340 aa   ATP
I3VV24 - 321 aa       ATP
A3DCA6 - 415 aa       PPi
A3DEW4 - 324 aa       ATP
A0A2K2FEL3 - 415 aa   PPi
A0A6M3ZEX6 - 319 aa   ATP
A0A0K6BYK5 - 326 aa   ATP
Q72H98 - 322 aa       ATP


-------------------------------------------------------------------------------------

## Alignment ###

ClustalW and MAFFT will be compared. 

### ClustalW:
clustalw -ALIGN -INFILE=pfk-samples.fasta -OUTFILE=PFK-aligned_CW.fasta -OUTPUT=FASTA




IUSEFOWEJNFOI VPOKME FCVIHJNSDS FEWPIOJ 

********** ADD INFO **************
LIUHS NDFCLIUESDNF CWESIDUHFCVJ)MEW FCIN


------------------------------------------------------------------------------------

## Downloading RAxML and IQ Tree: 

From this link: https://github.com/amkozlov/raxml-ng
I downloaded the compiled binary for RAxML-NG, which is a phylogenetic tree inference tool which uses maximum-likelihood (ML) optimality criterion

For now, I placed this: /Users/sydneyjohnston/Desktop/Botany_563/Project/PPiPFK/RAxML_NG_1.2.2
which is where it must be run.

To check if it works, I entered the command: 

```./raxml-ng -v``` 

After going to finder to open the executable with a right click and accepting the risk of online download. 




