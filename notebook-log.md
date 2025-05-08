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

16 sequences × 521 sites (aligned)

-------------------------------------------------------------------------------------

## Alignment ###

ClustalW and MAFFT will be compared. 

### ClustalW:
clustalw -ALIGN -INFILE=pfk-samples.fasta -OUTFILE=PFK-aligned_CW.fasta -OUTPUT=FASTA




### MAFFT:

install into conda environment:

conda activate alignment_env

conda install -c bioconda mafft


------------------------------------------------------------------------------
  MAFFT v7.526 (2024/Apr/26)
  https://mafft.cbrc.jp/alignment/software/
  MBE 30:772-780 (2013), NAR 30:3059-3066 (2002)
------------------------------------------------------------------------------

==> WARNING: A newer version of conda exists. <==
    current version: 25.1.1
    latest version: 25.3.1

"Please update conda by running"

    conda update -n base -c defaults conda

mafft --auto PFK-samples.fasta > PFK-aligned_MAFFT.fasta

input:PFK-samples.fasta
output: PFK-aligned_MAFFT.fasta



NEXT, see bottom:
iqtree -s PFK-aligned_MAFFT.fasta -m LG+F+G4 -bb 1000 -alrt 1000



------------------------------------------------------------------------------------

## Downloading RAxML and IQ Tree: 

From this link: https://github.com/amkozlov/raxml-ng
I downloaded the compiled binary for RAxML-NG, which is a phylogenetic tree inference tool which uses maximum-likelihood (ML) optimality criterion

For now, I placed this: /Users/sydneyjohnston/Desktop/Botany_563/Project/PPiPFK/RAxML_NG_1.2.2
which is where it must be run.

To check if it works, I entered the command: 

```
./raxml-ng -v
``` 

After going to finder to open the executable with a right click and accepting the risk of online download. 

To utilize the test data they provide: 
```
git clone https://github.com/amkozlov/ng-tutorial.git
``` 

Checking for MSA: 
```
./raxml-ng --check --msa ng-tutorial/bad.fa --model GTR+G
```

./raxml-ng --msa PFK-aligned_CW.fasta \
           --data-type AA \
           --model PROTGTR \
           --prefix T3 \
           --threads 2 \
           --seed 2



#### IQ Tree ####
In a new environment: 
  ```
  conda create -n phylogeny_env
  ```
  
I ran the install command for iqtree:
  ```
  conda install -c bioconda iqtree
  ```
  
### Running IQ Tree ### 
After checking iqtree -h (there are many modifiers), I ran a command on the ClustalW aligned file: 

```
iqtree -s PFK-aligned_CW.fasta
```

Analysis results written to: 
  IQ-TREE report:                PFK-aligned_CW.fasta.iqtree
  Maximum-likelihood tree:       PFK-aligned_CW.fasta.treefile
  Likelihood distances:          PFK-aligned_CW.fasta.mldist
  


### Installing MrBayes ###

Navigating to the desired directory (kdkdkdkdkd)
run the following set of commands: 

brew tap brewsci/bio
brew install mrbayes open-mpi
$which mb 

This will indicate the version. The output is: 

                            MrBayes 3.2.7a arm

                      (Bayesian Analysis of Phylogeny)

                             (Parallel version)
                         (1 processors available)

              Distributed under the GNU General Public License


               Type "help" or "help <command>" for information
                     on the commands that are available.

                   Type "about" for authorship and general
                       information about the program.


NEEDS TO BE NEXUS format INPUT
Fasta to Nexus file converter. 

seqmagick convert --output-format nexus --alphabet dna input.fasta output.nex
seqmagick convert --output-format nexus --alphabet prot input.fasta output.nex


_______________________________________________________________________________


I decided to run the alignment on MAFFT for comparison. See MAFFT. 

Then, to see that comparison better, I want to add bootstraps into the CW IQtree:
switch the environments: 

conda deactivate
conda activate phylogeny_env

iqtree -s PFK-aligned_CW.fasta -m LG+F+G4 -bb 1000 -alrt 1000

When attempting to rerun this IQtree for bootstraps, it is asking about overwriting the previously made files. 
I will have to do this, but in order to 'save' the previous results, I:
  1) updated git to my current version
  2) renamed the previous output files so they do not get overwritten. 
  
  PFK-aligned_CW.fasta.treefile  --> PFK-aligned_CW.firstrun.treefile
  PFK-aligned_CW.fasta.model.gz --> PFK-aligned_CW.firstrun.model.gz
  PFK-aligned_CW.fasta.mldist --> PFK-aligned_CW.firstrun.mldist
  PFK-aligned_CW.fasta.log  -->  PFK-aligned_CW.firstrun.log
  PFK-aligned_CW.fasta.iqtree  --> PFK-aligned_CW.firstrun.iqtree
  
  altered the command to:

iqtree -s PFK-aligned_CW.fasta -m LG+F+G4 -bb 1000 -alrt 1000 -redo

  Alignment has 16 sequences with 521 columns, 462 distinct patterns
  308 parsimony-informative, 106 singleton sites, 107 constant sites
                                Gap/Ambiguity  Composition  p-value
  Analyzing sequences: done in 1.90735e-05 secs using 47.19% CPU
     1  sp|P72830|PFKA1_SYNY3            30.71%    passed      9.78%
     2  sp|Q55988|PFKA2_SYNY3            26.30%    passed     22.93%
     3  tr|A0A081UD84|A0A081UD84_BACFG   35.51%    passed     98.52%
     4  sp|E0RTD9|PFKA_SPITD             29.75%    passed     37.52%
     5  sp|Q8NR14|PFKA_CORGL             33.59%    passed      6.49%
     6  tr|A0A0K6BYK5|A0A0K6BYK5_BACFG   37.43%    passed     99.93%
     7  sp|Q72H98|PFKA_THET2             38.20%    failed      0.13%
     8  tr|A0A0H2Z4Y4|A0A0H2Z4Y4_ECOK1   34.74%    passed     89.75%
     9  tr|A0A6M3ZEX6|A0A6M3ZEX6_BACSU   38.77%    passed     83.31%
    10  tr|I3VV24|I3VV24_THESW           38.39%    passed     68.83%
    11  tr|A3DEW4|A3DEW4_ACET2           37.81%    passed     93.56%
    12  tr|B5YEZ4|B5YEZ4_DICT6           38.77%    passed     11.38%
    13  tr|A3DCA6|A3DCA6_ACET2           20.35%    failed      3.87%
    14  tr|A0A2K2FEL3|A0A2K2FEL3_9CLOT   20.35%    failed      0.64%
    15  sp|O51669|PFKA_BORBU             14.20%    failed      0.02%
    16  sp|P06999|PFKB_ECOLI             40.69%    failed      0.56%
  ****  TOTAL                            32.22%  5 sequences failed composition chi2 test (p-value<5%; df=19)




Next, run that same command on the new MAFFT alignment. 

iqtree -s PFK-aligned_MAFFT.fasta -m LG+F+G4 -bb 1000 -alrt 1000

  Reading alignment file PFK-aligned_MAFFT.fasta ... Fasta format detected
  Reading fasta file: done in 0.000280857 secs using 40.59% CPU
  Alignment most likely contains protein sequences
  Alignment has 16 sequences with 608 columns, 499 distinct patterns
  294 parsimony-informative, 114 singleton sites, 200 constant sites
                                  Gap/Ambiguity  Composition  p-value
    Analyzing sequences: done in 0.00014019 secs using 7.846% CPU
       1  sp|Q8NR14|PFKA_CORGL             43.09%    passed      6.57%
       2  sp|P72830|PFKA1_SYNY3            40.62%    passed      9.80%
       3  sp|Q55988|PFKA2_SYNY3            36.84%    passed     23.14%
       4  tr|A0A081UD84|A0A081UD84_BACFG   44.74%    passed     98.52%
       5  sp|E0RTD9|PFKA_SPITD             39.80%    passed     37.67%
       6  sp|P06999|PFKB_ECOLI             49.18%    failed      0.56%
       7  sp|O51669|PFKA_BORBU             26.48%    failed      0.02%
       8  tr|B5YEZ4|B5YEZ4_DICT6           47.53%    passed     11.46%
       9  tr|A0A0H2Z4Y4|A0A0H2Z4Y4_ECOK1   44.08%    passed     89.84%
      10  tr|I3VV24|I3VV24_THESW           47.20%    passed     68.69%
      11  tr|A3DCA6|A3DCA6_ACET2           31.74%    failed      3.92%
      12  tr|A3DEW4|A3DEW4_ACET2           46.71%    passed     93.47%
      13  tr|A0A2K2FEL3|A0A2K2FEL3_9CLOT   31.74%    failed      0.65%
      14  tr|A0A6M3ZEX6|A0A6M3ZEX6_BACSU   47.53%    passed     83.17%
      15  tr|A0A0K6BYK5|A0A0K6BYK5_BACFG   46.38%    passed     99.92%
      16  sp|Q72H98|PFKA_THET2             47.04%    failed      0.13%
    ****  TOTAL                            41.92%  5 sequences failed composition chi2 test (p-value<5%; df=19)




Feature | ClustalW | MAFFT | Commentary/Comparison
Alignment Length | 521 columns | 608 columns | MAFFT aligned longer — added more gaps.
Distinct Patterns | 462 | 499 | MAFFT has slightly more unique patterns.
Informative Sites | 308 | 294 | ClustalW has slightly more parsimony-informative sites.
Singleton Sites | 106 | 114 | MAFFT has a few more "noise" sites.
Constant Sites | 107 | 200 | MAFFT has a lot more uniform positions.
Gap % | 32% | 42% | MAFFT introduced significantly more gaps.
Composition Failures | 5 | 5 | Same — no difference here.


Alignment | Strengths | Weaknesses
ClustalW :
    More parsimony-informative sites (strong signal for phylogeny). Lower gap %.
    Shorter alignment, maybe missing minor signal.
    
MAFFT :
    More total unique patterns. Longer alignment (could recover weak signals). 
    Higher gap %, slightly fewer informative sites, noisier singleton increase.



### Installing FigTree

To visualize the trees being created and compare the outputs of IQTree and MAFFT, I download FigTree. 
FigTree is a user-friendly program designed for visualizing, editing, and annotating phylogenetic trees,
making it ideal for quickly inspecting tree topology, branch support values, and evolutionary relationships

Fast and simple, local viewing. May use another program for the final product. 

https://github.com/rambaut/figtree/releases

FigTree v1.4.4

FigTree.v1.4.4.dmg
11.8 MB
Nov 25, 2018


Upon visualizing the tree with bootstraps, one of the E.coli samples consistently showed a bootstrap value of 0. 
Ultimately, this was because the sequence is for biologically distinct (pfkB isozyme 2), and is not homologous in the same functional class.

I am removing this sample and rerunning the IQtree command for better alignment 

conda activate phylogeny_env
iqtree -s PFK-aligned_MAFFT_noEcoli.fasta -m LG+F+G4 -bb 1000 -alrt 1000


new output:
  Reading alignment file PFK-aligned_MAFFT_noEcoli.fasta ... Fasta format detected
  Reading fasta file: done in 0.000324011 secs using 40.43% CPU
  Alignment most likely contains protein sequences
  WARNING: 36 sites contain only gaps or ambiguous characters.
  Alignment has 15 sequences with 608 columns, 452 distinct patterns
  290 parsimony-informative, 79 singleton sites, 239 constant sites
                                  Gap/Ambiguity  Composition  p-value
  Analyzing sequences: done in 0.000168085 secs using 6.544% CPU
       1  sp|Q8NR14|PFKA_CORGL             43.09%    passed      8.25%
       2  sp|P72830|PFKA1_SYNY3            40.62%    passed      7.74%
       3  sp|Q55988|PFKA2_SYNY3            36.84%    passed     20.68%
       4  tr|A0A081UD84|A0A081UD84_BACFG   44.74%    passed     98.69%
       5  sp|E0RTD9|PFKA_SPITD             39.80%    passed     38.37%
       6  sp|O51669|PFKA_BORBU             26.48%    failed      0.03%
       7  tr|B5YEZ4|B5YEZ4_DICT6           47.53%    passed     13.24%
       8  tr|A0A0H2Z4Y4|A0A0H2Z4Y4_ECOK1   44.08%    passed     93.15%
       9  tr|I3VV24|I3VV24_THESW           47.20%    passed     74.41%
      10  tr|A3DCA6|A3DCA6_ACET2           31.74%    failed      4.06%
      11  tr|A3DEW4|A3DEW4_ACET2           46.71%    passed     95.35%
      12  tr|A0A2K2FEL3|A0A2K2FEL3_9CLOT   31.74%    failed      0.77%
      13  tr|A0A6M3ZEX6|A0A6M3ZEX6_BACSU   47.53%    passed     87.42%
      14  tr|A0A0K6BYK5|A0A0K6BYK5_BACFG   46.38%    passed     99.94%
      15  sp|Q72H98|PFKA_THET2             47.04%    failed      0.07%
    ****  TOTAL                            41.44%  4 sequences failed composition chi2 test (p-value<5%; df=19)

ClustalW:
iqtree -s PFK-aligned_CW_noEcoli.fasta -m LG+F+G4 -bb 1000 -alrt 1000

  
  Reading alignment file PFK-aligned_CW_noEcoli.fasta ... Fasta format detected
  Reading fasta file: done in 0.000121832 secs using 73.87% CPU
  Alignment most likely contains protein sequences
  Alignment has 15 sequences with 521 columns, 451 distinct patterns
  305 parsimony-informative, 89 singleton sites, 127 constant sites
                                  Gap/Ambiguity  Composition  p-value
  Analyzing sequences: done in 5.31673e-05 secs using 16.93% CPU
       1  sp|P72830|PFKA1_SYNY3            30.71%    passed      7.72%
       2  sp|Q55988|PFKA2_SYNY3            26.30%    passed     20.50%
       3  tr|A0A081UD84|A0A081UD84_BACFG   35.51%    passed     98.69%
       4  sp|E0RTD9|PFKA_SPITD             29.75%    passed     38.23%
       5  sp|Q8NR14|PFKA_CORGL             33.59%    passed      8.17%
       6  tr|A0A0K6BYK5|A0A0K6BYK5_BACFG   37.43%    passed     99.94%
       7  sp|Q72H98|PFKA_THET2             38.20%    failed      0.07%
       8  tr|A0A0H2Z4Y4|A0A0H2Z4Y4_ECOK1   34.74%    passed     93.09%
       9  tr|A0A6M3ZEX6|A0A6M3ZEX6_BACSU   38.77%    passed     87.54%
      10  tr|I3VV24|I3VV24_THESW           38.39%    passed     74.53%
      11  tr|A3DEW4|A3DEW4_ACET2           37.81%    passed     95.41%
      12  tr|B5YEZ4|B5YEZ4_DICT6           38.77%    passed     13.16%
      13  tr|A3DCA6|A3DCA6_ACET2           20.35%    failed      4.02%
      14  tr|A0A2K2FEL3|A0A2K2FEL3_9CLOT   20.35%    failed      0.76%
      15  sp|O51669|PFKA_BORBU             14.20%    failed      0.02%
    ****  TOTAL                            31.66%  4 sequences failed composition chi2 test (p-value<5%; df=19)



