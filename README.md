# PPiPFK
Analyzing the distribution of pyrophosphate phosphofructokinases
This is the repository for Phylogenetics Analysis Class 563

Date: January 23, 2025 

enetic Analysis of Bacterial Phosphofructokinases (PFKs)

This repository contains the reproducible workflow for constructing a maximum likelihood phylogenetic tree of bacterial ATP- and PPi-dependent phosphofructokinases (PFKs). 
The goal is to reassess evolutionary relationships among PFK homologs using curated protein sequences referenced in Kuil et al. (2023).

---

## Project Structure

```
project/
│
├── data/
│   ├── PFK-samples.fasta               # Input FASTA of all 16 curated sequences
│   ├── PFK-aligned_CW.fasta            # ClustalW alignment output
│   ├── PFK-aligned_MAFFT.fasta         # MAFFT alignment output
│   ├── *_noEcoli.fasta                 # Same alignments with pfkB isozyme removed
│
├── trees/
│   ├── *.treefile                      # ML trees output by IQ-TREE
│   ├── *.iqtree                        # Summary reports
│   └── *.log                           # Run logs
│
├── scripts/
│   └── phylo_workflow.sh               # Reproducible workflow shell script
│
└── README.md
```

---

## Sequence Sources

All sequences were downloaded from UniProt based on accession IDs reported in Kuil et al. (2023). 
The initial dataset included 16 sequences (ATP- and PPi-dependent PFKs); the divergent pfkB isozyme from *E. coli* (P06999) was later removed.

---

## Dependencies

Install via conda:
```bash
conda create -n phylo_env -c bioconda iqtree mafft clustalw seqmagick
conda activate phylo_env
```

---

## Tree Inference Summary

- Alignments: ClustalW and MAFFT
- Tree Inference: IQ-TREE v2.3.2
- Model: LG+F+G4 (empirical matrix, observed frequencies, gamma-distributed rates)
- Bootstraps: 1000 UFBoot + 1000 SH-aLRT
- Visualization: FigTree / iTOL

---

## Final Tree

Final tree was constructed from ClustalW alignment (pfkB removed) using IQ-TREE with the LG+F+G4 model and 1000 bootstrap replicates. 
Visualization and clade annotation were performed using iTOL and Illustrator.

---

## Citation

If using or adapting this workflow, please reference the following for additional information:
- Nguyen et al., 2015 (IQ-TREE)
- Katoh & Standley, 2013 (MAFFT)
- Thompson et al., 1994 (ClustalW)
- Kozlov et al., 2019 (RAxML-NG)
