# _Mycobacteria_ MAB_0812 & MAB_0813c Genomic Analysis Files

## Introduction
This GitHub repository houses the files for the genomic analysis of the MAB_0812 and MAB_0813c genes in species across Mycobacteria. Analysis conducted for the de Pas lab at the University of Pittsburgh. 

MAB_0812 is often referred to as 03513, and MAB_0813c as 03514. 

Main files of interest are the figures in the Length Histograms, Logo Sequences, and Tree Files folders. Other folders/files contributed to the analysis that created these figures. 

## Directory Structure
- **03513**: houses the amino acid sequences of the Mycobacteria BLAST results for MAB_0812 (or 03513). Both raw sequences and aligned sequences, in FASTA format. 

- **03514**: houses the amino acid sequences of the Mycobacteria BLAST results for MAB_0813c (or 03514). Both raw sequences and aligned sequences, in FASTA format. 

- **Length Histograms**: houses the histograms showing the length distribution of 03513 and 03514 orthologs returned from BLAST. Lengths of mutated evolved 03513/03514 and wildtype 253a 03513/03514 gene sequences are marked on histograms.

- **Logo Sequences**: houses the logo sequences created from the multiple sequence alignments of 03513 and 03514, showing conserved vs. variable residues. 

- **Spreadsheet Analyses**: houses various spreadsheet analyses conducted for this project: histograms, sequence lengths, common vs. unique species in the BLAST results, M. smegmatis orthologs, and the information underlying the tree. M. abscessus spreadsheet pending.

- **Tree Files**: houses the tree file created from representative orthologs of 03513 (one ortholog per species, for the species that had >=3 orthologs in the BLAST results for 03513). Synteny of 03514 and isolation source information are annotated on the tree. 

- **Blast Sequence Length Analysis**: Jupyter Notebook file analyzing the lengths of the BLAST results. Creates the histograms in the Length Histograms folder.

- **fetchGenomesFromProteinAccessions**: Python script that pulls the genomes that encode proteins on NCBI. 

## Contact Information
Isabelle D'Amico
idamico@andrew.cmu.edu
