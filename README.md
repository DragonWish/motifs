motifs
======

Motif Generator: FIND_MOTIF
--------------

Finds a sequence motif from a set of provided peptide sequences.
Uses the Gibbs Motif Finding algorithm as outlined in the Bioinformatics Algorithms textbook:
https://stepic.org/Bioinformatics-Algorithms-2/


usage: `FIND_MOTIF.py inputfile outputfile k motif_file`

**inputfile**: collection of unaligned sequences belongign to the same protein family

**outputfile**: filename where the motif sequences from each original sequence (from inputfile) will be outputted. This file can then be uploaded to servers such as [Weblogo 3](http://weblogo.threeplusone.com/create.cgi) to create motif logos. 

**k**: length of the motif

**motif_file**: filename where the motif profile matrix will be outputted



Motif Finder from a proteome: find_most_probable.py
-----------

From a specified proteome or list of proteomes, rank all proteins in terms of the probability of the sequence best matching a given motif profile being a genuine motif. The boxplot function generates boxplots of the proteomes, with whiskers at 3IQR and outliers indicated. 
Linees 157-159 (commented out) also sorts and lists the top 25 proteins according to probability. 

Usage: `find_most_probable.py input proteomes output`

**input**: input motif file. Usually motif_file as outputted by FIND_MOTIF.py

**proteomes**: a list of the proteomes to be assessed. This should be a text file with the filenames of the fasta formatted proteins on each line. The boxplot uses these filenames for labeling the results.

**output**: output image file where the boxplot will be saved. 

Sample Files
-----------

**athaliana.fasta**: Sample proteome of the plant *Arabidopsis thaliana*

**spombe.fasta**: Sample proteome of the fungus *Schizosaccharomyces pombe*

**list_proteomes.txt**: Sample input file for `find_most_probable.py` listing the proteomes to be queried

**hpt_sequences.txt**: sample input reference set for `FIND_MOTIF.py`
