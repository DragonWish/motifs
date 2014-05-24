motifs
======

Motif Generator: FIND_MOTIF

Finds a sequence motif from a set of provided peptide sequences.
Uses the Gibbs Motif Finding algorithm as outlined in the Bioinformatics Algorithms textbook:
https://stepic.org/Bioinformatics-Algorithms-2/


Usage: FIND_MOTIF.py inputfile.txt outputfile.txt k



Motif Finder from a proteome: find_most_probable.py

From a specified proteome (line 162), rank all proteins in terms of the probability of the sequence best matching a given motif profile being a genuine motif. The boxplot function generates boxplots of the proteomes, with whiskers at 3IQR and outliers indicated. 
Linees 170-172 (commented out) also sorts and lists the top 25 proteins according to probability. 
