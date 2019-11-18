# SNP-to-alignments
R code for deriving alignments from GBS data

This code enables the derivation of short sequence alignments from genotyping-by-sequencing data.
Required inputs are (1) a filtered .vcf file containing variant sites and (2) binary alignment mapping (BAM) files.
Included are notes explaining each step, and three filters to remove multi-copy loci from the resulting alignments, 
described in detail in Shaik et al. (2020)
