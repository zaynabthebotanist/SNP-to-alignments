# SNP-to-alignments
R code for deriving sequence alignments from GBS data

This code enables the derivation of short sequence alignments from genotyping-by-sequencing data.
Required inputs are (1) a filtered .vcf file containing variant sites and (2) binary alignment mapping (BAM) files.
Included in the code are notes explaining each step, and three filters applicable to diploid taxa that remove multi-copy loci from the resulting alignments, as detailed in Shaik et al. (2020).

Note that Stacks v2.0 (Rochette et al., 2018) offers the same de novo haplotype assembly capability. Initial comparison of the two pipelines indicates that the R-based workflow presented here delivers more parsimony-informative site-rich sequences than Stacks

Further note - do NOT use the vcf2nexus converter offered by shankarkshakya/mypackage in R. The coding was in 0, 1 and ? rather than 0, 1, and 2 as it should be. I would recommend using something like PGDSpider for this conversion instead
