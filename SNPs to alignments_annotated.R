
# This script, used in Shaik et al. (2020), generates sequence alignments from GBS data 
# two inputs are required:

# (1) a filtered .vcf file including scaffold identifiers
# (2) BAM files for each individual containing reads

library(adegenet)
library(ape)
library(Biostrings)
library(DECIPHER)
library(GenomicAlignments)
library(haplotypes)
library(LEA)
library(msa)
library(msaR)
library(muscle)
library(pegas)
library(plyr)
library(rlist)
library(Rsamtools)
library(Rsamtools)
library(seqinr)
library(ssviz)
library(vcfR)


# read vcf
vcf <- read.vcfR("")

# include only relevant samples (by name) in c()
# RE that the gt element in a vcf includes a column called FORMAT, which we also want to keep
# so start like: c("FORMAT", "sample1", "sample2" etc...)

vcf <- vcf[, which(colnames(vcf@gt) %in% c())]

# exclude loci with any missing data. The rationale is that
# if we're going to use the variant scaffold identifiers to produce alignments, we want to 
# use those variants that have been captured across all individuals to reduce missingness 

queryMETA(vcf, "DP")
dp <- extract.gt(vcf, as.numeric = TRUE)
sum(is.na(dp[,1]))
myMiss <- apply(dp, MARGIN = 1, function(x){ sum( is.na(x) ) } )
myMiss <- myMiss / ncol(dp)
vcf <- vcf[myMiss < 0.05, ]
vcf

#query the scaffolds identifiers with reads containing variants mapped to them and assign to variable "tag"
head(vcf)
tag <- vcf@fix[ ,1] 
tag <- data.frame(tag)
names(tag) <- "name"
write.csv(tag, "tag.csv")

# next, let's import the binary allignment map (BAM) files using Rsamtools

# first we want to filter the BAM files to include only sequences containing variants
# These files are typically large and require lots of memory
# as such, individual import, filtering and export as .RData files is required
# the following three loops are for one BAM file

tag <- read.csv("tag.csv")[,2]

filenames <- as.list(list.files(pattern = "*.bam"))[1]  
n <- c(1:length(filenames))
bam <- list()
for (i in n) {
  bam[[i]] <- ssviz::readBam(filenames[[i]])
}

bamint <- list()
for (j in n) {
  bamint[[i]] <- bam[[i]][bam[[i]]$rname %in% unique(tag), ]
}

save(bamint, file = "bamtest.RData") 

#generate consensus sequences from short reads + align them

#load BAMs and set
filenames <- as.list(paste0("bamtest", ".RData"))  #this can be modified to load multiple files
load(filenames[[1]])

scaff <- bamint[[1]]
scaff$rname <- as.character(scaff$rname)
split <- split(scaff, f = scaff$rname) #splits read data by scaffold aligned to
split[[1]] # this is the first element in the list "split", note that the reads are all
           # aligned to the same scaffold, as we wanted

#####################
########LOOPS########
#####################

# at this point, we have a series of reads aligned to the same scaffold
# now let's haplotype them

## (1) create DNAStringSet objects
n <- c(1:length(split))
dna <- list()
for (i in n){
  dna[[i]] <- Biostrings::DNAStringSet(c(split[[i]]$seq))
}

n <- c(1:length(dna))
for (i in n) {
  names(dna[[i]]) <- paste0(split[[i]]$rname)
}

## (2) align sequences using muscle
n <- c(1:length(dna))
aln <- list()
for (j in n){
  aln[[j]] <- muscle::muscle(dna[[j]])
}

## (3) Convert from DNAMultipleAlignment to DNAbin
n <- c(1:length(aln))
string <- list()
for (l in n){
  string[[l]] <- as.DNAbin(DNAStringSet(aln[[l]]))
}

## (4) haplotyping
##but first exclude reads with depth 1 (cannot be haplotyped):
string <- string[lengths(string) > 1]  ##Cannot haplotype one read, so exclude these

n <- c(1:length(string))   ##Make haplotypes
hap <- list()
for (i in n) {
  hap[[i]] <- pegas::haplotype(string[[i]])   ##needs this format to sort:
}

##Now to sort by the most to the least common haplotypes:
##the pegas::sort function is acting up, run source code
##from (https://rdrr.io/cran/pegas/src/R/haplotype.R):

oc <- oldClass(hap[[1]])
from <- "string[[i]]"
what <- c("frequencies")

idx <- list()
n <- c(1:length(hap))
for (i in n) {
  idx[[i]] <- attr(hap[[i]], "index")
}

o <- list()
n <- c(1:length(hap))
for (i in n) {
  o[[i]] <- switch(what,
              frequencies = order(sapply(idx[[i]], length), decreasing = TRUE),
              labels = order(rownames(hap[[i]]), decreasing = TRUE))
}

hapa <- list()
for (i in n) {
  hapa[[i]] <- hap[[i]][o[[i]], ]
}

for (i in n) {
  attr(hapa[[i]], "index") <- idx[[i]][o[[i]]]  ##ordering the haplotypes from most to least freq
}

for (i in n){
  class(hapa[[i]]) <- oc
}

for (i in n) {
  attr(hapa[[i]], "from") <- from
}

#now query the rows within the first and second (most common) haplotype:
row1 <- list()
n <- c(1:length(hapa))
for (i in n) {
  row1[[i]] <- attr(hapa[[i]], "index")[[1]]
}

nonzero <- list() ##produce row2 by first querying haplotypes with >1 length (heterozygous haps only)
n <- c(1:length(hapa))
for (i in n) {
  nonzero[[i]] <- length(attr(hapa[[i]], "index"))
}

row2 <- vector("list", length = length(hapa))
n <- c(1:length(hapa))
for (i in n) {
  row2[[i]] <- ifelse(nonzero[[i]] > 1,
                      list(attr(hapa[[i]], "index")[[2]]), 1)  ##assign an arbitary value (in this case, 1) to second haplotypes so they don't stay NAs and create problems
}

row.2 <- list()
for (i in n) {
  row.2[[i]] <- unlist(row2[[i]])
}

haplotype1 <- row1  #rows in DNA in haplotype1
haplotype2 <- row.2 #rows in DNA in haplotype2

######################################
######################################
###########***ASIDE***################
######################################
######################################

##To interrogate phasing, query the 3rd and 4th haplotypes:

row3 <- vector("list", length = length(hapa))
n <- c(1:length(hapa))
for (i in n) {
  row3[[i]] <- ifelse(nonzero[[i]] > 2,
                      list(attr(hapa[[i]], "index")[[3]]), 1)  ##assign an arbitary value (in this case, 1) to third haplotypes so they don't stay NAs and create problems
}

row.3 <- list()
for (i in n) {
  row.3[[i]] <- unlist(row3[[i]])
}

row4 <- vector("list", length = length(hapa))
n <- c(1:length(hapa))
for (i in n) {
  row4[[i]] <- ifelse(nonzero[[i]] > 3,
                      list(attr(hapa[[i]], "index")[[4]]), 1)  ##assign an arbitary value (in this case, 1) to fourth haplotypes so they don't stay NAs and create problems
}

row.4 <- list()
for (i in n) {
  row.4[[i]] <- unlist(row4[[i]])
}

haplotype3 <- row.3 #rows in DNA in haplotype3
haplotype4 <- row.4 #rows in DNA in haplotype4

##******************PROBLEM:**********************

#Only the most "agregious" paralogues were depth-filtered out in dDocent
#and evidently many alignments with 3rd, 4th and even 5th haplotypes
#have not been filtered out yet. So now develop a
#list of haplotype depths for each individual, and filter out those
#alignments with high-frequency haplotypes exceeding 2 (assuming diploid species)


haplo1 <- lapply(haplotype1, function(x) length(x)) ##frequency (depth) of haplotype1

n <- c(1:length(haplo1))

for (i in n) {
  names(haplo1[[i]]) <- unique(names(dna[[i]]))
}

haplo2 <- lapply(haplotype2, function(x) length(x))
for (i in n) {
  names(haplo2[[i]]) <- unique(names(dna[[i]]))
}

haplo3 <- lapply(haplotype3, function(x) length(x))
for (i in n) {
  names(haplo3[[i]]) <- unique(names(dna[[i]]))
}

haplo4 <- lapply(haplotype4, function(x) length(x))
for (i in n) {
  names(haplo4[[i]]) <- unique(names(dna[[i]]))
}

haplo_1 <- list(haplo1, haplo2, haplo3, haplo4)
save(haplo_1, file = "haplotypes_1.RData")

rm(list=setdiff(ls(), c("dna", "haplotype1", "haplotype2", "row1", "row.2")))
#Clear memory of anything unneccessary

####################################################
#######Subsetting by 1st and 2nd haplotype##########
####################################################

## now subset the reads by haplotype1 and haplotype2 to exclude those containing errors
## NB: first assign the original names to seqs so that sequences can be traced back:

##Files are so large we have to do half and half

dna <- dna[lengths(dna) > 1] ##Before haplotyping, we had to filter out stringsets with length 1
                             ##so now this is necessary to make the dimensions of dna and string match
n <- c(1:length(dna))
hap1 <- list()
for (i in n) {
  hap1[[i]] <- dna[[i]][haplotype1[[i]]]
}

hap_1 <- hap1[lengths(hap1) > 0] ##Exclude those alignments of length zero
rm("hap1")

##Align using muscle and also convert to DNAStringSet:

###############
######1/2######
###############

n <- c(1:(length(hap_1)/2))
aln.hapA <- list()
for (i in n){
  aln.hapA[i] <- DNAStringSet(muscle::muscle(hap_1[[i]]))
}

n <- c(1:length(aln.hapA))
conhap1A <- list()
for (i in n) {
conhap1A[i] <- ConsensusSequence(aln.hapA[[i]], threshold = 0.05, minInformation = 0.75)
}
names1 <- lapply(hap_1[1:length(n)], function(x) unique(names(x)))
for (i in n) {
  names(conhap1A[[i]]) <- paste0(names1[[i]])
}

save(conhap1A, file = "conhap1A.RData")

rm(list=setdiff(ls(), c("hap_1", "dna", "haplotype1", "haplotype2", "row1", "row.2")))

##################
######2/2#########
##################

a <- ((length(hap_1)/2)+1)
b <- length(hap_1)
hapB <- hap_1[a:b]

n <- c(1:length(hapB))
aln.hapB <- list()
for (i in n){
  aln.hapB[[i]] <- DNAStringSet(muscle::muscle(hapB[[i]]))
}

n <- c(1:length(aln.hapB))
conhap1B <- list()
for (i in n) {
  conhap1B[i] <- ConsensusSequence(aln.hapB[[i]], threshold = 0.05, minInformation = 0.75)
}
names1 <- lapply(hapB, function(x) unique(names(x)))
for (i in n) {
  names(conhap1B[[i]]) <- paste0(names1[[i]])
}

save(conhap1B, file = "conhap1B.RData")
load("conhap1A.RData")
load("conhap1B.RData")

rm(list=setdiff(ls(), c("conhap1A", "conhap1B", "row1", "row.2", "dna", "haplotype2")))

conhap1 <- c(conhap1A, conhap1B)

##NB!!
#the consensus sequences based on read depth <3 need to be excluded:

conhap_1 <- vector("list", length = length(conhap1))
n <- c(1:length(conhap1))
for (i in n) {
  conhap_1[[i]] <- conhap1[[i]][length(row1[[i]]) >= 3]
}

##Exclude consensus seqs of length zero
conhap.1 <- conhap_1[lengths(conhap_1) > 0]

## (5) export the consensus sequences

nam1 <- list()
n <- c(1:length(conhap.1))
for (i in n) {
  nam1[[i]] <- names(conhap.1[[i]])
}
nam1 <- unlist(nam1)

write.fasta(sequences = conhap.1, names = paste0(nam1,"_hap1" ,"_1"), nbchar = 1000,
            file.out = "conhap1_1.fasta")

################
################
##rep for hap2##
################
################

n <- c(1:length(dna))
hap2 <- list()
for (i in n) {
  hap2[[i]] <- dna[[i]][haplotype2[[i]]]
}

hap_2 <- hap2[lengths(hap2) > 0] ##Exclude those alignments of length zero
rm("hap2")

###############
######1/2######
###############

n <- c(1:(length(hap_2)/2))
aln.hap2A <- list()
for (i in n) {
  aln.hap2A[[i]] <- DNAStringSet(muscle::muscle(hap_2[[i]]))
}

n <- c(1:length(aln.hap2A))
conhap2A <- list()
for (i in n) {
  conhap2A[i] <- ConsensusSequence(aln.hap2A[[i]], threshold = 0.05, minInformation = 0.75)
}
names2 <- lapply(hap_2[1:length(n)], function(x) unique(names(x)))
for (i in n) {
  names(conhap2A[[i]]) <- paste0(names2[[i]])
}

save(conhap2A, file = "conhap2A.RData")

rm(list=setdiff(ls(), c("hap_2", "row.2")))

###############
######2/2######
###############

a <- ((length(hap_2)/2)+1)
b <- length(hap_2)
hapB <- hap_2[a:b]

n <- c(1:length(hapB))
aln.hap2B <- list()
for (i in n) {
  aln.hap2B[[i]] <- DNAStringSet(muscle::muscle(hapB[[i]]))
}

n <- c(1:length(aln.hap2B))
conhap2B <- list()
for (i in n) {
  conhap2B[i] <- ConsensusSequence(aln.hap2B[[i]], threshold = 0.05, minInformation = 0.75)
}
names2 <- lapply(hapB, function(x) unique(names(x)))
for (i in n) {
  names(conhap2B[[i]]) <- paste0(names2[[i]])
}

save(conhap2B, file = "conhap2B.RData")

load("conhap2A.RData")
load("conhap2b.RData")

conhap2 <- c(conhap2A, conhap2B)

conhap_2 <- vector("list", length = length(conhap2))
n <- c(1:length(conhap2))
for (i in n) {
  conhap_2[[i]] <- conhap2[[i]][length(row.2[[i]]) >= 3]
}

conhap.2 <- conhap_2[lengths(conhap_2) > 0]

nam2 <- list()
n <- c(1:length(conhap.2))
for (i in n) {
  nam2[[i]] <- names(conhap.2[[i]])
}
nam2 <- unlist(nam2)

write.fasta(sequences = conhap.2, names = paste0(nam2,"_hap2" ,"_1"), nbchar = 1000,
            file.out = "conhap2_1.fasta")

################################################
###Read in the fastas and align them############
################################################

library(seqinr)
library(ShortRead)
library(stringr)
library(DescTools)

#read in the haplotype consensus sequence fastas
seqs <- readFasta("")
seq <-  seqs[order(seqs@id), ] ##sort by scaffold

names <- str_split_fixed(seq@id, "_", 2)
names <- unique(names[, 1])

nam <- seq@id

n <- c(1:length(names))
n1 <- list()
for (i in n) {
  n1[[i]] <- nam[nam %like% paste0(names[[i]], "%")]
}                                     ##Setting up the labels for the indiv con seqs

##only relevant to second-last loop

n <- c(1:length(names))
chrom <- list()
for (i in n) {
  chrom[[i]] <- (seq@id) %like% paste0(names[[i]], "%")}  ##query shared scaffold name

chroma <- list()
for (i in n) {
  chroma[[i]] <- seq[chrom[[i]] == TRUE]}   ##filter by shared scaffold name

c <- list()
for (i in n) {
  c[[i]] <- muscle::muscle(chroma[[i]]@sread)}  ##muscle align

cstr <- list()
for (i in n) {
  cstr[[i]] = as(c[[i]], "DNAStringSet")} ##convert to DNAStringSet

for (i in n) {
  names(cstr[[i]]) <- n1[[i]]
}

##Some manipulation to sort by inds rather than haplotype
##so that hap1 and hap2 for an individual appear in sequence

library("gtools")

unlist.names <- list()
n <- c(1:length(cstr))
for (i in n) {
  unlist.names[[i]] <- unlist(names(cstr[[i]]))
}

splitnames <- list()
n <- c(1:length(cstr))
for (i in n) {
  splitnames[[i]] <- as.array(strsplit(unlist.names[[i]], "_"))
}

n <- c(1:length(cstr))
n1 <- list()
for (i in n) {
  n1[[i]] <- lapply(splitnames[[i]], `[[`, 1)
}

n2 <- list()
for (i in n) {
  n2[[i]] <- lapply(splitnames[[i]], `[[`, 2)
}

n3 <- list()
for (i in n) {
  n3[[i]] <- lapply(splitnames[[i]], `[[`, 3)
}

names.ord <- list()
for (i in n) {
  names.ord[[i]] <- paste0(n1[[i]], "_", n3[[i]], "_",n2[[i]])
}

##assign the reformatted names in the same order
for (i in n) {
  names(cstr[[i]]) <- names.ord[[i]]
}

cstr1 <- list()
for (i in n) {
  cstr1[[i]] <- cstr[[i]][order(names(cstr[[i]])), ]
}

for (i in n) {
  writeXStringSet(cstr1[[i]], file= paste0(names[[i]], ".fasta"))
}

save(cstr1, file = "haplotype alignments.RData")

###############################
###############################
##****Manipulate fastas****####
###############################
###############################

# the resulting alignments clearly contain some very problematic sequences
# that are not worth retaining, so let's do some QC:

library(ips)
library(adegenet)
library(ape)
library(psych)
library(Biostrings)
library(seqinr)

files <- list.files(pattern = "*.fasta")
n <- c(1:length(files))
file <- list()
for (i in n) {
  file[[i]] <- readDNAStringSet(files[[i]], "fasta")
}

#********************************************************
## (1) assess no. of phylogenetically informative sites:
#********************************************************

n <- c(1:length(file))
fil <- list()
for (i in n) {
  fil[[i]] <- as.DNAbin(as.matrix(file[[i]]))  #DNAStringSet to matrix to DNAbin
}

n <- c(1:length(fil))
vars <- list()
for (i in n) {
  vars[[i]] <- pis(fil[[i]], what = "abs", use.ambiguities = FALSE)
}

files <- as.list(strsplit(files, ".fasta"))  ##assign names for the no. of PISs in each alignment
for (i in n) {
  names(vars[[i]]) <- files[[i]]
}

var <- data.frame(unlist(vars))
excl1 <- rownames(var)[var[,1] == 0]  ##Candidates for exclusion on account of 0 PIS

##visualise###
pis <- unlist(vars)
hist(pis, main = "Number of PIS per alignment", breaks = 50)
mean(pis)
median(pis)
abline(v = mean(pis), col="red", lwd = 2)
abline(v = median(pis), col = "blue", lwd = 2 )

excl2 <- rownames(var)[var[,1] > 18]  ##Candidates for exlcusion on account of >3 times the median number of PIS

##Wait a minute!
#At some point along this histogram, there's a transition from legitimately
#highly informative alignments to those that are trash. Do a quick
#neighbour-joining analysis to see if these alignments are at all informative

for (i in n){
  vars[[i]] <- vars[[i]][vars[[i]] > 18]
}

var <- unlist(vars)
names.var <- as.list(names(var))

##import only the high PIS files

n <- c(1:length(names.var))
file <- list()
for (i in n) {
  file[[i]] <- read.alignment(file = paste0(names.var[[i]], ".fasta"), format = "fasta")
}

dist <- list() ##make distance matrix for each alignment
for (i in n) {
  dist[[i]] <- dist.alignment(file[[i]], matrix = "identity")
}

nj <- nj(dist[[11]])
plot(nj, cex = 0.5)

#*********************************************************************
## (2) Flag alignments with significant frequency 3rd and 4rd
##haplotypes, then exclude them from the analysis as multicopy loci
#*********************************************************************

haplotypes <- list.files(pattern = "*.RData")
n <- c(1:length(haplotypes))
haplo <- list()
for (i in n) {
  haplo[[i]] <- load(paste0(haplotypes[[i]]))
}

#insert haplo file names:
dp <- list()


##Get the names of the high freq reads in hap3 and hap4
n <- c(1:length(dp))
hap3 <- list()
for (i in n) {
  hap3[[i]] <- unlist(dp[[i]][3])  ##subset to haplotype 3
}

hap3.list <- list()
for (i in n) {
  hap3.list[[i]] <- which(hap3[[i]] >= 3)   ####rows in haplotype 3 exceeding a depth of 3 (Andermann et al 2019 cutoff)
}                                           ####i.e. the alignment numbers which have haplotype 3 frequencies at or above 3

hap3.names <- list()
for (i in n) {
  hap3.names[[i]] <- names(hap3.list[[i]])
}

hap3names <- unlist(hap3.names)
unique(hap3names)               ##names of the alignments with high (>=3) frequency haplotype 3
length(unique(hap3names))
table(hap3names)                ##number of sequenced individuals with high-frequency 3rd haplotypes


hap4 <- list()
for (i in n) {
  hap4[[i]] <- unlist(dp[[i]][4])
}

hap4.list <- list()
for (i in n) {
  hap4.list[[i]] <- which(hap4[[i]] >= 3)
}

hap4.names <- list()
for (i in n) {
  hap4.names[[i]] <- names(hap4.list[[i]])
}

hap4names <- unlist(hap4.names)
unique(hap4names)              ##names of the alignments with high (>=3) frequency haplotype 4
length(unique(hap4names))
table(hap4names)               ##number of sequenced individuals with high-frequency 4th haplotypes

int <- intersect(unique(hap3names), unique(hap4names)) ##candidates for exclusion

#***************************************************************
##Assess alignments with:
# No PIS
# Unusually high (3X median) number of PIS
# High frequency 3rd and 4th haplotypes (paralogous alignments)
#***************************************************************

exclude <- c(excl1, excl2, int)
exclude <- unique(exclude)

load("haplotype alignments.RData")

n <- c(1:length(cstr1))
excl.list <- list()
for (i in n) {
  excl.list[[i]] <- gsub("_.*","", names(cstr1[[i]][1]))
}

diff <- as.list(setdiff(unlist(excl.list), exclude))

n <- c(1:915)
cstr2 <- list()
for (i in n) {
  cstr2[[i]] <- readDNAStringSet(paste0(diff[[i]], ".fasta"), "fasta")
}

save(cstr2, file = "filtered alignments.RData")
for (i in n) {
  writeXStringSet(cstr2[[i]], file= paste0(diff[[i]], ".fasta"))
}

load("filtered alignments.RData")

#######################################################
#######################################################
#####Format the haplotype alignments for DISSECT#######
#######################################################
#######################################################

##include a duplicate of hap1 for "hap2" for homozygous haplotypes, as required by DISSECT
library(stringi)

load("filtered alignments.RData")  #read in the haplotype consensus sequence fastas

hap.nomen <- list()                 ##First produce a list of lists, containing hap names
n <- c(1:length(cstr2))
for (i in n) {
  hap.nomen[[i]] <- names(cstr2[[i]])
}

nomen.un <- list()                  ##Break names to remove hap names for manipulation
for (i in n) {
  nomen.un[[i]] <- str_split_fixed(hap.nomen[[i]], "_h", 2)
}

cstr3 <- cstr2                      ##alignment files with new names
for (i in n) {
  names(cstr3[[i]]) <- nomen.un[[i]][, 1]
}

nomen.uniq <- list()                ##all of the individuals with sequences for each alignment
for (i in n) {
  nomen.uniq[[i]] <- unique(nomen.un[[i]][, 1])
}

dups <- list()
for (i in n) {
  dups[[i]] <- data.frame(table(nomen.un[[i]][, 1]))
}

###Produce lists of names for heterozygotes and homozygotes:

sing <- list()                     ##the individuals with 1 haplotype (i.e. homozygotes)
for (i in n) {
  sing[[i]] <- dups[[i]][dups[[i]]$Freq == 1,]
}

dup <- list()                      ##the individuals with 2 haplotypes (i.e. heterozygotes)
for (i in n) {
  dup[[i]] <- dups[[i]][dups[[i]]$Freq == 2,]
}

singles <- list()                  ##alignments of homozygotes only
for (i in n) {
  singles[[i]] <- cstr3[[i]][sing[[i]]$Var1]
}

twins <- list()                    ##alignments of heterozygotes only
for (i in n) {
  twins[[i]] <- cstr2[[i]][names(cstr3[[i]]) %like% dup[[i]]$Var1]
}

##Duplicate the homozygotes and name the duplicates "hap2"

singles1 <- singles

for (i in n) {
  singles1[[i]] <- unname(singles[[i]])
}

nems <- list()
for (i in n) {
  nems[[i]] <- names(singles[[i]])
}

nems1 <- list()
for (i in n) {
  nems1[[i]] <- paste0(nems[[i]], "_hap2")
}

for (i in n){
  names(singles1[[i]]) <- nems1[[i]]
}

######Add these dups to the original dataset:

calign <- list()
n <- c(1:length(cstr2))
for (i in n) {
  calign[[i]] <- c(cstr2[[i]], singles1[[i]])
}

lst <- list()
for (i in n) {
  lst[[i]] = calign[[i]][order(names(calign[[i]]))]
}

##export

names <- list()
for (i in n) {
  names[[i]] <- str_split_fixed(names(lst[[i]][1]), "_", 2)[, 1]
}

lst2 <- lst
for (i in n) {
  names(lst2[[i]]) <- paste0("Ind", str_split_fixed(names(lst[[i]]), "_", 2)[, 2])
}

for (i in n) {
  writeXStringSet(lst2[[i]], file= paste0(names[[i]], ".fasta"))
}

save(lst2, file = "DISSECT ready.RData")

####Because DISSECT requires that all files contain identical no. of seqs,
####Insert seqs of "---" for inds without captured seqs

n <- c(1:length(lst2))
check <- list()
for (i in n) {
  check[[i]] <- str_split_fixed(names(lst2[[i]]), "_", 2)[, 1]
}

names.uniq <- unique(unlist(check))

names.missing <- list()
for (i in n){
  names.missing[[i]] <- setdiff(names.uniq, check[[i]])
}

n1 <- list()
for (i in n){
  n1[[i]] <- paste0(names.missing[[i]], "_hap1")
}

n2 <- list()
for (i in n){
  n2[[i]] <- paste0(names.missing[[i]], "_hap2")
}

missing <- list()
for (i in n){
  missing[[i]] <- c(n1[[i]], n2[[i]])
}

##make a generic DNAStringSet of arbitrary length, to be trimmed to specific length
seq.len <- list()
n <- c(1:length(lst2))
for (i in n) {
  seq.len[[i]] <- lengths(lst[[i]])[1]
}

dnastring <- DNAStringSet("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")

incl.seq <- list()
for (i in n) {
  incl.seq[[i]] <- rep(dnastring, length(missing[[i]]))
}

##Trim the sequences to specific lengths
seq.incl <- list()
for (i in n) {
  seq.incl[[i]] <- DNAStringSet(strtrim(incl.seq[[i]], seq.len[[i]]))
}

for (i in n) {
  names(seq.incl[[i]]) <- missing[[i]]
}

##Produce the final dataset with equal numbers of seqs

lst.new <- list()
for (i in n) {
  lst.new[[i]] <- c(lst2[[i]], seq.incl[[i]])
}

for (i in n){
  lst.new[[i]] <- lst.new[[i]][names(lst.new[[i]]) != "_hap1"]
}

for (i in n){
  lst.new[[i]] <- lst.new[[i]][names(lst.new[[i]]) != "_hap2"]
}

##Check lengths to ensure correct

lengths(lst.new) #yep

save(lst.new, file = "haplotype alignments_DISSECT ready_with blank sequences.RData")

for (i in n) {
  writeXStringSet(lst.new[[i]], file= paste0(names[[i]], ".fasta"))
}

#############################
###Querying matching SNPs####
#############################

load("haplotype alignments_DISSECT ready_with blank sequences.RData")
names <- list.files(pattern = "*.fasta")
names <- gsub('.{6}$', '', names)

library("LEA")
library("plyr")
library("vcfR")
library("adegenet")
library("SNPRelate")
library("stringr")
library("dartR")
library("Biostrings")
library("ips")

vcf <- read.vcfR("DP3g95p5maf05.fil5.FIL.prim.vcf")

#calculate missingness

f <- as.list(vcf@fix[,8])
ns <- list()
n <- c(1:length(f))
for (i in n) {
  ns[[i]] <- strsplit(f[[i]], ";")
}

nsd <- list()
for (i in n) {
  nsd[[i]] <- data.frame(ns[[i]])
}

nsf <- list()
for (i in n) {
  nsf[[i]] <- nsd[[i]][18, ]
}

df <- data.frame(unlist(nsf))
char <- as.list(as.character(df[,1]))

fin <- list()
for (i in n) {
  fin[[i]] <- substring(char[[i]], 4)
}
  
df <- data.frame(unlist(fin))

c <- as.character(df[,1])

a <- c[nchar(c) == 3]
b <- c[nchar(c) > 3]
b <- substring(b, 1, 3)

c <- data.frame(as.numeric(c(a, b)))
mean(c[,1])
  
#####

vcf <- vcf[, which(colnames(vcf@gt) %in% c())]

vcfr <- subset(vcf, vcf@fix[,1] %in% c(names))

vcfR::write.vcf(vcfr, file = "filter.vcf")

##filter out indels and >2 state SNPs

#import .vcf:
vcf <- vcfR::read.vcfR("filter.vcf")
vcf_sansindel <- extract.indels(vcf, return.indels = FALSE)

##which ones did contain indels?
vcf_test <- extract.indels(vcf, return.indels = TRUE)
vcf_test@fix[,1]

##there's no automatic tri and quad-snp removal, so do some
##querying of allele types

maf1 <- maf(vcf_sansindel, element = 1) #first allele
maf2 <- maf(vcf_sansindel, element = 2) #second allele
maf3 <- maf(vcf_sansindel, element = 3) #third allele ##to exclude

vcf_sansindel@fix[,3] <- rownames(maf1)

multisnp <- rownames(maf3)[maf3[,3] > 0] ##3- and 4-state SNP names
chroms <- unique(str_split_fixed(multisnp, "_", 2)[,1]) ##chroms that contain multi-state SNPs

l <- unique(vcf_sansindel@fix[,1])
incl.chrom <- setdiff(l, chroms)  ##chroms w/o 3/4-state SNPs
vcf_IDMS <- subset(vcf_sansindel, vcf_sansindel@fix[,1] %in% c(incl.chrom))       ##InDels and MultiSnps excluded

##Great! now select a single biallelic SNP from each chrom for BFD (reduce linkage)

df <- cbind(data.frame(vcf_IDMS@fix[,1]), data.frame(vcf_IDMS@fix[,3]))
names(df) <- c("chrom", "snp")
df_split <- split(df, df$chrom)

##now we have a list of chroms, along with the snps they contain
##select one snp randomly using sample()

sample <- list()
n <- c(1:length(df_split))
for (i in n) {
  sample[[i]] <- sample(df_split[[i]][,2], 1)
}

sample <- lapply(sample, function(x) as.character(x))
sample <- unlist(sample)

###Subset VCF

vcf_sub <- subset(vcf_IDMS, vcf_IDMS@fix[,3] %in% c(sample))
vcfR::write.vcf(vcf_sub, file = "vcf_indel_biallelic.vcf")

###produce the fastas that correspond with these snps

sample_chroms <- as.list(vcf_sub@fix[,1])

n <- c(1:length(sample_chroms))
file <- list()
for (i in n) {
  file[[i]] <- readDNAStringSet(paste0(sample_chroms[[i]], ".fasta"))
}

for (i in n) {
  writeXStringSet(file[[i]], file= paste0(sample_chroms[[i]], ".fasta"))
}

save(file, file = "fastas_indel_biallelic.RData")


####Removing paralogues

library(Biostrings)
library(phangorn)
library(ape)
library(Rcpp)
library(stats)
library(DECIPHER)
library(gdata)
library(seqinr)

##First query whether the 1st and 2nd haps for an ind fall in the same
##clade in a upgma tree. We expect paralogs to fall in different clades  

files <- list.files(pattern = "*.fasta")[c(1:810, 812:816)]
n <- c(1:length(files))
file <- list()
for (i in n) {
  file[[i]] <- read.phyDat(files[[i]], format="fasta", type = "DNA")
}

dist <- list()
n <- c(1:length(file))
for (i in n) {
  dist[[i]] <- dist.ml(file[[i]])
}

tree <- list()
n <- c(1:length(dist))
for (i in n) {
  tree[[i]] <- upgma(dist[[i]])
}

plot(tree[[1]], cex = 0.4)  
nodelabels()
tree[[1]]$edge
a <- extract.clade(tree[[1]], 74)
plot(a)   ##Query the two clades on either side of the root node

node <- list()
n <- c(1:length(tree))
for (i in n) {
  node[[i]] <- data.frame(tree[[i]]$edge)
}

length <- list()
n <- c(1:length(node))
for (i in n) {
  length[[i]] <- dim(node[[i]])[1]  ##number of items in node element of tree
}

select <- list()
for (i in n) {
  select[[i]] <- node[[i]][rownames(node[[i]]) == length[[i]], ]
}

childnodes <- list()
for (i in n) {
  childnodes[[i]] <- node[[i]][node[[i]][,1] == select[[i]][1,1], ] 
}

##Clade 1 
c1 <- lapply(childnodes, function(x) x[1,2])

##Clade 2
c2 <- lapply(childnodes, function(x) x[2,2])

##Query the second clade (first clade frequently has one tip - the outgroup - which cannot be queried)

clade2 <- list()
n <- c(1:length(tree))
for (i in n) {
  clade2[[i]] <- extract.clade(tree[[i]], c2[[i]])
}

##the extract.clade function has rearranged the node numbers
##but otherwise correct

plot(tree[[1]], cex = 0.4)
nodelabels()
plot(clade2[[1]], cex = 0.4)
nodelabels()    ##looks correct

##Since all the trees have the exact same taxa, use setdiff
##to query the tiplabels in clade1 that we couldn't get before

names <- tree[[1]]$tip.label

ndiff <- list()
n <- c(1:length(clade2))
for (i in n) {
  ndiff[[i]] <- setdiff(names, clade2[[i]]$tip.label)
}

nint <- lapply(clade2, function(x) x$tip.label)


ndiff  ##Nice!! this is a list of names in clade 1
nint   ##list of names in clade 2

##Break the names to exclude _hap1 and _hap2

haps <- c("hap1", "hap2")

splitdiff <- list()
n <- c(1:length(clade2))
for (i in n) {
  splitdiff[[i]] <- unlist(strsplit(ndiff[[i]], "_"))
}
splitdiff <- lapply(splitdiff, function(x) setdiff(x, haps))

splitint <- list()
n <- c(1:length(clade2))
for (i in n) {
  splitint[[i]] <- unlist(strsplit(nint[[i]], "_"))
}
splitint <- lapply(splitint, function(x) setdiff(x, haps))

clad1 <- splitdiff
clad2 <- splitint

int <- list()
for (i in n) {
  int[[i]] <- intersect(clad1[[i]], clad2[[i]])
}

##item "int" above is the names of individuals 
##whose 1st and 2nd haplotypes fall into two differenct clades

para <- lengths(int)
df <- cbind(data.frame(para), data.frame(files))

##Exclude alignments with 0 overlaps - "above suspicion" for paralogy

df <- df[df$para > 0, ]  ##661 alignments remain to be interrogated
hist(df$para, breaks = 50, main = "", ylab = "No. of alignments (/815)" ,xlab = "No. of inds with non-monophyletic hap1 and hap2")

df10 <- df[df$para > 9, ]  ##highly suspicious
df7 <- df[df$para > 6, ]  ##moderately suspicious
df5 <- df[df$para > 4, ]  ##moderately suspicious
df3 <- df[df$para > 2, ]  ##somewhat suspicious

##import and export excl. df5 (lenient)
gooddf5 <- df[df$para < 5, ]
files <- as.list(as.character(gooddf5$files))
n <- c(1:length(files))
file <- list()
for (i in n) {
  file[[i]] <- readDNAStringSet(file = files[[i]], format = "fasta")
}

for (i in n) {
  writeXStringSet(file[[i]], file = files[[i]])
}

##import and export excl. df3 (strict)
gooddf3 <- df[df$para < 3, ]
files <- as.list(as.character(gooddf3$files))
n <- c(1:length(files))
file <- list()
for (i in n) {
  file[[i]] <- readDNAStringSet(file = files[[i]], format = "fasta")
}

for (i in n) {
  writeXStringSet(file[[i]], file = files[[i]])
}

##filter vcf by alignments obtained via UPGMA filtering:

# original vcf name 
vcf <- read.vcfR()

filenames <- as.list(list.files(pattern = "*.fasta"))  ##have to subset
files <- unlist(strsplit(unlist(filenames), ".fasta"))

vcf_strict <- vcf[which(vcf@fix[,1] %in% files)]
bi <- is.biallelic(vcf_strict)
vcf_strict <- vcf_strict[which(bi == TRUE), ]
vcfR::write.vcf(vcf_strict, file = "vcf_strict.vcf")

strict_excl <-  files[which(bi == FALSE)]

filenames <- as.list(list.files(pattern = "*.fasta"))  ##have to subset
files <- unlist(strsplit(unlist(filenames), ".fasta"))

vcf_lenient <- vcf[which(vcf@fix[,1] %in% files)]
bi <- is.biallelic(vcf_lenient)
vcf_lenient <- vcf_lenient[which(bi == TRUE), ]
vcfR::write.vcf(vcf_lenient, file = "vcf_lenient.vcf")

lenient_excl <- files[which(bi == FALSE)]

##Convert to nexus
vcf2SNAPP(vcf_lenient, file = "lenient.nex")

