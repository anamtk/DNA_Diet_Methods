#Installation Needs####
#This was drawn from both Happy Belly Informatics 
#https://astrobiomike.github.io/amplicon/dada2_workflow_ex

#First, make sure that Xcode (mac C++ functionality) was updated in this first
#in terminal: xcode-select --install

#installs BiocManager from BioConductor 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.10")

#I had to load these two packages manually after I realized that Xcode was not up to date. 
#However if you are sure that Xcode is up to date first, you can probably just run 
#the BiocManager updater and the dada2 installer below without this issue. 
#BiocManager::install("XML", update = TRUE, ask = FALSE)
#BiocManager::install("RCurl", update = TRUE, ask = FALSE)

#this code checks that all packages are up  to date
BiocManager::valid()

#installs dada2
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2")

#load dada2 and get package version 
library(dada2); packageVersion("dada2")
library(tidyverse) #for data cleaning and manipulation
#remotes::install_github("HuntsmanCancerInstitute/hciR")
library(hciR) #as_matrix code
library(ShortRead)
library(here)
#from Happy Belly Bioinformatics code chunk:
list.files() # make sure what we think is here is actually here

#dada2 on dec2019 run####
## first we're setting a few variables we're going to use ##
# one with all sample names, by scanning our "samples" file we made earlier
#if you don't have this variable, you can make it in the terminal with this
#code, just make sure you're in that working directory
#ls *_R1_*.fq.gz | cut -f1-2 -d "_" > samples
getwd()
setwd(here("data", "trimmed", "dec2019"))
#now we can create the samples variable, this should be 60 samples
samples <- scan("samples", what = "character")

#now we create two files hodling forward and reverse reads
# one holding the file names of all the forward reads
forward_reads <- paste0(samples, "_sub_R1_trimmed.fq.gz")
# and one with the reverse
reverse_reads <- paste0(samples, "_sub_R2_trimmed.fq.gz")

# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads <- paste0(samples, "_sub_R1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_sub_R2_filtered.fq.gz")

#plotQualityProfile(forward_reads)
#plotQualityProfile(reverse_reads)
# and just plotting the last 4 samples of the reverse reads
plotQualityProfile(reverse_reads[17:20])
plotQualityProfile(forward_reads[17:20])

#Happy Belly uses the info above to use the truncLen argument
#in the following filterAndTrim function. I am not going to do that just yet
#the original code looked like:
#filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
#                              reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
#                              rm.phix=TRUE, minLen=175, truncLen=c(250,200))

#i'm going to do maxEE at 1, but may reconsider later 
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads, maxEE=c(1,1),
                              rm.phix=TRUE, minLen = 100, multithread = TRUE,
                              matchIDs = TRUE)

filtered_out
#one of my negatives got filtered to zero so I need to delete that for the next
#step - it's throwing it off

#this would probably look better if I included the truncLen argument in the filtering
#step, however, I don't want to do that just yet. Also - My bp region is ~363 bp long,
#so would have to truncate >= to 182 - I feel like truncating at all may be a bad idea
plotQualityProfile(filtered_reverse_reads[17:20])

#delete samples with zero outputs post filtering for the sequence table below
filtered_out1 <- filtered_out %>%
  as_tibble(rownames = "sample") %>%
  filter(reads.out > 0) %>%
  as_matrix() #from the hciR package for this purpose!
filtered_out1

#need to find the samples with zero values in reads.out. This pipe does that:
empty <- filtered_out %>%
  as_tibble(rownames = "sample") %>%
  filter(reads.out <= 0) %>%
  dplyr::select(sample) %>%
  as_vector() %>%
  unname() %>%
  str_sub(end=-22)

empty #which samples are empty 

#now remove them from forward and reverse, first by creating a vector with the name of that sample
empty_forward <- paste0(empty, "_sub_R1_filtered.fq.gz")
empty_reverse <- paste0(empty, "_sub_R2_filtered.fq.gz")

#then removing from the filtered forward reads file
filtered_forward_reads1 <- filtered_forward_reads[! filtered_forward_reads %in% empty_forward]

#and reverse reads
filtered_reverse_reads1 <- filtered_reverse_reads[! filtered_reverse_reads %in% empty_reverse]

#and samples file
samples1 <- samples[! samples %in% empty]
#error model of forward and reverse reads
#the truncate length above may alter this error structure, something to consider
err_forward_reads <- learnErrors(filtered_forward_reads1, randomize = TRUE,
                                 multithread = TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads1, randomize = TRUE,
                                 multithread = TRUE)

plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)



#dereplicate the forward and reverse reads now that we have removed the sample with zero filter
derep_forward <- derepFastq(filtered_forward_reads1, verbose=TRUE)
names(derep_forward) <- samples1 # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtered_reverse_reads1, verbose=TRUE)
names(derep_reverse) <- samples1

#this is the dada ASV assigning step, which incorporates both abundance and error
dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread = TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread = TRUE)

#this is the merging step, which merges the forward and reverse with default overlap of 12bp. I can't
#remember what my overlap should be but the original is here and I'm removing the minOverlap command
#merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
#                               derep_reverse, trimOverhang=TRUE, minOverlap=170)

merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, trimOverhang=TRUE)

# this object holds a lot of information that may be the first place you'd want to look if you want to start poking under the hood
class(merged_amplicons) # list
length(merged_amplicons) # 20 elements in this list, one for each of our samples
names(merged_amplicons) # the names() function gives us the name of each element of the list 

#these did not work for me, not sure why
#class(merged_amplicons$B1) # each element of the list is a dataframe that can be accessed and manipulated like any ordinary dataframe
#names(merged_amplicons$B1) # the names() function on a dataframe gives you the column names
# "sequence"  "abundance" "forward"   "reverse"   "nmatch"    "nmismatch" "nindel"    "prefer"    "accept"

#this is a sequence table including chimeras, or OTU table
seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab) # matrix
dim(seqtab) # 59 245

#this detects chimeras, and removes them
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T) # Identified 31 bimeras out of 245 input sequences.

#this command shows how much read abundance we lost by removing chimeras
# though we only lost 54 sequences, we don't know if they held a lot in terms of abundance, this is one quick way to look at that
sum(seqtab.nochim)/sum(seqtab) # 0.999612 # good, we barely lost any in terms of abundance

# set a little function
getN <- function(x) sum(getUniques(x))

#NEED TO UPDATE TO REMOVE THE ONE NEGATIVE
# making a little table
summary_tab <- data.frame(row.names=samples1, dada2_input=filtered_out1[,1],
                          filtered=filtered_out1[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out1[,1]*100, 1))

summary_tab

#assign taxonomy in dada2. i think this option requires downloading some fastas from internet
#taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", tryRC=T)
# taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=TRUE, tryRC=T)

#extracting some useful things from dada2
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa") #this goes to the samples folder - will need to move it

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA) #also goes to samples folder, need to move it!

# tax table:
#asv_tax <- taxa
#row.names(asv_tax) <- sub(">", "", asv_headers)
#write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)
#dada2 on both runs####
## first we're setting a few variables we're going to use ##
# one with all sample names, by scanning our "samples" file we made earlier
samples <- scan(here("data", "samples"), what = "character")

#move us into that data folder for next steps, since that seems like 
#the only way these will work
setwd(here("data"))
# one holding the file names of all the forward reads
forward_reads <- paste0(samples, "_sub_R1_trimmed.fq.gz")
# and one with the reverse
reverse_reads <- paste0(samples, "_sub_R2_trimmed.fq.gz")

# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads <- paste0(samples, "_sub_R1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_sub_R2_filtered.fq.gz")

#plotQualityProfile(forward_reads)
#plotQualityProfile(reverse_reads)
# and just plotting the last 4 samples of the reverse reads
plotQualityProfile(reverse_reads[17:20])
plotQualityProfile(forward_reads[17:20])

#Happy Belly uses the info above to use the truncLen argument
#in the following filterAndTrim function. I am not going to do that just yet
#the original code looked like:
#filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
#                              reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
#                              rm.phix=TRUE, minLen=175, truncLen=c(250,200))

#i'm going to do maxEE at 1, but may reconsider later 
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads, maxEE=c(1,1),
                              rm.phix=TRUE, minLen = 100, multithread = TRUE,
                              matchIDs = TRUE)

filtered_out
#one of my negatives got filtered to zero so I need to delete that for the next
#step - it's throwing it off

#this would probably look better if I included the truncLen argument in the filtering
#step, however, I don't want to do that just yet. Also - My bp region is ~363 bp long,
#so would have to truncate >= to 182 - I feel like truncating at all may be a bad idea
plotQualityProfile(filtered_reverse_reads[17:20])

#error model of forward and reverse reads
#the truncate length above may alter this error structure, something to consider
err_forward_reads <- learnErrors(filtered_forward_reads, multithread = TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread = TRUE)

plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)


#delete negative for the sequence table below
filtered_out1 <- filtered_out %>%
  as_tibble(rownames = "sample") %>%
  filter(reads.out > 0) %>%
  as_matrix() #from the hciR package for this purpose!

#filtered_out %>%
#  as_tibble(rownames = "sample")

#need to find the samples with zero values in reads.out. This pipe does that:
empty <- filtered_out %>%
  as_tibble(rownames = "sample") %>%
  filter(reads.out <= 0) %>%
  select(sample) %>%
  as_vector() %>%
  unname() %>%
  str_sub(end=-19)

empty #check the output

#create vector of empty sample names
empty_forward <- paste0(empty, "_sub_R1_filtered.fq.gz")
empty_reverse <- paste0(empty, "_sub_R2_filtered.fq.gz")

#remove them from the filtered foward and reverse reads
filtered_forward_reads1 <- filtered_forward_reads[! filtered_forward_reads %in% empty_forward]

filtered_reverse_reads1 <- filtered_reverse_reads[! filtered_reverse_reads %in% empty_reverse]

#remove these samples from sample file too
samples1 <- samples[! samples %in% empty]
#dereplicate the forward and reverse reads now that we have removed the sample with zero filter
derep_forward <- derepFastq(filtered_forward_reads1, verbose=TRUE)
names(derep_forward) <- samples1 # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtered_reverse_reads1, verbose=TRUE)
names(derep_reverse) <- samples1

#this is the dada ASV assigning step, which incorporates both abundance and error
dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread = TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread = TRUE)

#this is the merging step, which merges the forward and reverse with default overlap of 12bp. I can't
#remember what my overlap should be but the original is here and I'm removing the minOverlap command
#merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
#                               derep_reverse, trimOverhang=TRUE, minOverlap=170)

merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, trimOverhang=TRUE)

# this object holds a lot of information that may be the first place you'd want to look if you want to start poking under the hood
class(merged_amplicons) # list
length(merged_amplicons) # 20 elements in this list, one for each of our samples
names(merged_amplicons) # the names() function gives us the name of each element of the list 

#these did not work for me, not sure why
#class(merged_amplicons$B1) # each element of the list is a dataframe that can be accessed and manipulated like any ordinary dataframe
#names(merged_amplicons$B1) # the names() function on a dataframe gives you the column names
# "sequence"  "abundance" "forward"   "reverse"   "nmatch"    "nmismatch" "nindel"    "prefer"    "accept"

#this is a sequence table including chimeras, or OTU table
seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab) # matrix
dim(seqtab) # 97 353

#this detects chimeras, and removes them
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T) # Identified 54 bimeras out of 353 input sequences.

#this command shows how much read abundance we lost by removing chimeras
# though we only lost 54 sequences, we don't know if they held a lot in terms of abundance, this is one quick way to look at that
sum(seqtab.nochim)/sum(seqtab) # 0.999612 # good, we barely lost any in terms of abundance

# set a little function
getN <- function(x) sum(getUniques(x))

#NEED TO UPDATE TO REMOVE THE ONE NEGATIVE
# making a little table
summary_tab <- data.frame(row.names=samples1, dada2_input=filtered_out1[,1],
                          filtered=filtered_out1[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out1[,1]*100, 1))

summary_tab

#assign taxonomy in dada2. i think this option requires downloading some fastas from internet
#taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", tryRC=T)
# taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=TRUE, tryRC=T)

#extracting some useful things from dada2
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
#asv_tax <- taxa
#row.names(asv_tax) <- sub(">", "", asv_headers)
#write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)