##  ###################################################  ##
##  Processing raw ITS seqs into a phyloseq object       ##
##  This script processes 2210KMI-0004                   ##
##                                                       ##
##  Mangrove bacteria and archaea microbiome             ##
##                                                       ##
##  Author: Ming Sheng - 12 Nov, 2022                    ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.2.0                                            ##
##  dada2 v 1.24.0                                       ##
##  purrr v 0.3.4                                        ##
##  tidyverse v 1.3.1                                    ##
##  readxl v 1.4                                         ##
##  decontam v 1.16                                      ##
##  phyloseq v 1.40.0                                    ##
##  ShortRead v 1.54.0                                   ##
##  ###################################################  ##


# Load packages ####
library(dada2); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(tidyverse); packageVersion("tidyverse")
library(readxl); packageVersion("readxl")
library(ShortRead); packageVersion("ShortRead")

# Load metadata ####
full_meta <- readxl::read_xlsx("16S_MetaData.xlsx")

# Remove extra cols
names(full_meta)[1:1] <- 'SampleID'

# Add column identifying PCR negatives
full_meta$PCR_Negative <- FALSE
full_meta$PCR_Negative[grep("BLANK",full_meta$SampleID)] <- TRUE

# Find raw fastq files and prepare workspace ####
path <- "./2210KMI-0004"
fqs <- list.files(path, full.names = TRUE, recursive = FALSE, pattern = "_R1.fastq.gz$|_R2.fastq.gz$")

# Parse fwd and rev reads
fnFs <- fqs[grep("_R1.fastq.gz",fqs)]
fnRs <- fqs[grep("_R2.fastq.gz",fqs)]

# Get Sample Names
sample.names <- str_remove(basename(fnFs),"_R1.fastq.gz")

# subset metadata to this run's samples
meta <- full_meta %>% filter(SampleID %in% sample.names)
# quick cleanup of environment

# Peek at quality profiles
plotQualityProfile(fnFs[c(1,30)]) # fwd reads # drops below Q30 at around 250
plotQualityProfile(fnRs[c(10,150)]) # rev reads # drops below Q30 at around 150

# Make filtered outfile names
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# make new directory for filtered files
if(!dir.exists(file.path(path,"filtered"))){
  dir.create(file.path(path,"filtered"))
}

# check for duplicated sample names
sum(duplicated(sample.names))

# CHECK FOR AND REMOVE PRIMER SITES WITH CUTADAPT ####
FWD <- "GTGYCAGCMGCCGCGGTAA" # Sequence of 515F FWD primer 
REV <- "GGACTACNVGGGTWTCTAAT" # Sequence of 806R REV primer

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients; REV.orients


# Prefilter to remove reads with ambiguous (N) bases ####
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE) # on Windows, set multithread = FALSE

# Discover primer matches, regardless of orientation ####
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[10]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[10]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[10]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[10]]))

# Run cutadapt ####
# If the following command returns an error, you do not have cutadapt installed correctly
cutadapt <- "./cutadapt" # download the executable file from https://github.com/marcelm/cutadapt/releases and run it directly
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "--minimum-length 100", # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# sanity check
# This should show no occurrences in any of the orientations now
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[10]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[10]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[10]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[10]]))

# Filter and trim ####

# cut fwd reads at 250 and rev reads at 150
out <- filterAndTrim(fnFs.cut, filtFs, fnRs.cut, filtRs, truncLen=c(250,150),
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE,verbose = TRUE)

#Check trimmed sequence quality
plotQualityProfile(filtFs[c(20,30)]) 
plotQualityProfile(filtRs[c(20,30)])

saveRDS(out,"./Output/out.RDS") # save object
out <- readRDS("./Output/out.RDS") # reload point
##HG06 has 9 reads.in but 0 reads.out

filtpath <- file.path(path,"filtered")

## get sample names again
sample.names <- str_remove(basename(fnFs),"_R1.fastq.gz")

# reassign filts for any potentially lost samples
filtFs <- list.files("./2210KMI-0004/filtered", pattern = "_F_filt", full.names = TRUE)
filtRs <- list.files("./2210KMI-0004/filtered", pattern = "_R_filt", full.names = TRUE)

# Learn error rates ####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errF,"./Output/errF.RDS") # save object
saveRDS(errR,"./Output/errR.RDS") # save object

errF <- readRDS("./Output/errF.RDS") # reload point
errR <- readRDS("./Output/errR.RDS") # reload point
out <- readRDS("./Output/out.RDS")

# plot error rates for sanity
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Infer sequence variants ####

# add names to filts
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Dereplication, sample inferrence, and merging ####

# loop through each pair of fwd and rev reads, one file at a time
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}

saveRDS(mergers,"./Output/mergers_1_15.RDS") #mergers were conducted in HPC

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table ####
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab,"./Output/seqtab.RDS") # save object
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Save progress
saveRDS(seqtab.nochim,"./Output/seqtab.nochim.RDS")
seqtab.nochim = readRDS("./Output/seqtab.nochim.RDS")

# Assign taxonomy #### CONSTAXv2 was used instead https://github.com/liberjul/CONSTAXv2
taxa <- assignTaxonomy(seqtab.nochim,"./silva_nr99_v138.1_train_set.fa.gz", minBoot = 80,multithread = TRUE)
saveRDS(taxa,"./Output/taxa.RDS")

taxa <- addSpecies(taxa, "./silva_species_assignment_v138.1.fa.gz")
saveRDS(taxa,"./Output/taxa_with_spp.RDS")


# rename seqtab object samples
seqtab.df <- as.data.frame(seqtab.nochim)
row.names(seqtab.df)

# Create phyloseq object ####

# subset to remove missing samples
in.meta <- which(names(seqtab.nochim[,1]) %in% full_meta$SampleID == TRUE)
seqtab.nochim <- seqtab.nochim[in.meta,]
dim(seqtab.nochim)
in.seqtab <- which(full_meta$SampleID %in% names(seqtab.nochim[,1]))
meta <- meta[in.seqtab,]

# re-order
full_meta <- full_meta[order(full_meta$SampleID),]
row.names(full_meta) <- full_meta$SampleID
seqtab.nochim <- (seqtab.nochim[row.names(full_meta),])
identical(row.names(seqtab.nochim), as.character(full_meta$SampleID))

# make phyloseq object
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
met <- sample_data(full_meta)
tax <- tax_table(taxa)

sample_names(met) <- met$SampleID
ps <- phyloseq(otu,met,tax)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# save it
saveRDS(ps, "./Output/raw_ps_object.RDS")


# Identify and remove contaminants ####
library(decontam)
ps <- readRDS("./Output/raw_ps_object.RDS")

blanks = which(ps@sam_data$PCR_Negative == TRUE)
contamdf.prev <- isContaminant(ps, neg=blanks, threshold = 0.01)
table(contamdf.prev$contaminant)

ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps) #no contam

# save contam-free phyloseq object
saveRDS(ps.noncontam, "./Output/noncontam_ps_object.RDS")

ps_sp <- readRDS("./Output/noncontam_ps_object.RDS")

# Sequencing reads produced and processed ####
out %>% as.data.frame() %>% select(reads.in) %>% sum()

no.chim <- readRDS("./Output/seqtab.nochim.RDS")
run <- out %>% 
  as.data.frame() %>% 
  mutate(no.chim = no.chim %>% rowSums(),sampleID = row.names(.))

write_csv(run, "Output/Stats/ReadStats.csv")

## LULU algorithm https://github.com/tobiasgf/lulu
devtools::install_github('tobiasgf/lulu')
library(lulu)

ps <- readRDS('./Output/noncontam_ps_object.RDS')

# Export all sequences to FASTA to be converted to .txt file
library(ShortRead)
fasta_dir <- file.path(getwd(), 'Output/refs/')
outfile <- file.path(dirname(fasta_dir), "all.fasta")
writeFasta(ps@refseq, outfile, mode = 'a')

# Export OTU table as text file (Samples as columns, ASVs as rows): to produce match list with blastn
write.table(t(ps@otu_table), './Output/otu_table.txt', sep = '\t', quote = F)

# run LULU algorithm
otutab <- as.data.frame(t(ps@otu_table))
matchlist <- read.table("./Output/sed_match_list.txt", header = FALSE, as.is = T, stringsAsFactors = FALSE)

curated_result <- lulu(otutab, matchlist)

curated_result$curated_table 

new_tax_table <- ps@tax_table %>% as.data.frame %>% .[c(curated_result$curated_otus),] %>% tax_table()
rownames(new_tax_table) <- ps@tax_table %>% as.data.frame %>% .[c(curated_result$curated_otus),] %>% rownames()
colnames(new_tax_table) <- ps@tax_table %>% as.data.frame %>% .[c(curated_result$curated_otus),] %>% colnames()

new_otu_tab <- otu_table(t(curated_result$curated_table), taxa_are_rows = FALSE)
met <- sample_data(full_meta)
sample_names(met) <- met$`Sample ID`

ps_object <- phyloseq(new_otu_tab,met,new_tax_table)

saveRDS(ps_object, './Output/lulu_raw_ps_object.RDS')
