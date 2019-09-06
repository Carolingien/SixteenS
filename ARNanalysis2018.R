# repertoire principal
rep = "C:/Users/caroline.doose/Google Drive/2-Th année 2018/ManipBiof-2018/résultats/Stat/"

#sur ordi portable
rep = "C:/Users/Lenovo/Google Drive/2-Th année 2018/ManipBiof-2018/résultats/Stat/"

BiocManager::install("BiocStyle")
BiocManager::install(c("ShortRead", "phyloseq"))
BiocManager::install('dada2')
BiocManager::install("DECIPHER")
BiocManager::install("phangorn")

library('dada2')
library("knitr")
library("BiocStyle")
library("phangorn")
library("DECIPHER")

# Load packages into session, and print package version

sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

set.seed(100)

miseq_path <- "./MiSeq_SOP" # CHANGE to the directory containing the fastq files after unzipping.
list.files(miseq_path)


# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
fnFs[1:3]


plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

## define the filenames for the filtered fastq.gz files:

filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,210),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF)
plotErrors(errR)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]


mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
dim(seqtabAll)

table(nchar(getSequences(seqtabAll)))

seqtabNoC <- removeBimeraDenovo(seqtabAll)

### Assign taxonomy

fastaRef <- "C:/Users/Lenovo/Desktop/ARN16s-Stat/silva_nr_v132_train_set.fa"
taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)
taxTab <- addSpecies(taxTab, "C:/Users/Lenovo/Desktop/ARN16s-Stat/silva_species_assignment_v132.fa.gz")
unname(head(taxTab))

seqs <- getSequences(seqtabNoC)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)


samples.out <- rownames(seqtabNoC)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
Time <- substr(subject,4,5)
Concentration <- substr(subject,7,8)
#day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(subject=subject,Time=Time, Concentration=Concentration)
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxTab))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample


dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

PRCH=plot_richness(ps, x="Time", measures=c("Shannon", "Simpson"), color="Concentration")
PRCH+theme_classic()

svg("ARN16SDiv2018")
PRCH+scale_color_manual(values = c("#99CCCC","#009999", "#00333F"))
dev.off()

is.na(ps.prop)
na.omit(ps.prop)


# Transform data to proportions as appropriate for Bray-Curtis distances ##### MARCHE PASSSS
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")


### Bar plot:
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:100]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
PBF=plot_bar(ps.top20, x="Time", fill="Family") + facet_wrap(~Concentration, scales="free_x")

svg("ARN16SFamille20181002")
PBF+theme_classic()
dev.off()
