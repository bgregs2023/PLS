# Install Bioconductor manager if not already installed
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")  # Installs BiocManager if not already present

# Install packages from Bioconductor
BiocManager::install(c("dada2", "phyloseq", "ShortRead"))  # Install bioinformatics packages

# Install CRAN packages
install.packages(c("tidyverse","vegan", "ggplot2", "dplyr"))  # Install data manipulation and plotting tools

# Load all the necessary packages
library(dada2)       # For amplicon sequence processing
library(ShortRead)   # For reading and handling FASTQ files
library(phyloseq)    # For microbial ecology analysis and plotting
library(ggplot2)     # For data visualization
library(dplyr)       # For data manipulation

path <- "~/Documents/PLS/16S/"  # CHANGE ME to the directory containing the 16S fastq files
list.files(path)  # Check the contents of the directory

# Forward and reverse fastq filenames have format
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))  # Forward reads
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))  # Reverse reads

# Plot quality profile for first natural sample (forward and reverse)
plotQualityProfile(fnFs[1])  # Visualize quality of the first forward read
plotQualityProfile(fnRs[1])  # Visualize quality of the first reverse read

# Plot quality profile for first urban sample (forward and reverse)
plotQualityProfile(fnFs[4])  # Visualize quality of a later sample (forward)
plotQualityProfile(fnRs[4])  # Visualize quality of a later sample (reverse)

# Split filenames at underscores to extract sample names
parts <- strsplit(basename(fnFs), "_")
sample.names <- sapply(parts, function(x) paste(x[1:3], collapse = "_"))  # Create sample names
sample.names  # View sample names

# Create paths for your filtered outputs to be saved on your machine
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))  # Filtered forward reads
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))  # Filtered reverse reads
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Perform filtering and trimming based on quality thresholds
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

# Dereplicate forward and reverse reads (collapse identical sequences)
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Define file paths to the first filtered forward and reverse reads
filt_file_F <- filtFs[1]
filt_file_R <- filtRs[1]

# Set a seed for reproducibility
set.seed(42)

# Subsample Forward Reads (5%) to speed up error learning
fq_F <- readFastq(filt_file_F)
subsampled_fq_F <- fq_F[sample(1:length(fq_F), length(fq_F) * 0.05)]
subsampled_path_F <- tempfile(fileext = ".fastq.gz")  # Create temporary output path
writeFastq(subsampled_fq_F, subsampled_path_F, compress = TRUE)

# Subsample Reverse Reads (5%) to speed up error learning
fq_R <- readFastq(filt_file_R)
subsampled_fq_R <- fq_R[sample(1:length(fq_R), length(fq_R) * 0.05)]
subsampled_path_R <- tempfile(fileext = ".fastq.gz")
writeFastq(subsampled_fq_R, subsampled_path_R, compress = TRUE)

# Learn Error Rates from Subsampled Data
errF <- learnErrors(subsampled_path_F, multithread = TRUE)
errR <- learnErrors(subsampled_path_R, multithread = TRUE)

# Load precomputed dada2 denoised objects to save time
dadaFs <- readRDS("dadaFs.rds")
dadaRs <- readRDS("dadaRs.rds")
# Or pick the file manually
dadaFs <- readRDS(file.choose())  # for forward reads
dadaRs <- readRDS(file.choose())  # for reverse reads

# Merge paired reads based on overlap
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Construct the ASV sequence table (samples as rows, ASVs as columns)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)  # Check dimensions of the table

# Load precomputed chimera-filtered table
seqtab.nochim <- readRDS("seqtab.nochim.rds")

# Assign taxonomy using a pre-trained SILVA classifier (can take a few minutes)
taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/PLS/16S/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)

# Load precomputed taxonomy assignments if preferred
taxa <- readRDS("taxa.rds")
taxa <- readRDS(file.choose())

# Inspect the results
taxa_print <- taxa
rownames(taxa_print) <- NULL
head(taxa_print)

# Read in sample metadata (CSV file with sample info)
sample_metadata <- read.csv("~/Documents/PLS/16S/metadata.csv", row.names=1)

# Create a phyloseq object from ASV table, taxonomy, and metadata
otu <- otu_table(seqtab.nochim, taxa_are_rows=FALSE)  # ASVs as columns
tax <- tax_table(taxa)
sam <- sample_data(sample_metadata)
ps <- phyloseq(otu, tax, sam)
ps  # View summary of the phyloseq object

# Alpha diversity: plot observed richness and Shannon diversity by Condition
plot_richness(ps, measures=c("Observed", "Shannon"), x="Condition", color="Condition")

# ---- Taxonomic Composition Plots ----

# Transform counts to relative abundance per sample
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Agglomerate ASVs at the Phylum level
ps.phylum <- tax_glom(ps.rel, taxrank = "Phylum")

# Convert phyloseq object to long-format dataframe (sample-taxon pairs)
ps.melt <- psmelt(ps.phylum)

# Average relative abundance by Condition and Phylum
avg_abundance <- ps.melt %>%
  group_by(Condition, Phylum) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

# Identify top 10 phyla
top_phyla <- avg_abundance %>%
  group_by(Phylum) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

# Group remaining phyla as "Other"
avg_abundance$Phylum <- as.character(avg_abundance$Phylum)
avg_abundance$Phylum[!(avg_abundance$Phylum %in% top_phyla)] <- "Other"
avg_abundance$Phylum <- factor(avg_abundance$Phylum, levels = c(sort(top_phyla), "Other"))

# Plot stacked bar chart of phylum composition by condition
ggplot(avg_abundance, aes(x = Condition, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Relative Abundance (%)") +
  xlab("Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        panel.grid = element_blank())

# Repeat for genus-level composition
ps.genus <- tax_glom(ps.rel, taxrank = "Genus")
ps.melt <- psmelt(ps.genus)
avg_abundance <- ps.melt %>%
  group_by(Condition, Genus) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

top_genera <- avg_abundance %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

avg_abundance$Genus <- as.character(avg_abundance$Genus)
avg_abundance$Genus[!(avg_abundance$Genus %in% top_genera)] <- "Other"
avg_abundance$Genus <- factor(avg_abundance$Genus, levels = c(sort(top_genera), "Other"))

ggplot(avg_abundance, aes(x = Condition, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Relative Abundance (%)") +
  xlab("Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        panel.grid = element_blank())



