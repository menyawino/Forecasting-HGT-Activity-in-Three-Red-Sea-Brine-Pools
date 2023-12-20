library(dada2)
library(DECIPHER)
library(phangorn)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(dplyr)


path <- "/home/menyawino/btt/Project/Datasets"
list.files(path)


# Fastq filenames
fnFs <- sort(list.files(path, pattern=".fastq", full.names = TRUE))
fnFs

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)
sample.names <- sapply(strsplit(sample.names, ".fa"), `[`, 1)
sample.names


# Inspect read quality profiles
initial_plot <- plotQualityProfile(fnFs[1:5])
ggsave(filename = "/home/menyawino/btt/Project/Datasets/plots/initialQ.png"
       , plot = initial_plot, dpi = 300)


# Assign the filenames for the filtered fastq.gz files
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names,
       "_filt.fastq.gz"))
names(filtFs) <- sample.names
filtFs


# Filter and trim
out <- filterAndTrim(fnFs, filtFs, truncLen=c(450), maxN=0, 
       maxEE=2, truncQ=2, rm.phix=TRUE, compress=TRUE, 
       multithread=TRUE)
head(out)


# Check Quality Profiles Again
trimQ_plot <- plotQualityProfile(filtFs[1:5])
ggsave(filename = "/home/menyawino/btt/Project/Datasets/plots/trimQ.png"
       , plot = trimQ_plot, dpi = 300)


# Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)


# visualize the estimated error rates
error_plot <- plotErrors(errF, nominalQ=TRUE) 
ggsave(filename = "/home/menyawino/btt/Project/Datasets/plots/error.png"
       , plot = error_plot, dpi = 300)


# dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
names(derepFs) <- sample.names
derepFs


# Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE,
       HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
dadaFs[[1]]


# Construct sequence table (ASV)
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
         multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab) #get the percent of non-chimeras


# Track the reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "nonchim")
rownames(track) <- sample.names
head(track)


# Assign taxonomy to the genus level
taxa <- assignTaxonomy(seqtab.nochim, 
                      "/home/menyawino/btt/practical/tax/RefSeq-RDP16S_v3_May2018.fa.gz", 
                      multithread=TRUE, tryRC=TRUE)


# Inspect the taxonomic assignments
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa.print, "/home/menyawino/btt/Project/ASV_taxonomy_RDP.csv")


# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}


# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "/home/menyawino/btt/Project/ASVs.fa")


# count table
asv_tax <- t(seqtab.nochim)
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "/home/menyawino/btt/Project/ASVs_counts.tsv", 
            sep="\t", quote=F, col.names=NA)


head(asv_tax)


asv_IDs <- taxa.print
rownames(asv_IDs) <- rownames(asv_tax)
asv_df <- as.table(asv_IDs)


write.csv(asv_df, file = "/home/menyawino/btt/Project/asv_class.csv")

##############################################################
###################   DADA2 ENDS HERE   ######################
###################      PHANGORN       ######################
##############################################################



# Align sequences 
# Extract sequences from DADA2 output
sequences <- getSequences(seqtab.nochim)
names(sequences) <- sequences

# Run seq. alignment using DECIPHER
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)

# Construct a phylogeny tree-
# change seq alignment output into a phyDat structure
phang.align <- phyDat(as(alignment,"matrix"), type ="DNA")

# create distance matrix
dm <-dist.ml(phang.align)

# perform neighbor joining 
treeNJ <- NJ(dm)

# internal max. likelihood
fit=pml(treeNJ,data=phang.align)

# negative edge length changed into 01 
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <-optim.pml(fitGTR, model= "GTR",optInv = TRUE, optGamma = TRUE
       ,rearrangement = "stochastic",control= pml.control(trace=0))


##############################################################
################### PHANGORN ENDS HERE  ######################
###################      PHYLOSEQ       ######################
##############################################################


# Import data into phyloseq, mapping file =metadata
# since it is not already prepared, we are going to prepare it
# We’ll also add the small amount of metadata we have 
# (the pool from which each sample was taken) to the phyloseq object.
samples.out <- rownames(seqtab.nochim)
section <- c("Atlantis II Deep S3", "Atlantis II Deep S5",
         "Chain Deep S1", "Discovery Deep S1", "Discovery Deep S3")
pool <- c("Atlantis II Deep", "Atlantis II Deep", "Chain Deep",
         "Discovery Deep", "Discovery Deep")
mapdf <- data.frame(Pool=pool, Section=section)
rownames(mapdf) <- samples.out
mapdf


# construct a phyloseq object directly from the dada2 outputs
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(mapdf), phy_tree(fitGTR$tree),
               tax_table(taxa))


estimate_richness(ps, split = TRUE, measures = NULL)
richness_plot <- plot_richness(ps, x="Section", measures=c("Shannon"), color="Pool")
ggsave(filename = "/home/menyawino/btt/Project/Datasets/plots/richness.png"
       , plot = richness_plot, dpi = 300)



# drop asv_df first six columns
asv_df_gs <- asv_df[,6:ncol(asv_df)]


# get species and genus name from asv_df and concatenate into strings in one vector
concatenated_strings <- apply(asv_df_gs, 1, function(row) paste(row, collapse = " "))


# 'concatenated_strings' contains the concatenated strings
head(concatenated_strings)


# add the asv name to the genus and species name in concatenated_strings
concatenated_strings <- paste0(taxa_names(ps), " ", concatenated_strings)
head(concatenated_strings)


# store the DNA sequences of our ASVs in the refseq slot of the 
# phyloseq object, and then rename our taxa to a short string.
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0(concatenated_strings)
ps


# Transform data to proportions for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
bray_plot <- plot_ordination(ps.prop, ord.nmds.bray, color="Pool", title="Bray NMDS")
ggsave(filename = "/home/menyawino/btt/Project/Datasets/plots/bray.png"
       , plot = bray_plot, dpi = 300)


# root phylogeny tree
set.seed(711)
phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps),1), resolve.root= TRUE)
is.rooted(phy_tree(ps))

# Plot the tree
tree_plot <- plot_tree(ps, color = 'Section', size = "abundance", min.abundance = Inf,
  label.tips = "taxa_names", text.size = NULL, sizebase = 5,
  base.spacing = 0.02, ladderize = FALSE, plot.margin = 0.2,
  title = 'Phylogenetic Tree', treetheme = NULL, justify = "jagged")
ggsave(filename = "/home/menyawino/btt/Project/Datasets/plots/tree.png"
       , plot = tree_plot, dpi = 300)


# Subset and plot trees for each sample
AtlantisII_S3 <- subset_samples(ps, Section == "Atlantis II Deep S3")
AtlantisII_S5 <- subset_samples(ps, Section == "Atlantis II Deep S5")
Chain_S1 <- subset_samples(ps, Section == "Chain Deep S1")
Discovery_S1 <- subset_samples(ps, Section == "Discovery Deep S1")
Discovery_S3 <- subset_samples(ps, Section == "Discovery Deep S3")

AtlantisII_S3_tree_plot <- plot_tree(Atlantis_II_S3, color = NULL , size = "abundance", 
  min.abundance = Inf,  label.tips = "taxa_names", text.size = NULL, sizebase = 5,
  base.spacing = 0.02, ladderize = FALSE, plot.margin = 0.2,
  title = 'Phylogenetic Tree', treetheme = NULL, justify = "jagged")

ggsave(filename = "/home/menyawino/btt/Project/Datasets/plots/AtlantisII_S3_tree.png"
       , plot = AtlantisII_S3_tree_plot, dpi = 300)

AtlantisII_S5_tree_plot <- plot_tree(Atlantis_II_S5, color = NULL , size = "abundance", 
  min.abundance = Inf,  label.tips = "taxa_names", text.size = NULL, sizebase = 5,
  base.spacing = 0.02, ladderize = FALSE, plot.margin = 0.2,
  title = 'Phylogenetic Tree', treetheme = NULL, justify = "jagged")

ggsave(filename = "/home/menyawino/btt/Project/Datasets/plots/AtlantisII_S5_tree.png"
       , plot = AtlantisII_S5_tree_plot, dpi = 300)

Chain_S1_tree_plot <- plot_tree(Chain_S1, color = NULL , size = "abundance",
       min.abundance = Inf,  label.tips = "taxa_names", text.size = NULL, sizebase = 5,
       base.spacing = 0.02, ladderize = FALSE, plot.margin = 0.2,
       title = 'Phylogenetic Tree', treetheme = NULL, justify = "jagged")

ggsave(filename = "/home/menyawino/btt/Project/Datasets/plots/Chain_S1_tree.png")
       , plot = Chain_S1_tree_plot, dpi = 300)

Discovery_S1_tree_plot <- plot_tree(Discovery_S1, color = NULL , size = "abundance",
       min.abundance = Inf,  label.tips = "taxa_names", text.size = NULL, sizebase = 5,
       base.spacing = 0.02, ladderize = FALSE, plot.margin = 0.2,
       title = 'Phylogenetic Tree', treetheme = NULL, justify = "jagged")

ggsave(filename = "/home/menyawino/btt/Project/Datasets/plots/Discovery_S1_tree.png"
       , plot = Discovery_S1_tree_plot, dpi = 300)

Discovery_S3_tree_plot <- plot_tree(Discovery_S3, color = NULL , size = "abundance",
       min.abundance = Inf,  label.tips = "taxa_names", text.size = NULL, sizebase = 5,
       base.spacing = 0.02, ladderize = FALSE, plot.margin = 0.2,
       title = 'Phylogenetic Tree', treetheme = NULL, justify = "jagged")

ggsave(filename = "/home/menyawino/btt/Project/Datasets/plots/Discovery_S3_tree.png"
       , plot = Discovery_S3_tree_plot, dpi = 300)


# ggtree plot

mergedps <- rarefy_even_depth(ps, rngseed=394582)
mergedps <- tax_glom(mergedps,"Order")
melt_simple <- psmelt(mergedps) %>%
               filter(Abundance < 120) %>%
               select(OTU, val=Abundance)
p <- ggtree(mergedps, layout="fan", open.angle=10) +
       geom_tippoint(mapping=aes(color=Phylum), 
                         size=1.5,
                         show.legend=FALSE)
p <- rotate_tree(p, -90)

p <- p +
       geom_fruit(
           data=melt_simple,
           geom=geom_boxplot,
           mapping = aes(
                       y=OTU,
                       x=val,
                       group=label,
                       fill=Phylum,
                     ),
           size=.2,
           outlier.size=0.5,
           outlier.stroke=0.08,
           outlier.shape=21,
           axis.params=list(
                           axis       = "x",
                           text.size  = 1.8,
                           hjust      = 1,
                           vjust      = 0.5,
                           nbreak     = 3,
                       ),
           grid.params=list()
       )

p <- p +
       scale_fill_discrete(
           name="Phyla",
           guide=guide_legend(keywidth=0.8, keyheight=0.8, ncol=1)
       ) +
       theme(
           legend.title=element_text(size=9), 
           legend.text=element_text(size=7) 
       )

ggsave(filename = "/home/menyawino/btt/Project/Datasets/plots/ggtree.png" 
       , plot = p, dpi = 300)



# Bar plot for the top 20 taxa
  
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]

ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))

ps.top20 <- prune_taxa(top20, ps.top20)
top20_plot <- plot_bar(ps.top20, x="Section", fill="Family")
ggsave(filename = "/home/menyawino/btt/Project/Datasets/plots/top20.png"
       , plot = top20_plot, dpi = 300)
       
make_biom(data, sample_metadata = NULL, observation_metadata = NULL,
id = NULL, matrix_element_type = "int")

b <-    make_biom(
        data = seqtab.nochim,
    )

write_biom(b, "/home/menyawino/btt/Project/asv_counts.biom")

save.image("/home/menyawino/btt/Project/dada2_final_project.RData")
