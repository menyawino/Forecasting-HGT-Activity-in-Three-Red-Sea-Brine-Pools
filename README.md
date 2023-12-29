# Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools

# Introduction
Microbial communities have taught us a lot about the dynamics of life, including lessons about ourselves. Many of today's clinical and industrial applications are inspired by studying the interaction between microorganisms. Typically, the most unique sources to discover microbiomes are in extreme environments, such as hydrothermal vents, hypersaline lakes, and brine pools. Over the past decade, there has been a growing research activity to discover the Red Sea brine pools. Brine pools such as Atlantis II, Discovery Deep, and Chain Deep showed interesting features that enabled further research to extract anti-cancer and anti-microbial proteins and chemicals and understand the effects of shared phenotypes on microbial communities. Metagenomics studies are the primary tool that enables the investigation of the diversity of life in such environments. In this study, we aim to investigate the diversity of Horizontal Gene Transfer (HGT) systems in three brine pools: Atlantis II, Discovery Deep, and Chain Deep. The primary objective of this literature review is to identify previously published research that can provide relevant input for the project. In this report, we conduct a thorough review of the literature, identify potential publicly available datasets, propose a research question, and finally outline an associated experimental design. It is important to note that the project will rely exclusively on publicly available data, and no new sequencing will be conducted.
## Microbial Adaptation and Horizontal Gene Transfer
Microbes are microscopic organisms that include bacteria and archaea (prokaryotes), viruses, fungi, and protists. These organisms, especially prokaryotes and viruses, are known to rapidly adapt to novel biotic and abiotic environmental changes by rapidly altering their genomes, resulting in beneficial phenotypes. Understanding the genetic bases for these adaptations was made possible by studying the genomes of isolated organisms of interest (Harding et al., 2017; Parkhill et al., 2003). New genes can be acquired through three processes: gene duplication, gain of a novel gene (e.g., in previously noncoding DNA), gene loss, and horizontal gene transfer (HGT). HGT, or lateral gene transfer, transfers genetic material outside parent–offspring inheritance. Unlike the other processes mentioned above, HGT enables the rapid transfer of genes between distantly related species and is thought to be especially important for adapting microbes to novel environments (Ochman et al., 2000). 
### Red Sea Brine Pools
Due to their unique characteristics, Red Sea brine pools have become a hotspot for metagenomic studies. Characterized by extreme salinity and high temperatures, these environments house microorganisms with exceptional adaptation mechanisms. Microbial communities' density and composition are formed through variations in temperature, salinity, nutrient availability, and dissolved oxygen (DO) levels across different depths. Recent metagenomic analyses of the Discovery Deep and Atlantis II brine pools revealed intriguing insights. Although similar microbial communities were identified in the surface waters (20 - 50 m), clearly distinct communities were observed in the deeper columns (200 m - 1500 m). This observation suggests a clear vertical stratification of microbial communities depending on depth (Qian et al., 2011). Interestingly, similar depths at different locations exhibited different microbial profiles. Upper layers featured higher microbial density and diversity, likely due to factors such as increased light, temperature, DO, and primary productivity. Cyanobacteria dominated the bacterial population in the surface waters, whereas Proteobacteria, particularly the gamma subdivision, became prevalent in the deeper layers (200 m -1500 m). Halophilic Archaea Haplobacteriales dominated the upper layers, while the deeper layers exhibited an enrichment of anaerobic sulfur-dependent heterotrophic Archaea Desulfurococcales. The presence of Desulfurococcales in the deeper layers can be attributed to the high sulfur concentration (Qian et al., 2011) and lower DO levels, where these Archaea likely play a crucial role in sulfur reduction and cycling.
### Datasets
The data used in this study is adapted from the samples collected at various sections within the Atlantis II, Discovery, and Chain Deep brine pools in the Red Sea. The samples were sequenced using the 454 GS FLX Titanium instrument platform at The American University in Cairo (AUC). The amplicons used in these samples are V4-V6, which was shown to be the best option for precise microorganism identification (López-Aladid et al., 2023). Each dataset represents a specific section of the brine pool sediments, allowing for an in-depth exploration of the diverse prokaryotic consortia thriving at varying depths and locations within these extreme environments.
### Research Question 
In this study, we aim to explore the microbial communities thriving in the extreme environments of Atlantis II, Discovery, and Chain Deep brine pools. Our research focuses on characterizing the microbial composition within these unique ecosystems and understanding HGT effects in brine pools. Metagenomics analysis will be performed to profile the microbial diversity and abundance in the brine pools. This analysis will serve as the foundation for our subsequent investigations. In parallel, we will investigate whether HGT is crucial in enabling functional independence from microbial community compositions. We will examine the conservation of functions across different microbial communities to evaluate the proposed hypothesis.

# Methods
To investigate the microbial diversity among the three brine pools, 5 samples were fetched: two sections from Atlantis II, two sections from Discovery Deep, and one section from Chain Deep. Metagenomic analyses were done separately on each sample, then the results were combined to draw conclusions from. This workflow and the datasets used are available at https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools.git 
### Quality Control
Quality control was done using FastQC, MultiQC, and the DADA2 metagenomics package quality plot functionality (Tully, 2016; Ewels et al., 2016; Ombrello, 2020). More than one tool were used to give access to largest amount of data on the samples. The visualized data allowed us to proceed with read filtering and trimming.
### Constructing ASV Abundance Table and Species Identification
To construct ASV abundance table, DADA2 R package was used (Ombrello, 2020). First of all, trimming was performed on the reads based on the quality analysis results in the previous step and the actual length of the amplicons used (V4-V6). The truncation length was set to 450 bases, with read filtration if more than two expected errors within the same read. There was no need to remove Phix genomic reads as this step is only applicable in Illumina sequencing techniques. The filtered samples were then dereplicated to remove duplicate sequences Error rate model was then constructed for the samples and sanity check error plots were produced. The error model was used, along with the filtered reads to infer and denoise ASVs. ASV table was constructed and chimeras were removed as well. The final ASVs were then aligned to RDP dataset (2018 version) to assign species to the genus or species level. 
### Constructing Phylogeny Trees
To construct phylogeny trees for the resulting ASVs, multiple sequence alignment was performed using Phangorn R package (Schliep, 2010). The alignment was used to calculate the distance between species and construct the tree. Phyloseq R package was used to perform estimates on the alpha diversity using Shannon index and beta diversity using bray-curtis index and plot the phylogeny tree (McMURDIE & HOLMES, 2011). ggtree R package was used to perform further phylogeny tree plotting. Finally, the top 20 species in terms of abundance were selected to investigate the metagenomic diversity between the five samples. 
### Predicting Functional Pathway Abundance
To analyze the Horizontal Gene Transfer (HGT) activity in the three brine pools selected, PICRUSt 2.0 R package was used first to align the 16s data to the reference trees: EPA-NG and GAPPA (Douglas et al., 2020). The alignment was used to infer gene family copy number of ASVs for each sample independently. The predicted sample gene family profiles were mapped to their pathways and pathway abundances for each sample were calculated. The pathway abundance table was then passed to ggpicrust2 R package for further analysis and visualization (Yang et al., 2023). Analysis of the differential pathway abundance was performed on Atlantis II and Discovery Deep pools only as computational power was restricting differential abundance analysis for more than two groups. 

# Results

![Alt text](https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/q_pretrim.png)
### Figure 1 – Quality profiles per sample before trimming 
graphs show a general quality drop after at location 450.


![Alt text](https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/q_posttrim.png)
### Figure 2 – Quality profiles per sample after trimming
Reads were truncated at location 450 to remove low quality bases.


![Alt text](https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/err.png)
### Figure 3 – Estimated error rates per sample
Error rates for all samples show normal deviation from the optimal error.


| Sample       | Input  | Filtered | Denoised | Nonchim |
|--------------|-------:|---------:|---------:|--------:|
| AtlantisS3   | 16758  | 15262    | 14978    | 11833   |
| AtlantisS5   | 34243  | 26043    | 25529    | 21927   |
| ChainS1      | 21373  | 16313    | 13396    | 10811   |
| DiscoveryS1  | 25767  | 20279    | 15086    | 12530   |
| DiscoveryS7  | 26907  | 24059    | 23561    | 19218   |
### Table 1 – Number of Reads Through Processing
The processing steps show no significant loss of the data.


![Alt text](https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/shannon.png)
### Figure 4 – Alpha diversity measure using Shannon index.
Samples show difference in Shannon index.


![Alt text](https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/bray.png)
### Figure 5 – Beta diversity measure using Bray-Curtis index.
Samples show difference in Bray-Curtis index.


![Alt text](https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/phylo.png)
Figure 6 – Phylogeny Tree with Abundance
The tree shows qualitative diversity as many phyla are present in the samples and quantitative diversity as phyla abundance is varied.


![Alt text](https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/family_box.png)
### Figure 7 – Top Families Identified per Sample
The boxplot shows families dominating different samples.


![Alt text](https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/genus_box.png)
### Figure 8 – Top Genus Identified per Sample
The boxplot shows Genus dominating different samples.


| Rank | Pathway                             | Rank | Pathway                             |
|------|-------------------------------------|------|-------------------------------------|
| 1    | Photosynthesis                      | 6    | Taurine and hypotaurine metabolism |
| 2    | Lipopolysaccharide biosynthesis     | 7    | Photosynthesis - antenna proteins  |
| 3    | Flagellar assembly                  | 8    | Carotenoid biosynthesis             |
| 4    | Tetracycline biosynthesis           | 9    | Arachidonic acid metabolism         |
| 5    | Ascorbate and aldarate metabolism   | 10   | Biofilm formation                   |

### Table 2 – Highest Abundant Pathways in Discovery Deep section 7
Highest active pathways after filtration and annotation


![Alt text](https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/pca.png)
### Figure 9 – Pathway Differential Abundance PCA
Samples show a clear distinction and can be separated by pool


![Alt text](https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/heatmap.png)
Figure 10 – Pathway Differential Abundance Heatmap
Different samples within the same pool show different pathway abundance

