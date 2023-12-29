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


<p float="left">
  <img src="https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/q_pretrim.png" width="49.5%" />
  <img src="https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/q_posttrim.png" width="49.5%" />
</p>

### Figure 1 – Quality profiles per sample before trimming 
graphs show a general quality drop after at location 450.

### Figure 2 – Quality profiles per sample after trimming
Reads were truncated at location 450 to remove low quality bases.

---


![Alt text](https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/err.png)
### Figure 3 – Estimated error rates per sample
Error rates for all samples show normal deviation from the optimal error.

---


| Sample       | Input  | Filtered | Denoised | Nonchim |
|--------------|-------:|---------:|---------:|--------:|
| AtlantisS3   | 16758  | 15262    | 14978    | 11833   |
| AtlantisS5   | 34243  | 26043    | 25529    | 21927   |
| ChainS1      | 21373  | 16313    | 13396    | 10811   |
| DiscoveryS1  | 25767  | 20279    | 15086    | 12530   |
| DiscoveryS7  | 26907  | 24059    | 23561    | 19218   |
### Table 1 – Number of Reads Through Processing
The processing steps show no significant loss of the data.

---


<p float="left">
  <img src="https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/shannon.png" width="49.5%" />
  <img src="https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/bray.png" width="49.5%" />
</p>

### Figure 4 – Alpha diversity measure using Shannon index.
Samples show difference in Shannon index.

### Figure 5 – Beta diversity measure using Bray-Curtis index.
Samples show difference in Bray-Curtis index.

---

![Alt text](https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/phylo.png)
### Figure 6 – Phylogeny Tree with Abundance
The tree shows qualitative diversity as many phyla are present in the samples and quantitative diversity as phyla abundance is varied.

---

<p float="left">
  <img src="https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/family_box.png" width="49.5%" />
  <img src="https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/genus_box.png" width="49.5%" />
</p>

### Figure 7 – Top Families Identified per Sample
The boxplot shows families dominating different samples.

### Figure 8 – Top Genus Identified per Sample
The boxplot shows Genus dominating different samples.

---

| Rank | Pathway                             | Rank | Pathway                             |
|------|-------------------------------------|------|-------------------------------------|
| 1    | Photosynthesis                      | 6    | Taurine and hypotaurine metabolism |
| 2    | Lipopolysaccharide biosynthesis     | 7    | Photosynthesis - antenna proteins  |
| 3    | Flagellar assembly                  | 8    | Carotenoid biosynthesis             |
| 4    | Tetracycline biosynthesis           | 9    | Arachidonic acid metabolism         |
| 5    | Ascorbate and aldarate metabolism   | 10   | Biofilm formation                   |

### Table 2 – Highest Abundant Pathways in Discovery Deep section 7
Highest active pathways after filtration and annotation

---

<img src="https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/pca.png" width="100%" />

### Figure 9 – Pathway Differential Abundance PCA
Samples show a clear distinction and can be separated by pool

---

<img src="https://github.com/menyawino/Forecasting-HGT-Activity-in-Three-Red-Sea-Brine-Pools/blob/main/results/heatmap.png" width="100%" />

### Figure 10 – Pathway Differential Abundance Heatmap
Different samples within the same pool show different pathway abundance

---

# Discussion
The quality plots shown in figure 1 show that the quality of the reads generally drop below a Phred score of 20  at location 450 in all datasets, hence why the value 450 was chosen as the trimming point for all the reads. Figure 2 shows the quality plots for the datasets after trimming, which depicts an overall improvement in the quality of the reads. Figure 3 shows a graphic representation of the error model produced based on the trimmed reads, which shows a normal deviation from the optimal error rate. Tracing the reads between processes showed no significant drop in the amount of data, which means the data is of high quality. The alpha diversity measure, shown in figure 4, shows a higher diversity in Atlantis II pool compared to Discovery Deep pool. Atlantis II pool scored 4.36, meaning that around 4 species dominate the sample, while Discovery Deep pool scored 4.12 which was not very different than Atlantis II. The beta diversity analysis with Bray-Curtis index in figure 5 was performed between the three pools showed differences between samples, with one sample from Discovery Deep close to Atlantis II samples in diversity.

The phylogeny tree in figure 6 shows the main phyla present in the samples. A qualitative as well as quantitative diversity are shown through the tree as there are various phyla in the samples with different abundances. Thaumarchaeota appears to be the most dominant across the tree, with the highest number of species and their respective abundances. The main reason behind Thaumarchaeota abundance is their ability to adapt to extreme environments. Red Sea brine pools are a highly saline and anoxic environment with high pressure and temperature. Thaumarchaeota are known to be the most abundant archaea on earth and have been found in various extreme environments. They display several adaptations to extreme environments, including the ability to survive high salinity and low ammonia conditions. Previous comparative genomics studies have revealed unique adaptations of Thaumarchaeota that allow them to survive harsh conditions, such as the presence of specific membrane lipids and conserved signature indels (Schouten et al., 2008). They are also capable of ammonia oxidation, which is a key metabolic feature that allows them to survive in environments with low ammonia (Zhang et al., 2020). This result is also consistent with other metagtenomic studies that showed that the abundance of Thaumarchaeota is particularly high in various Red Sea brine pools samples, indicating their specialization in these environments (Ziko et al., 2019). Euryarchaeota also appear in relatively high abundance in the tree. This can be due to several reasons, of which the most important is their osmoregulation, which is crucial to survive in highly saline environments (Castro‐Fernandez et al., 2017). Some Euryarchaeota are also methanogenic, they utilize unique enzymes and coenzymes to produce methane from carbon dioxide and hydrogen, which can serve as an energy source for other microorganisms (Huber, 2006).

The filtered box plots in figures 7 and 8 came consistent with the phylogeny tree for Atlantis II samples and Discovery section 7 sample which was consistent with the data obtained from Bray-Curtis diversity index in figure 6, yet Chain deep and Discovery Deep section 3 samples show very different distribution than the rest. The top family in the first mentioned three was Nitrosopumilaceae, which belongs to phylum Thaumarchaeota. Nitrosopumilaceae developed adaptations including the use of alternative electron acceptors to yield energy, such as iron, manganese, sulfate, elemental sulfur, carbon dioxide, nitrite, and nitrate (Zhang et al., 2020). The other family that showed noticeable abundance was Methanomassiliicoccales, which belongs to phylum Euryarchaeota. Species from this family can be found in various habitats, including hydrothermal vents and anoxic brine pools. They can tolerate up to 30% salinity levels and grow optimally in hypersaline conditions (around 15%) (Bueno de Mesquita et al., 2021). For Chain Deep and Discovery Deep section 1, SAR11 (Pelagibacterales) was the only abundant family. This was surprising as SAR11 are more adapted to oxic environments (Bougouffa et al., 2013). The genus in figure 8 shows the same trend of contrast between the samples which might be due to the spatial heterogeneity caused by high salinity levels in these brine pools.

To investigate HGT activity, PCA and heatmap were constructed to select the richest and most representative sample for pathway annotation. Figure 9 shows the difference between samples with principal components, the points were inseparable by pool, mostly due to the big difference in the metagenomic makeup of Discovery Deep section 1 sample. Section 7 sample was closer to Atlantis II samples which makes sense since they were similar in metagenomic composition. The heatmap in figure 10 also showed a difference in pathway abundance between all the samples, mostly in Discovery Deep section 1. The pathway analysis was done on Discovery Deep section 7 only since it conveyed the criteria mentioned above. The pathway analysis showed that HGT is not present in the sample, however, it showed other information on the pathways present. The highest abundant pathways in this environment were related to photosynthesis, Lipopolysaccharide biosynthesis, and flagellar assembly. The high abundance of photosynthesis-related pathways, such as photosystem I, II, and cytochrome b6/f complex, suggests that photosynthesis is a crucial process in this ecosystem despite its depth under the sea. Genes involved in lipopolysaccharide biosynthesis were highly abundant, which suggests that lipopolysaccharide production is an essential process in the microorganisms inhabiting this brine pool.



# Conclusion
While the main goal of this study was to explore the role of HGT in shaping microbial communities within the Red Sea brine pool environments, the direct detection of HGT associated pathways was not achieved. However, the comprehensive analysis yielded insights into the diverse microbial communities and functional characteristics present in Atlantis II, Discovery, and Chain Deep brine pools. The analysis showed the high dominance of Thaumarchaeota and Euryarchaeota in Atlantis II samples and one sample from Discovery Deep. The spatial heterogeneity observed within the brine pools resulted in the variability in microbial compositions in the two Discovery Deep samples. The prevalence of pathways related to photosynthesis and lipopolysaccharide biosynthesis highlighted essential biological processes within these extreme environments. With only five samples analyzed, this information serves as a first glance in the comparison between the differences in microbial communities and pathways present in the three brine pools analyzed. 
# Limitations
Of the primary limitations of this study is the low number of samples analyzed from the three brine pools. The small sample size often restricts the representation of the microbial communities in these complex environments. The findings might not represent the real landscape present in the analyzed pools. The study utilized PICRUSt2 R package to predict the pathways based on alignment of 16S rRNA data. While this approach offers insights into the functional potential, it relies on inference and prediction rather than direct genomic sequencing. This predictive nature could easily introduce inaccuracies which can potentially limit the validity of the predicted functional profiles. To ensure robustness and accuracy, the use of shotgun sequencing techniques would be indispensable for an accurate assessment of the functional potential.
# References

- Bougouffa, S., Yang, J. K., Lee, O. O., Wang, Y., Batang, Z., Al-Suwailem, A., & Qian, P. Y. (2013). Distinctive microbial community structure in highly stratified deep-sea brine water columns. Applied and Environmental Microbiology, 79(11), 3425–3437. [Link](https://doi.org/10.1128/aem.00254-13)

- Bueno de Mesquita, C. P., Zhou, J., Theroux, S. M., & Tringe, S. G. (2021). Methanogenesis and salt tolerance genes of a novel halophilic Methanosarcinaceae metagenome-assembled genome from a former Solar Saltern. Genes, 12(10), 1609. [Link](https://doi.org/10.3390/genes12101609)

- Elbehery, A. H. A., Beason, E., & Siam, R. (2023). Metagenomic profiling of antibiotic resistance genes in Red Sea brine pools. Archives of Microbiology, 205(5). [Link](https://doi.org/10.1007/s00203-023-03531-x)

- Harding, T., Roger, A. J., & Simpson, A. G. B. (2017). Adaptations to high salt in a halophilic protist: Differential expression and gene acquisitions through duplications and gene transfers. Frontiers in Microbiology, 8(MAY). [Link](https://doi.org/10.3389/fmicb.2017.00944)

- Huber, H. (2006). Euryarchaeota. Encyclopedia of Life Sciences. [Link](https://doi.org/10.1038/npg.els.0004243)

- López-Aladid, R., Fernández-Barat, L., Alcaraz-Serrano, V., Bueno-Freire, L., Vázquez, N., Pastor-Ibáñez, R., Palomeque, A., Oscanoa, P., & Torres, A. (2023). Determining the most accurate 16S rrna hypervariable region for taxonomic identification from respiratory samples. Scientific Reports, 13(1). [Link](https://doi.org/10.1038/s41598-023-30764-z)

- McMURDIE, P. J., & HOLMES, S. (2011). PHYLOSEQ: A bioconductor package for handling and analysis of high-throughput phylogenetic sequence data. Biocomputing 2012. [Link](https://doi.org/10.1142/9789814366496_0023)

- Ochman, H., Lawrence2, J. G., & Groisman3, E. A. (2000). Lateral gene transfer and the nature of bacterial innovation. In NATURE (Vol. 405). [Link](www.nature.com)

- Ombrello, A. K. (2020). DADA2. Encyclopedia of Medical Immunology, 251–257. [Link](https://doi.org/10.1007/978-1-4614-8678-7_118)

- Parkhill, J., Sebaihia, M., Preston, A., Murphy, L. D., Thomson, N., Harris, D. E., Holden, M. T. G., Churcher, C. M., Bentley, S. D., Mungall, K. L., Cerdeño-Tárraga, A. M., Temple, L., James, K., Harris, B., Quail, M. A., Achtman, M., Atkin, R., Baker, S., Basham, D., … Maskell, D. J. (2003). Comparative analysis of the genome sequences of bordetella pertussis, bordetella parapertussis and bordetella bronchiseptica. Nature Genetics, 35(1), 32–40. [Link](https://doi.org/10.1038/ng1227)

- Qian, P. Y., Wang, Y., Lee, O. O., Lau, S. C. K., Yang, J., Lafi, F. F., Al-Suwailem, A., & Wong, T. Y. H. (2011). Vertical stratification of microbial communities in the Red Sea revealed by 16S rDNA pyrosequencing. ISME Journal, 5(3), 507–518. [Link](https://doi.org/10.1038/ismej.2010.112)

- Schliep, K. P. (2010). Phangorn: Phylogenetic analysis in R. Bioinformatics, 27(4), 592–593. [Link](https://doi.org/10.1093/bioinformatics/btq706)

- Schouten, S., Hopmans, E. C., Baas, M., Boumann, H., Standfest, S., Könneke, M., Stahl, D. A., & Sinninghe Damsté, J. S. (2008). Intact membrane lipids of candidatus Nitrosopumilus Maritimus,” a cultivated representative of the cosmopolitan Mesophilic Group I Crenarchaeota. Applied and Environmental Microbiology, 74(8), 2433–2440. [Link](https://doi.org/10.1128/aem.01709-07)

- Tully, B. (2016). Quality Assessment: FASTQC V1. [Link](https://doi.org/10.17504/protocols.io.fa3bign)

- Yang, C., Mai, J., Cao, X., Burberry, A., Cominelli, F., & Zhang, L. (2023). GGPICRUST2: An R package for picrust2 predicted functional profile analysis and visualization. Bioinformatics, 39(8). [Link](https://doi.org/10.1093/bioinformatics/btad470)

- Ziko, L., Adel, M., Malash, M. N., & Siam, R. (2019). Insights into Red Sea brine pool specialized metabolism gene clusters encoding potential metabolites for biotechnological applications and Extremophile Survival. Marine Drugs, 17(5), 273. [Link](https://doi.org/10.3390/md17050273)

- Zhang, X.-H., Zhong, H., Lehtovirta-Morley, L., Liu, J., Zheng, Y., Lin, H., Song, D., Todd, J. D., & Tian, J. (2020). Novel Insights into the Thaumarchaeota in the Deepest Oceans: Their Metabolism and Potential Adaptation Mechanisms. [Link](https://doi.org/10.21203/rs.3.rs-162



