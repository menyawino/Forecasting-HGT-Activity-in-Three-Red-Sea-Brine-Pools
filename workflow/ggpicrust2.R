library(ggpicrust2)

custom_colors <- c("red", "blue")


# Load KEGG pathway abundance
metadata <- read_delim("D:/metadata2.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
kegg_abundance <- ko2kegg_abundance("D:/pred_metagenome_unstrat.tsv") 


# Perform pathway differential abundance analysis (DAA) using DESeq2 method
daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, 
                              group = "Pool", daa_method = "DESeq2", select = NULL, reference = NULL)


# Filter results for DESeq2 test method
daa_results_df <- na.omit(daa_results_df)


# Annotate pathway results using KO to KEGG conversion
daa_annotated_results_df <- pathway_annotation(pathway = "KO", 
                              daa_results_df = daa_results_df, ko_to_kegg = TRUE)

daa_annotated_results_df <- daa_annotated_results_df[!is.na(daa_annotated_results_df$pathway_name),]
daa_annotated_results_df$p_adjust <- round(daa_annotated_results_df$p_adjust,5)
low_p_feature <- daa_annotated_results_df[order(daa_annotated_results_df$p_adjust), ]

# Generate pathway error bar plot
pathway_errorbar(abundance = kegg_abundance,
                      daa_results_df = daa_annotated_results_df,
                      Group = metadata$Pool,
                      ko_to_kegg = TRUE,
                      p_values_threshold = 0.05,
                      order = "pathway_class",
                      select = low_p_feature,
                      p_value_bar = TRUE,
                      colors = NULL,
                      x_lab = "pathway_name")


kegg_abundance$feature <- rownames(kegg_abundance)
joined_df <- left_join(low_p_feature, kegg_abundance, by = "feature")
selected_columns <- c("AtlantisS3", "AtlantisS5", "DiscoveryS1", "DiscoveryS7", "pathway_name")
joined_df_clean <- joined_df[, selected_columns, drop = FALSE]
joined_df_clean <- na.omit(joined_df_clean)
rownames(joined_df_clean) <- joined_df_clean[, "pathway_name"]
joined_df_clean$pathway_name <- NULL


# Generate pathway heatmap
pathway_heatmap(abundance = joined_df_clean, metadata = metadata, group = "Pool", colors = custom_colors)


# Generate pathway PCA
pathway_pca(abundance = joined_df_clean, metadata = metadata, "Pool")
