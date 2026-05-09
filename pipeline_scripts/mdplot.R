# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

# Set directories
bpmd_file <- "05_BPMD_Validation/bpmd_validation_results.csv"
final_file <- "05_BPMD_Validation/final_validated_interactions.csv"
int_file <- "02_Interaction_Analysis/interaction_frequency_summary.csv"
output_dir <- "04_Visualizations"
dir.create(output_dir, showWarnings = FALSE)

publication_theme <- theme_bw() +
  theme(panel.grid.minor = element_blank(), text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

# --- DATA PREP: Merging Glide and BPMD ---
if (file.exists(int_file) && file.exists(final_file)) {
  df_int <- read.csv(int_file, check.names = FALSE)
  df_final <- read.csv(final_file)
  
  df_val <- df_final %>%
    separate_rows(Validated_Interactions, sep = ",") %>%
    mutate(Residue_Val = trimws(Validated_Interactions)) %>%
    filter(Residue_Val != "Fingerprint Missing" & Residue_Val != "No specific interactions detected") %>%
    mutate(MatchKey = gsub("[()]", "", Residue_Val)) %>%
    mutate(MatchKey = gsub(" ", "_", MatchKey)) %>%
    mutate(BPMD_Confirmed = TRUE) %>%
    select(Ligand, MatchKey, BPMD_Confirmed)

  df_int$Freq <- df_int[[grep("Frequency", names(df_int), value = TRUE)[1]]]
  df_int$MatchKey <- paste(df_int$Residue, df_int$Interaction_Type, sep="_")
  
  df_merged <- df_int %>%
    left_join(df_val, by = c("Ligand", "MatchKey")) %>%
    mutate(BPMD_Confirmed = ifelse(is.na(BPMD_Confirmed), FALSE, TRUE))

  # --- PLOT 1: Combined Integrated Validation (Faceted) ---
  df_top_merged <- df_merged %>% 
    group_by(Ligand) %>% 
    slice_max(order_by = Freq, n = 8, with_ties = FALSE)
  
  p_combined <- ggplot(df_top_merged, aes(x = reorder(Residue, -Freq), y = Freq)) +
    geom_bar(aes(fill = "Docking Noise (Glide)"), stat = "identity", position = "stack", alpha = 0.2) +
    geom_bar(data = filter(df_top_merged, BPMD_Confirmed == TRUE), 
             aes(fill = Interaction_Type), stat = "identity", position = "stack") +
    facet_wrap(~Ligand, scales = "free_x") +
    publication_theme +
    scale_fill_manual(values = c("Docking Noise (Glide)" = "grey70", 
                                "H-Bond" = "#E41A1C", "Salt-Bridge" = "#377EB8", 
                                "Pi-Stacking" = "#4DAF4A", "Pi-Cation" = "#984EA3",
                                "Contact" = "#FF7F00", "Pi-Stacking (Face-to-Face)" = "#4DAF4A",
                                "Pi-Stacking (Edge-to-Face)" = "#4DAF4A")) +
    labs(title = "Integrated Validation Profiles: All Ligands",
         subtitle = "Comparison of Glide Predictions (Grey) vs. BPMD Validated Anchors (Color)",
         x = "Residue", y = "Frequency (%)", fill = "Status")
  
  ggsave(file.path(output_dir, "integrated_validation_combined.png"), p_combined, width = 16, height = 10)

  # [Individual plots logic remains same...]
  for (lig in unique(df_merged$Ligand)) {
    df_lig <- df_merged %>% filter(Ligand == lig) %>% slice_max(order_by = Freq, n = 10, with_ties = FALSE)
    p_val <- ggplot(df_lig, aes(x = reorder(Residue, -Freq), y = Freq)) +
      geom_bar(aes(fill = "Docking Noise (Glide)"), stat = "identity", position = "stack", alpha = 0.2) +
      geom_bar(data = filter(df_lig, BPMD_Confirmed == TRUE), 
               aes(fill = Interaction_Type), stat = "identity", position = "stack") +
      publication_theme +
      scale_fill_manual(values = c("Docking Noise (Glide)" = "grey70", 
                                  "H-Bond" = "#E41A1C", "Salt-Bridge" = "#377EB8", 
                                  "Pi-Stacking" = "#4DAF4A", "Pi-Cation" = "#984EA3",
                                  "Contact" = "#FF7F00", "Pi-Stacking (Face-to-Face)" = "#4DAF4A",
                                  "Pi-Stacking (Edge-to-Face)" = "#4DAF4A")) +
      labs(title = paste("Integrated Validation Profile:", lig), x = "Residue", y = "Frequency (%)")
    safe_name <- gsub("[^[:alnum:]]", "_", lig)
    ggsave(file.path(output_dir, paste0("integrated_validation_", safe_name, ".png")), p_val, width = 9, height = 6)
  }
}

# --- PLOT 2, 3, 4: Consensus, Landscape, Refinement ---
if (file.exists(final_file)) {
  df_final <- read.csv(final_file)
  df_final$Short_Label <- sapply(strsplit(as.character(df_final$Validated_Interactions), ","), function(x) paste(trimws(head(x, 2)), collapse="\n"))
  p_con <- ggplot(df_final, aes(x = reorder(Ligand, -BPMD_PoseScore), y = BPMD_PoseScore, fill = BPMD_Persistence)) +
    geom_bar(stat = "identity", width = 0.7) + geom_text(aes(label = Short_Label), vjust = -0.5, size = 3, lineheight = 0.8) +
    scale_fill_gradient(low = "yellow", high = "darkgreen") + publication_theme +
    labs(title = "Final Pipeline Consensus Summary", x = "Ligand", y = "BPMD PoseScore") + ylim(0, max(df_final$BPMD_PoseScore) * 1.3)
  ggsave(file.path(output_dir, "final_consensus_summary.png"), p_con, width = 12, height = 7)
}
if (file.exists(bpmd_file)) {
  df_bpmd <- read.csv(bpmd_file)
  p_land <- ggplot(df_bpmd, aes(x = PoseScore, y = Persistence, color = Ligand)) +
    annotate("rect", xmin = 1.5, xmax = Inf, ymin = 1.0, ymax = Inf, alpha = 0.1, fill = "green") +
    geom_point(size = 5, alpha = 0.7) + publication_theme + labs(title = "BPMD Stability Landscape")
  ggsave(file.path(output_dir, "scientific_validation_landscape.png"), p_land, width = 8, height = 6)
}

message("Combined Integrated Validation Plot generated.")
