# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Set directories
bpmd_file <- "05_BPMD_Validation/bpmd_validation_results.csv"
int_file <- "02_Interaction_Analysis/interaction_frequency_summary.csv"
output_dir <- "04_Visualizations"
dir.create(output_dir, showWarnings = FALSE)

publication_theme <- theme_bw() +
  theme(panel.grid.minor = element_blank(), text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

# --- PLOT 1: Interaction Profiles (Individual & Combined) ---
if (file.exists(int_file)) {
  df_int <- read.csv(int_file, check.names = FALSE)
  freq_col <- grep("Frequency", names(df_int), value = TRUE)
  
  if (length(freq_col) > 0) {
    df_int$Freq <- df_int[[freq_col[1]]]
    
    # A. Combined Faceted Plot
    df_top_all <- df_int %>% group_by(Ligand) %>% slice_max(order_by = Freq, n = 6, with_ties = FALSE)
    p_all <- ggplot(df_top_all, aes(x = reorder(Residue, -Freq), y = Freq, fill = Interaction_Type)) +
      geom_bar(stat = "identity", position = "stack") + facet_wrap(~Ligand, scales = "free_x") +
      publication_theme + labs(title = "Combined Interaction Profiles", x = "Residue", y = "Frequency (%)") +
      scale_fill_brewer(palette = "Set2")
    ggsave(file.path(output_dir, "ligand_interaction_profiles_combined.png"), p_all, width = 14, height = 8)
    
    # B. Individual Plots for each Ligand (Replacing Heatmaps)
    for (lig in unique(df_int$Ligand)) {
      df_lig <- df_int %>% filter(Ligand == lig) %>% slice_max(order_by = Freq, n = 10, with_ties = FALSE)
      
      p_lig <- ggplot(df_lig, aes(x = reorder(Residue, -Freq), y = Freq, fill = Interaction_Type)) +
        geom_bar(stat = "identity", position = "stack") +
        publication_theme +
        labs(title = paste("Interaction Profile:", lig),
             subtitle = "Top-10 stabilized residues identified across docking ensemble",
             x = "Residue", y = "Interaction Frequency (%)", fill = "Type") +
        scale_fill_brewer(palette = "Set2")
      
      # Use a safe filename
      safe_name <- gsub("[^[:alnum:]]", "_", lig)
      ggsave(file.path(output_dir, paste0("interaction_profile_", safe_name, ".png")), p_lig, width = 8, height = 6)
    }
  }
}

# --- PLOT 2: Validation Landscape ---
if (file.exists(bpmd_file)) {
  df_bpmd <- read.csv(bpmd_file)
  df_bpmd$Status <- ifelse(df_bpmd$PoseScore > 1.5 & df_bpmd$Persistence > 1.0, "High Fidelity", "Stable/Other")
  p2 <- ggplot(df_bpmd, aes(x = PoseScore, y = Persistence, color = Ligand)) +
    annotate("rect", xmin = 1.5, xmax = Inf, ymin = 1.0, ymax = Inf, alpha = 0.1, fill = "green") +
    geom_point(aes(shape = Status), size = 5) + publication_theme +
    labs(title = "BPMD Stability Landscape", x = "PoseScore", y = "Persistence")
  ggsave(file.path(output_dir, "scientific_validation_landscape.png"), p2, width = 8, height = 6)
}

# --- PLOT 3: Docking Refinement ---
if (file.exists(bpmd_file)) {
  p3 <- ggplot(df_bpmd, aes(x = reorder(Full_Title, -CompScore))) +
    geom_segment(aes(xend = reorder(Full_Title, -CompScore), y = GlideScore, yend = CompScore), color = "grey") +
    geom_point(aes(y = GlideScore, color = "Glide"), size = 3) +
    geom_point(aes(y = CompScore, color = "BPMD"), size = 4) +
    coord_flip() + publication_theme +
    labs(title = "Energy Refinement: Glide vs. BPMD", x = "Pose", y = "Score")
  ggsave(file.path(output_dir, "docking_refinement_comparison.png"), p3, width = 10, height = 7)
}

message("Individual interaction profiles generated (Heatmaps replaced).")
