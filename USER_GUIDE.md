# Integrated Docking & BPMD Analysis Pipeline

This pipeline automates the post-processing of Schrödinger Glide docking results and validates them using Desmond Binding Pose Metadynamics (BPMD) scores. It provides a structured, theory-driven workflow to identify high-fidelity binding modes.

## 🚀 Quick Start

To run the entire analysis on your data, use the master command:

```bash
$SCHRODINGER/run run_docking_pipeline.py DATA/your_docking_results.mae
```

*Note: If no file is provided, the script will automatically detect any `.mae` file in the `DATA/` folder.*

---

## 📂 Project Organization

Results are automatically organized into numbered directories to maintain a clean, logical workflow:

### 00 | Pose Filtering (Significance Filter)
- **Theory**: Uses **Emodel** to select the top 10 most physically stable poses per ligand.
- **Output**: `filtered_significant_poses.mae`. This file is the input for all downstream steps.

### 01 | Score & Summary Analysis
- **Contents**: Basic docking statistics, GlideScore distributions, and Ligand Efficiency (LE).
- **Key File**: `score_distribution.csv` (Sorted by ligand name).

### 02 | Interaction Fingerprinting
- **Contents**: Residue-level contact maps (H-Bonds, Salt-Bridges, Pi-Stacking).
- **Key File**: `key_hotspot_residues.csv`. Identifies "Structural Anchors" based on consensus frequency.

### 03 | Clustering Analysis
- **Contents**: Groups poses by RMSD (2.0 Å threshold) using Schrödinger's library.
- **Purpose**: Distinguishes between unique binding orientations.

### 04 | Visualizations (Scientific Reports)
- **Heatmaps**: Detailed residue interaction patterns for each ligand.
- **Validation Landscape**: Scatter plot correlating Stability (BPMD) vs. Persistence.
- **Docking Refinement**: Comparison of initial GlideScore vs. refined CompScore.
- **Stability Summary**: Bar chart of validated "Structural Anchors."

### 05 | BPMD Validation (The Final Answer)
- **Contents**: Metadynamics stability metrics (PoseScore, Persistence).
- **Key File**: `final_validated_interactions.csv`. **THIS IS YOUR FINAL ANSWER.** It lists the residues that are "for sure" interacting in the most stable dynamic pose.

---

## 🧬 Scientific Interpretation

### How to Validate a "Hit"
1. Check **`04_Visualizations/scientific_validation_landscape.png`**.
2. Look for ligands in the **Green "High Fidelity" Zone** (PoseScore > 1.5, Persistence > 1.0).
3. These ligands have binding modes that are physically stable and reproducible.

### Resolving "Fuzzy" Docking
If docking shows many different clusters (e.g., **Analog1**), refer to **`05_BPMD_Validation/final_validated_interactions.csv`**. The BPMD simulation acts as a "tie-breaker" to identify which of those many orientations is actually stable over time.

---

## 📊 Metric Analysis & Glossary

To help you interpret the plots and CSV results, here is a breakdown of every metric used in the pipeline:

### 1. Docking Metrics (Module 01)
*   **GlideScore (X-axis in Affinity plots)**: An empirical value predicting binding affinity. **Lower (more negative) is better.** It accounts for H-bonds, lipophilic contacts, and penalties.
*   **Emodel Rank (Y-axis in Heatmaps)**: Schrödinger’s primary metric for choosing the best pose for a single ligand. It is heavily weighted by force-field interactions (vdW and Electrostatics).
*   **Ligand Efficiency (LE)**: Calculated as `GlideScore / Heavy Atom Count`. It measures how much each atom contributes to the binding. High LE means a small molecule that "punches above its weight."

### 2. BPMD Metadynamics Metrics (Module 05)
*   **PoseScore (X-axis in Validation Landscape)**: Measures how well the ligand stayed in the starting orientation during force-based simulations. 
    *   **Value > 1.0**: Physically stable.
    *   **Value > 1.5**: Extremely stable (High Fidelity).
*   **Persistence (Y-axis in Validation Landscape)**: Measures the "lifetime" of specific H-bonds and Salt-Bridges.
    *   **Value > 1.0**: The interaction was maintained for the entire simulation on average.
*   **CompScore**: A "Composite Score" that combines Docking and BPMD results. **Lower is generally more confident.**

### 3. Binding Mode Classifications
*   **High Fidelity Anchor**: A "bioactive" candidate. It is both stable (High PoseScore) and has persistent specific interactions (High Persistence).
*   **Stable Binding Mode**: A "reasonable" candidate. It fits well in the pocket but might not have a single dominant, permanent interaction.
*   **Metastable / Artifact**: Poses that dissociate or drift significantly. These should be excluded from final structural interpretations.

---

## 🛠 Maintenance
- **Scripts**: All logic is stored in `pipeline_scripts/`.
- **Adding Data**: Place new `.mae` files or BPMD `-out.maegz` files in `DATA/` and re-run the pipeline.
