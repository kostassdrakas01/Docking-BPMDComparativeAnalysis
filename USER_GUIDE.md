# User Guide: Metric Glossary & Interpretation

This guide provides detailed definitions for the terms, axes, and metrics used throughout the Docking-BPMD pipeline.

---

## 📊 Metric Glossary

### 1. Docking Metrics (Initial Screening)
*   **GlideScore**: An empirical value predicting binding affinity. **Lower (more negative) is better.** It accounts for H-bonds, lipophilic contacts, and steric penalties.
*   **Emodel**: Schrödinger’s primary metric for selecting the best binding pose *within* a set for a single ligand. It is heavily weighted by force-field interactions (vdW and Electrostatics).
*   **Ligand Efficiency (LE)**: Calculated as `GlideScore / Heavy Atom Count`. It measures how much each atom contributes to the binding. High LE indicates a highly optimized molecule.

### 2. BPMD Validation Metrics (Refinement)
*   **PoseScore (X-axis in Validation Landscape)**: Measures global geometric stability during force-biased metadynamics simulations.
    *   **> 1.0**: Stable fit (The ligand stayed in the pocket).
    *   **> 1.5**: **High Fidelity** (The orientation is extremely resilient).
*   **Persistence (Y-axis in Validation Landscape)**: Measures the life-time of specific chemical interactions (H-bonds, salt bridges).
    *   **> 1.0**: On average, at least one key interaction was maintained for the entire simulation across all trials.
*   **CompScore**: A "Composite Score" that mathematically combines initial docking energy with dynamic stability. **Lower is better.**

---

## 📈 Plot Interpretation

### Validation Landscape (`scientific_validation_landscape.png`)
*   **Green Zone (Top Right)**: High-fidelity binders. These poses are both physically stable and chemically consistent.
*   **Bottom Left**: Docking artifacts or transient binders that dissociate under force.

### Structural Anchor Profiles (`interaction_profile_*.png`)
*   **X-Axis**: Residue name (e.g., ASP125).
*   **Y-Axis**: Interaction Frequency (%).
*   **Stacking**: Shows the contribution of different interaction types (H-Bond, Salt-Bridge, Pi-Stacking).
*   **Interpretation**: Residues with high, stacked bars across many poses are your **Structural Anchors**.

---

## 🧬 Scientific Classifications

*   **High Fidelity Anchor** (`PoseScore > 1.5` AND `Persistence > 1.0`)
    *   **Meaning**: A bioactive candidate. The ligand is geometrically locked in place and its specific interaction network (H-bonds/Salt-bridges) is strong enough to resist metadynamics unbinding forces throughout the simulation.
*   **Stable Binding Mode** (`PoseScore > 1.0`)
    *   **Meaning**: A physically reasonable orientation. The ligand fits well in the pocket and does not clash, but it may lack the rigid, long-lived specific interactions required for High Fidelity.
*   **Metastable / Artifact** (`PoseScore < 1.0`)
    *   **Meaning**: Poses that drift significantly or dissociate early in the simulation. These are likely artifacts of the docking grid and should not be used for drug design decisions.
