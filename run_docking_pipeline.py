import os
import sys
import subprocess
import glob
import shutil
from schrodinger import structure

def run_cmd(cmd):
    print(f"\nExecuting: {cmd}")
    result = subprocess.run(cmd, shell=True)
    return result.returncode == 0

def cleanup_results():
    print("\n>>> CLEANING UP PREVIOUS RESULTS...")
    folders = ["00_Significant_Poses", "01_Score_Analysis", "02_Interaction_Analysis", 
               "03_Clustering_Analysis", "04_Visualizations", "05_BPMD_Validation"]
    for folder in folders:
        if os.path.exists(folder): shutil.rmtree(folder)
    print("Cleanup complete.")

def main():
    schrodinger_run = os.environ.get('SCHRODINGER', '') + "/run"
    scripts_dir = "pipeline_scripts"
    
    if "--clean" in sys.argv:
        cleanup_results()
        sys.argv.remove("--clean")

    all_mae = glob.glob("DATA/*.mae") + glob.glob("DATA/*.maegz")
    input_files = [f for f in all_mae if "metadynamics" not in f.lower()]
    
    if not input_files: return
    
    # --- PHASE 1: COLLATING SIGNIFICANT POSES (RECEPTOR-SAFE) ---
    all_filtered_mae = "00_Significant_Poses/all_significant_poses.mae"
    os.makedirs("00_Significant_Poses", exist_ok=True)
    
    merged_ligands = []
    master_receptor = None
    
    for i, input_file in enumerate(input_files):
        print(f"\n>>> PROCESSING DOCKING FILE: {input_file}")
        # Run module0 but we'll manually handle the output to avoid receptor duplication
        run_cmd(f"{schrodinger_run} {scripts_dir}/module0_filter.py {input_file}")
        
        filtered_part = "00_Significant_Poses/filtered_significant_poses.mae"
        if os.path.exists(filtered_part):
            reader = list(structure.StructureReader(filtered_part))
            if not master_receptor: master_receptor = reader[0]
            # Structures 1+ are the significant ligand poses
            merged_ligands.extend(reader[1:])
            os.remove(filtered_part)

    if master_receptor:
        print(f"\n>>> SAVING {len(merged_ligands)} POSES TO MASTER FILE (Excluding extra receptors)")
        writer = structure.StructureWriter(all_filtered_mae)
        writer.append(master_receptor)
        for lig in merged_ligands:
            writer.append(lig)
        writer.close()
    
    # --- PHASE 2: INTEGRATED ANALYSIS ---
    if os.path.exists(all_filtered_mae):
        run_cmd(f"{schrodinger_run} {scripts_dir}/analyze_docking.py {all_filtered_mae}")
        run_cmd(f"{schrodinger_run} {scripts_dir}/module1_scores.py {all_filtered_mae}")
        run_cmd(f"{schrodinger_run} {scripts_dir}/module2_interactions.py {all_filtered_mae}")
        run_cmd(f"{schrodinger_run} {scripts_dir}/module2_hotspots.py")
        run_cmd(f"{schrodinger_run} {scripts_dir}/module3_clustering.py {all_filtered_mae}")
        run_cmd(f"{schrodinger_run} {scripts_dir}/module4_visualize.py")

    # --- PHASE 3: BPMD VALIDATION ---
    bpmd_files = glob.glob("DATA/*metadynamics*.maegz")
    if bpmd_files:
        run_cmd(f"{schrodinger_run} {scripts_dir}/module5_bpmd.py")
        run_cmd(f"{schrodinger_run} {scripts_dir}/module6_consensus.py")
        run_cmd(f"Rscript {scripts_dir}/mdplot.R")

    print("\nPIPELINE EXECUTION COMPLETE")

if __name__ == "__main__":
    main()
