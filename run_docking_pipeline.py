import os
import sys
import subprocess
import glob
import shutil

def run_cmd(cmd):
    print(f"\nExecuting: {cmd}")
    result = subprocess.run(cmd, shell=True)
    return result.returncode == 0

def cleanup_results():
    print("\n>>> CLEANING UP PREVIOUS RESULTS...")
    folders = [
        "00_Significant_Poses", "00_Initial_Docking_Summary", 
        "01_Score_Analysis", "02_Interaction_Analysis", 
        "03_Clustering_Analysis", "04_Visualizations", 
        "05_BPMD_Validation"
    ]
    for folder in folders:
        if os.path.exists(folder):
            shutil.rmtree(folder)
            print(f"  Removed: {folder}")
    print("Cleanup complete.")

def main():
    schrodinger_run = os.environ.get('SCHRODINGER', '') + "/run"
    scripts_dir = "pipeline_scripts"
    
    if "--clean" in sys.argv:
        cleanup_results()
        sys.argv.remove("--clean")

    # 1. Detect Docking Files
    input_files = []
    if len(sys.argv) >= 2:
        input_files = [sys.argv[1]]
    else:
        all_mae = glob.glob("DATA/*.mae") + glob.glob("DATA/*.maegz")
        input_files = [f for f in all_mae if "metadynamics" not in f.lower()]
    
    if not input_files:
        print("Warning: No docking files found. Checking for BPMD...")
    
    # --- PHASE 1: COLLATING SIGNIFICANT POSES ---
    all_filtered_mae = "00_Significant_Poses/all_significant_poses.mae"
    os.makedirs("00_Significant_Poses", exist_ok=True)
    
    # We'll use a temporary folder to store intermediate filtered files
    temp_filter_dir = "00_Significant_Poses/temp_parts"
    os.makedirs(temp_filter_dir, exist_ok=True)
    
    filtered_parts = []
    for i, input_file in enumerate(input_files):
        print(f"\n>>> FILTERING SIGNIFICANT POSES FROM: {input_file}")
        part_mae = f"{temp_filter_dir}/part_{i}.mae"
        # Run module0_filter.py but redirect output (requires module0 update or manual move)
        # Actually, I'll just run it and then move the result
        run_cmd(f"{schrodinger_run} {scripts_dir}/module0_filter.py {input_file}")
        if os.path.exists("00_Significant_Poses/filtered_significant_poses.mae"):
            shutil.move("00_Significant_Poses/filtered_significant_poses.mae", part_mae)
            filtered_parts.append(part_mae)

    # Merge all parts if any
    if filtered_parts:
        print("\n>>> MERGING ALL SIGNIFICANT POSES...")
        # Use Schrodinger's structcat or just concatenate if they are MAE
        # For MAE, simple concat works but structcat is safer.
        structcat = os.environ.get('SCHRODINGER', '') + "/utilities/structcat"
        cmd = f"{structcat} -omae {all_filtered_mae} {' '.join(filtered_parts)}"
        run_cmd(cmd)
    
    # --- PHASE 2: INTEGRATED ANALYSIS ---
    if os.path.exists(all_filtered_mae):
        print("\n>>> STARTING INTEGRATED DOCKING ANALYSIS")
        run_cmd(f"{schrodinger_run} {scripts_dir}/analyze_docking.py {all_filtered_mae}")
        run_cmd(f"{schrodinger_run} {scripts_dir}/module1_scores.py {all_filtered_mae}")
        run_cmd(f"{schrodinger_run} {scripts_dir}/module2_interactions.py {all_filtered_mae}")
        run_cmd(f"{schrodinger_run} {scripts_dir}/module2_hotspots.py")
        run_cmd(f"{schrodinger_run} {scripts_dir}/module3_clustering.py {all_filtered_mae}")
        run_cmd(f"{schrodinger_run} {scripts_dir}/module4_visualize.py")
    else:
        print("\n>>> SKIPPING DOCKING ANALYSIS (No significant poses found)")

    # --- PHASE 3: BPMD VALIDATION ---
    bpmd_files = glob.glob("DATA/*metadynamics*.maegz")
    if bpmd_files:
        print("\n>>> STARTING BPMD VALIDATION PHASE")
        run_cmd(f"{schrodinger_run} {scripts_dir}/module5_bpmd.py")
        if os.path.exists("02_Interaction_Analysis"):
            run_cmd(f"{schrodinger_run} {scripts_dir}/module6_consensus.py")
        run_cmd(f"Rscript {scripts_dir}/mdplot.R")
    else:
        print("\n>>> SKIPPING BPMD VALIDATION (No files found)")

    # Cleanup temp parts
    if os.path.exists(temp_filter_dir):
        shutil.rmtree(temp_filter_dir)

    print("\n" + "="*60)
    print("PIPELINE EXECUTION COMPLETE (All data integrated)")
    print("="*60)

if __name__ == "__main__":
    main()
