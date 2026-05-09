import os
import sys
import pandas as pd
import numpy as np
from schrodinger import structure
from schrodinger.structutils import rmsd
from rdkit.ML.Cluster import Butina

def cluster_poses(poses, threshold=2.0):
    n = len(poses)
    if n <= 1: return [0] * n
    
    # Calculate distance matrix (RMSD)
    dists = []
    for i in range(n):
        for j in range(i):
            # Use Schrödinger's RMSD calculation
            # This handles atom mapping if needed, but for poses of same ligand, 
            # we usually assume same atom order.
            try:
                val = rmsd.calculate_rmsd(poses[i], poses[j], use_symmetry=True)
            except:
                val = 100.0 # High value if fails
            dists.append(val)
            
    clusters = Butina.ClusterData(dists, n, threshold, isDistData=True)
    
    cluster_ids = [0] * n
    for cluster_id, cluster in enumerate(clusters):
        for pose_idx in cluster:
            cluster_ids[pose_idx] = cluster_id
    return cluster_ids

def main(input_file):
    output_dir = "03_Clustering_Analysis"
    os.makedirs(output_dir, exist_ok=True)
    
    reader = list(structure.StructureReader(input_file))
    if not reader: return
    ligands = reader[1:]
    
    ligand_groups = {}
    for i, lig in enumerate(ligands):
        if lig.title not in ligand_groups: ligand_groups[lig.title] = []
        ligand_groups[lig.title].append(lig)
        
    all_results = []
    
    for title, poses in ligand_groups.items():
        print(f"Clustering {title}...")
        if not poses: continue
        
        cluster_ids = cluster_poses(poses, threshold=2.0)
        
        for i, (pose, cid) in enumerate(zip(poses, cluster_ids)):
            all_results.append({
                'Ligand': title,
                'Pose_ID': i + 1,
                'Cluster_ID': cid,
                'GlideScore': pose.property.get('r_i_glide_gscore')
            })
            
    df = pd.DataFrame(all_results)
    # Sort by Ligand to group results
    df = df.sort_values(by=['Ligand', 'Pose_ID'])
    output_csv = os.path.join(output_dir, "clusters.csv")
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df.to_csv(output_csv, index=False)
    print(f"Module 3 complete. Results in {output_dir}")

if __name__ == "__main__":
    main(sys.argv[1])
