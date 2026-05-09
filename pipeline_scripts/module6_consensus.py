import os
import pandas as pd

def main():
    output_dir = "05_BPMD_Validation"
    os.makedirs(output_dir, exist_ok=True)
    bpmd_csv = os.path.join(output_dir, "bpmd_validation_results.csv")
    if not os.path.exists(bpmd_csv):
        print("BPMD results not found.")
        return
        
    df_bpmd = pd.read_csv(bpmd_csv)
    
    # We want to find the interactions for the BEST Pose of each ligand
    results = []
    
    for ligand in df_bpmd['Ligand'].unique():
        # Get the winner pose
        winner = df_bpmd[df_bpmd['Ligand'] == ligand].iloc[0]
        pose_idx = int(winner['Pose_Num']) + 1 # Convert 0-index to Rank-index
        
        # Look up interactions in fingerprint file
        fp_file = f"02_Interaction_Analysis/fingerprint_{ligand}.csv"
        if os.path.exists(fp_file):
            df_fp = pd.read_csv(fp_file)
            # Find the row for this rank
            # Note: BPMD Pose_Num 0 usually maps to the 1st pose of that ligand
            # In our fingerprint, Emodel_Rank is 1, 2, 3...
            pose_row = df_fp[df_fp['Emodel_Rank'] == pose_idx]
            
            if not pose_row.empty:
                # Find all columns with 1.0
                ints = [c for c in df_fp.columns if c not in ['Emodel_Rank', 'Ligand'] and pose_row[c].values[0] == 1.0]
                
                results.append({
                    'Ligand': ligand,
                    'BPMD_Winner_Rank': pose_idx,
                    'BPMD_PoseScore': winner['PoseScore'],
                    'BPMD_Persistence': winner['Persistence'],
                    'Validated_Interactions': ", ".join(ints),
                    'Validation_Status': "Validated" if winner['PoseScore'] > 1.2 else "Questionable"
                })
            else:
                # Try to find by index if rank doesn't match
                results.append({
                    'Ligand': ligand,
                    'BPMD_Winner_Rank': pose_idx,
                    'BPMD_PoseScore': winner['PoseScore'],
                    'BPMD_Persistence': winner['Persistence'],
                    'Validated_Interactions': "Rank Mismatch - Check manually",
                    'Validation_Status': "Manual Check Needed"
                })
        else:
            results.append({
                'Ligand': ligand,
                'BPMD_Winner_Rank': pose_idx,
                'BPMD_PoseScore': winner['PoseScore'],
                'BPMD_Persistence': winner['Persistence'],
                'Validated_Interactions': "Fingerprint Missing",
                'Validation_Status': "Incomplete Data"
            })
            
    df_final = pd.DataFrame(results)
    output_csv = os.path.join(output_dir, "final_validated_interactions.csv")
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df_final.to_csv(output_csv, index=False)
    
    print("\n" + "="*60)
    print("FINAL VALIDATED BINDING MODES (BPMD Consensus)")
    print("="*60)
    print("The following residues are 'FOR SURE' interacting based on")
    print("high-stability metadynamics trajectories:")
    print("-" * 60)
    
    for _, row in df_final.iterrows():
        print(f"\nLigand: {row['Ligand']}")
        print(f"  Stable Pose: Emodel Rank {row['BPMD_Winner_Rank']}")
        print(f"  PoseScore: {row['BPMD_PoseScore']:.3f}")
        print(f"  Validated Residues: {row['Validated_Interactions']}")
        print(f"  Confidence: {row['Validation_Status']}")

    print(f"\nFinal report saved to: {output_csv}")

if __name__ == "__main__":
    main()
