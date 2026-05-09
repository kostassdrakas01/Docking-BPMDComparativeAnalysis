import os
import sys
from schrodinger import structure

def main(input_file):
    output_dir = "00_Significant_Poses"
    os.makedirs(output_dir, exist_ok=True)
    output_mae = os.path.join(output_dir, "filtered_significant_poses.mae")
    os.makedirs(os.path.dirname(output_mae), exist_ok=True)
    
    reader = list(structure.StructureReader(input_file))
    if not reader:
        print("No structures found.")
        return
        
    receptor = reader[0]
    ligands = reader[1:]
    
    # Group by ligand
    ligand_groups = {}
    for lig in ligands:
        if lig.title not in ligand_groups:
            ligand_groups[lig.title] = []
        ligand_groups[lig.title].append(lig)
        
    significant_poses = [receptor]
    
    print("\n" + "="*50)
    print("POSE FILTERING (Theory: Emodel-based significance)")
    print("="*50)
    print("According to Glide theory, Emodel is the most reliable metric")
    print("for selecting the best binding pose within a set of results.")
    print("-" * 50)
    
    for title, poses in ligand_groups.items():
        # Sort by Emodel (lower is better)
        # Note: r_i_glide_emodel is usually negative
        poses_sorted = sorted(poses, key=lambda x: x.property.get('r_i_glide_emodel', 0))
        
        # Keep top 10 significant poses per ligand
        top_n = 10
        filtered = poses_sorted[:top_n]
        
        print(f"Ligand {title}: Kept {len(filtered)} significant poses (out of {len(poses)})")
        significant_poses.extend(filtered)
        
    # Write the filtered poses
    writer = structure.StructureWriter(output_mae)
    for st in significant_poses:
        writer.append(st)
    writer.close()
    
    print(f"\nFiltered poses saved to: {output_mae}")
    return output_mae

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: $SCHRODINGER/run module0_filter.py <input.mae>")
    else:
        main(sys.argv[1])
