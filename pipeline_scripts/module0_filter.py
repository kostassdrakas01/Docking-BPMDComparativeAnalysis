import os
import sys
from schrodinger import structure

def main(input_file):
    output_dir = "00_Significant_Poses"
    output_mae = os.path.join(output_dir, "filtered_significant_poses.mae")
    os.makedirs(output_dir, exist_ok=True)
    
    reader = list(structure.StructureReader(input_file))
    if not reader: return
        
    receptor = reader[0]
    raw_ligands = reader[1:]
    ligands = []
    
    print("\n>>> Filtering out potential receptor structures from ligand list...")
    for lig in raw_ligands:
        # ATOM-BASED EXCLUSION: If it has > 500 atoms, it's a protein/receptor, not a ligand
        if lig.atom_total > 500:
            print(f"Skipping receptor-like structure: {lig.title} ({lig.atom_total} atoms)")
            continue
            
        t = lig.title.lower()
        if "4r67" in t or "minimized" in t or "receptor" in t:
            print(f"Skipping title-matched receptor: {lig.title}")
            continue
            
        ligands.append(lig)
        
    ligand_groups = {}
    for lig in ligands:
        if lig.title not in ligand_groups: ligand_groups[lig.title] = []
        ligand_groups[lig.title].append(lig)
        
    significant_poses = [receptor]
    
    print("\n" + "="*50)
    print("POSE FILTERING (Theory: Emodel-based significance)")
    print("="*50)
    
    for title, poses in ligand_groups.items():
        poses_sorted = sorted(poses, key=lambda x: x.property.get('r_i_glide_emodel', 0))
        top_n = 10
        filtered = poses_sorted[:top_n]
        print(f"Ligand {title}: Kept {len(filtered)} significant poses (out of {len(poses)})")
        significant_poses.extend(filtered)
        
    writer = structure.StructureWriter(output_mae)
    for st in significant_poses:
        writer.append(st)
    writer.close()
    print(f"\nFiltered poses saved to: {output_mae}")

if __name__ == "__main__":
    if len(sys.argv) >= 2:
        main(sys.argv[1])
