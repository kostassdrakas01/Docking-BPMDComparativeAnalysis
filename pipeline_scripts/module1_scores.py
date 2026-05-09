import os
import sys
import pandas as pd
from schrodinger import structure

def main(input_file):
    output_dir = "01_Score_Analysis"
    os.makedirs(output_dir, exist_ok=True)
    
    reader = list(structure.StructureReader(input_file))
    if not reader: return
    
    receptor = reader[0]
    ligands = reader[1:]
    
    data = []
    for i, lig in enumerate(ligands):
        # Extract properties
        row = {
            'Pose_ID': i + 1,
            'Title': lig.title,
            'GlideScore': lig.property.get('r_i_glide_gscore'),
            'Emodel': lig.property.get('r_i_glide_emodel'),
            'Lipo': lig.property.get('r_i_glide_lipo'),
            'Hbond': lig.property.get('r_i_glide_hbond'),
            'HeavyAtoms': lig.atom_total # Simplification for LE
        }
        
        # Count heavy atoms properly
        heavy_atoms = len([a for a in lig.atom if a.element != 'H'])
        row['HeavyAtoms'] = heavy_atoms
        
        # Calculate LE
        if row['GlideScore'] is not None and heavy_atoms > 0:
            row['LE'] = row['GlideScore'] / heavy_atoms
        else:
            row['LE'] = None
            
        data.append(row)
        
    df = pd.DataFrame(data)
    # Sort by Title to group results
    df = df.sort_values(by='Title')
    df.to_csv(os.path.join(output_dir, "score_distribution.csv"), index=False)
    print(f"Module 1 complete. Results in {output_dir}")

if __name__ == "__main__":
    main(sys.argv[1])
