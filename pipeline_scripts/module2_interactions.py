import os
import sys
import pandas as pd
import numpy as np
from schrodinger import structure
from schrodinger.structutils import interactions
from schrodinger.structutils.interactions import hbond

def get_res_info(atom):
    return f"{atom.pdbres.strip()}{atom.resnum}"

def get_res_from_centroid(centroid, st, rec_count):
    for idx in centroid.atoms:
        if idx <= rec_count:
            return get_res_info(st.atom[idx])
    return None

def analyze_interactions(receptor, ligand):
    res_interactions = []
    rec_count = receptor.atom_total
    combined = receptor.copy()
    combined.extend(ligand)
    rec_atoms = list(range(1, rec_count + 1))
    lig_atoms = list(range(rec_count + 1, combined.atom_total + 1))
    
    # Hydrogen bonds
    try:
        hb_list = hbond.get_hydrogen_bonds(combined, atoms1=rec_atoms, st2=combined, atoms2=lig_atoms)
        for hb in hb_list:
            at1, at2 = (hb.atom1, hb.atom2) if hasattr(hb, 'atom1') else (hb[0], hb[1])
            res = get_res_info(at1) if at1.index <= rec_count else get_res_info(at2)
            res_interactions.append((res, "H-Bond"))
    except: pass

    # Salt bridges
    try:
        sb_list = interactions.get_salt_bridges(combined, group1=rec_atoms, struc2=combined, group2=lig_atoms)
        for sb in sb_list:
            at1, at2 = (sb.atom1, sb.atom2) if hasattr(sb, 'atom1') else (sb[0], sb[1])
            res = get_res_info(at1) if at1.index <= rec_count else get_res_info(at2)
            res_interactions.append((res, "Salt-Bridge"))
    except: pass

    # Pi-stacking
    try:
        ps_list = interactions.find_pi_pi_interactions(combined, struct2=combined, atoms1=rec_atoms, atoms2=lig_atoms)
        for ps in ps_list:
            res = None
            for attr in ['ring1', 'ring2', 'ring1_centroid', 'ring2_centroid']:
                if hasattr(ps, attr):
                    obj = getattr(ps, attr)
                    if hasattr(obj, 'atoms'):
                        r = get_res_from_centroid(obj, combined, rec_count)
                        if r: res = r; break
            if res: res_interactions.append((res, "Pi-Stacking"))
    except: pass

    # Pi-cation
    try:
        pc_list = interactions.find_pi_cation_interactions(combined, struct2=combined, atoms1=rec_atoms, atoms2=lig_atoms)
        for pc in pc_list:
            res = None
            for attr in ['cation_centroid', 'pi_centroid', 'cation', 'ring']:
                if hasattr(pc, attr):
                    obj = getattr(pc, attr)
                    if hasattr(obj, 'atoms'):
                        r = get_res_from_centroid(obj, combined, rec_count)
                        if r: res = r; break
                    elif hasattr(obj, 'index'):
                        if obj.index <= rec_count: res = get_res_info(obj); break
            if res: res_interactions.append((res, "Pi-Cation"))
    except: pass

    return res_interactions

def main(input_file):
    output_dir = "02_Interaction_Analysis"
    os.makedirs(output_dir, exist_ok=True)
    
    reader = list(structure.StructureReader(input_file))
    if not reader: return
    receptor = reader[0]
    ligands = reader[1:]
    
    ligand_groups = {}
    for i, lig in enumerate(ligands):
        if lig.title not in ligand_groups: ligand_groups[lig.title] = []
        ligand_groups[lig.title].append((i, lig))
        
    summary_data = []
    
    for title, poses_info in ligand_groups.items():
        print(f"Analyzing interactions for {title}...")
        num_poses = len(poses_info)
        
        # We'll calculate interaction frequency across ALL poses for the summary
        freq_counts = {}
        
        # For the fingerprint CSV (top 10), we'll store them here
        fp_data = []
        
        for rank, (idx, lig) in enumerate(poses_info):
            ints = analyze_interactions(receptor, lig)
            unique_ints = set(ints) # Set of (res, type)
            
            # Update frequency counts
            for item in unique_ints:
                freq_counts[item] = freq_counts.get(item, 0) + 1
            
            # Store in fingerprint data if it's top 10
            if rank < 10:
                row = {'Emodel_Rank': rank + 1, 'Ligand': title}
                for res, itype in unique_ints:
                    row[f"{res}_{itype}"] = 1
                fp_data.append(row)
        
        # Save fingerprint CSV for top 10
        if fp_data:
            df_fp = pd.DataFrame(fp_data).fillna(0)
            df_fp.to_csv(os.path.join(output_dir, f"fingerprint_{title}.csv"), index=False)
            
        # Add to summary data
        for (res, itype), count in freq_counts.items():
            summary_data.append({
                'Ligand': title,
                'Residue': res,
                'Interaction_Type': itype,
                'Count': count,
                'Total_Poses': num_poses,
                'Frequency (%)': round((count / num_poses) * 100, 1)
            })
            
    if summary_data:
        df_summary = pd.DataFrame(summary_data)
        df_summary = df_summary.sort_values(by=['Ligand', 'Frequency (%)'], ascending=[True, False])
        df_summary.to_csv(os.path.join(output_dir, "interaction_frequency_summary.csv"), index=False)
        
    print(f"Module 2 complete. Results in {output_dir}")

if __name__ == "__main__":
    main(sys.argv[1])
