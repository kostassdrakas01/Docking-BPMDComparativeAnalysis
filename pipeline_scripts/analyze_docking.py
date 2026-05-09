import os
import sys
import pandas as pd
from schrodinger import structure
from schrodinger.structutils import interactions
from schrodinger.structutils.interactions import hbond

def get_res_info(atom):
    return f"{atom.pdbres.strip()}{atom.resnum}{atom.inscode.strip()}"

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
            # Try to find which ring belongs to the receptor
            # If ps has ring1/ring2 attributes...
            res = None
            for attr in ['ring1', 'ring2', 'ring1_centroid', 'ring2_centroid']:
                if hasattr(ps, attr):
                    obj = getattr(ps, attr)
                    if hasattr(obj, 'atoms'):
                        r = get_res_from_centroid(obj, combined, rec_count)
                        if r: res = r; break
            
            if res:
                itype = "Pi-Stacking"
                if hasattr(ps, 'interaction_type'): itype = f"Pi-Stacking ({ps.interaction_type})"
                res_interactions.append((res, itype))
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
                    elif hasattr(obj, 'index'): # Single atom
                        if obj.index <= rec_count:
                            res = get_res_info(obj); break
            
            if res:
                res_interactions.append((res, "Pi-Cation"))
    except: pass

    return res_interactions

def main(input_file):
    if not os.path.exists(input_file):
        print(f"File not found: {input_file}")
        return

    reader = list(structure.StructureReader(input_file))
    if not reader:
        print("No structures found.")
        return

    receptor = reader[0]
    ligands = reader[1:]
    
    print(f"Receptor: {receptor.title}")
    print(f"Number of ligands/poses: {len(ligands)}")
    
    ligand_groups = {}
    for lig in ligands:
        title = lig.title
        if title not in ligand_groups:
            ligand_groups[title] = []
        ligand_groups[title].append(lig)
    
    results = []
    
    for title, poses in ligand_groups.items():
        print(f"Processing ligand: {title} ({len(poses)} poses)")
        valid_poses = [p for p in poses if 'r_i_glide_gscore' in p.property]
        if not valid_poses: valid_poses = poses
            
        best_pose = min(valid_poses, key=lambda x: x.property.get('r_i_glide_gscore', 0))
        best_score = best_pose.property.get('r_i_glide_gscore', 0)
        
        all_interactions = []
        for pose in poses:
            ints = analyze_interactions(receptor, pose)
            all_interactions.append(set(ints))
            
        freq = {}
        num_poses = len(poses)
        for pose_ints in all_interactions:
            for item in pose_ints:
                freq[item] = freq.get(item, 0) + 1
        
        best_ints = analyze_interactions(receptor, best_pose)
        best_res_types = [f"{res} ({itype})" for res, itype in best_ints]
        
        freq_parts = []
        sorted_freq = sorted(freq.items(), key=lambda x: x[1], reverse=True)
        for (res, itype), count in sorted_freq:
            percent = (count / num_poses) * 100
            freq_parts.append(f"{res} ({itype}): {percent:.1f}%")
            
        results.append({
            'Ligand': title,
            'Best Glide Score': best_score,
            'Num Poses': num_poses,
            'Best Pose Interactions': "; ".join(best_res_types),
            'Interaction Frequency (%)': "; ".join(freq_parts)
        })
        
    df = pd.DataFrame(results)
    # Sort by Ligand to group results
    df = df.sort_values(by='Ligand')
    output_csv = "00_Initial_Docking_Summary/docking_analysis_results.csv"
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df.to_csv(output_csv, index=False)
    
    print("\n" + "="*50)
    print("DOCKING INTERACTION ANALYSIS REPORT")
    print("="*50)
    for res in results:
        print(f"\nLigand: {res['Ligand']}")
        print(f"  Best Glide Score: {res['Best Glide Score']:.3f}")
        print(f"  Best Pose Interactions:")
        if res['Best Pose Interactions']:
            seen = set()
            for i in res['Best Pose Interactions'].split("; "):
                if i not in seen:
                    print(f"    - {i}")
                    seen.add(i)
        else:
            print("    - None detected")
        print(f"  Interaction Frequencies (across {res['Num Poses']} poses):")
        if res['Interaction Frequency (%)']:
            for f in res['Interaction Frequency (%)'].split("; "):
                print(f"    - {f}")
        else:
            print("    - None detected")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: $SCHRODINGER/run analyze_docking.py <input.mae>")
    else:
        main(sys.argv[1])
