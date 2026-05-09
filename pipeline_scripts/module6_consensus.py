import os
import sys
import pandas as pd
from schrodinger import structure
from schrodinger.structutils import interactions
from schrodinger.structutils.interactions import hbond
import glob

# --- ROBUST LIVE INTERACTION ANALYSIS ---
def get_res_info(atom):
    try:
        return f"{atom.pdbres.strip()}{atom.resnum}"
    except:
        return f"RES{atom.resnum}"

def analyze_live_interactions(st):
    res_ints = []
    mols = list(st.molecule)
    if len(mols) < 2: return []
    mol_sizes = [len(m.atom) for m in mols]
    lig_idx = mol_sizes.index(min(mol_sizes))
    lig_atoms = list(mols[lig_idx].getAtomIndices())
    rec_atoms = []
    for i, m in enumerate(mols):
        if i != lig_idx: rec_atoms.extend(list(m.getAtomIndices()))
    rec_set = set(rec_atoms)
    
    # Check H-Bonds and contacts
    try:
        hb_list = hbond.get_hydrogen_bonds(st, atoms1=rec_atoms, st2=st, atoms2=lig_atoms)
        for hb in hb_list:
            at1 = st.atom[hb.atom1]; at2 = st.atom[hb.atom2]
            res = get_res_info(at1) if at1.index in rec_set else get_res_info(at2)
            res_ints.append(f"{res}_H-Bond")
    except: pass
    
    if not res_ints:
        try:
            nearby = set()
            for l in lig_atoms:
                for r in rec_atoms:
                    if st.measure(l, r) < 3.5:
                        nearby.add(get_res_info(st.atom[r]))
            for r in nearby: res_ints.append(f"{r}_Contact")
        except: pass
    return sorted(list(set(res_ints)))

def main():
    bpmd_csv = "05_BPMD_Validation/bpmd_validation_results.csv"
    if not os.path.exists(bpmd_csv): return
        
    df_bpmd = pd.read_csv(bpmd_csv)
    winners = df_bpmd.groupby('Ligand').first().reset_index()
    
    final_results = []
    for _, row in winners.iterrows():
        ligand = row['Ligand']
        source_file = f"DATA/{row['Source_File']}"
        full_title = row['Full_Title']
        
        validated_res = "Fingerprint Missing"
        confidence = "Incomplete Data"
        
        # --- FIX: Mapping BPMD 0-based Pose_Num to Glide 1-based Rank ---
        bpmd_pose = int(row['Pose_Num'])
        glide_rank = bpmd_pose + 1 # Offset correction
        
        fp_file = f"02_Interaction_Analysis/fingerprint_{ligand}.csv"
        if os.path.exists(fp_file):
            df_fp = pd.read_csv(fp_file)
            match = df_fp[df_fp['Emodel_Rank'] == glide_rank]
            if not match.empty:
                cols = [c for c in df_fp.columns if c not in ['Emodel_Rank', 'Ligand']]
                active_ints = [c for c in cols if match.iloc[0][c] == 1]
                if active_ints:
                    validated_res = ", ".join(active_ints)
                    confidence = "Validated"
        
        # Fall-back if still missing
        if validated_res == "Fingerprint Missing" and os.path.exists(source_file):
            try:
                for st in structure.StructureReader(source_file):
                    if st.title == full_title:
                        live = analyze_live_interactions(st)
                        if live:
                            validated_res = ", ".join(live)
                            confidence = "Validated (Live Analysis)"
                        break
            except: pass
            
        final_results.append({
            'Ligand': ligand,
            'BPMD_Winner_Rank': glide_rank,
            'BPMD_PoseScore': row['PoseScore'],
            'BPMD_Persistence': row['Persistence'],
            'Validated_Interactions': validated_res,
            'Validation_Status': confidence
        })
        
    df_final = pd.DataFrame(final_results)
    output_csv = "05_BPMD_Validation/final_validated_interactions.csv"
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df_final.to_csv(output_csv, index=False)
    print("Consensus validation complete with offset correction.")

if __name__ == "__main__":
    main()
