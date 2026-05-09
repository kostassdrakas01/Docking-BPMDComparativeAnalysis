import os
import sys
import pandas as pd
from schrodinger import structure
import glob

def main():
    output_dir = "05_BPMD_Validation"
    os.makedirs(output_dir, exist_ok=True)
    
    bpmd_files = glob.glob("DATA/*metadynamics*out.maegz")
    if not bpmd_files:
        print("No BPMD files found.")
        return
        
    bpmd_data = []
    for f in bpmd_files:
        print(f"Extracting BPMD results from {f}...")
        try:
            for st in structure.StructureReader(f):
                full_title = st.title
                
                # Standardized Ligand Naming Logic
                ligand = full_title
                if "analog1" in full_title.lower(): ligand = "analog1"
                elif "5oh" in full_title.lower() or "hydroxystaurosporine" in full_title.lower(): 
                    ligand = "5’-hydroxystaurosporine"
                elif "staurosporine" in full_title.lower(): ligand = "staurosporine"
                elif "midostaurin" in full_title.lower(): ligand = "midostaurin"
                    
                row = {
                    'Ligand': ligand,
                    'Full_Title': full_title,
                    'Pose_Num': st.property.get('i_psp_BindingPoseMetadynamics_Pose_Number'),
                    'PoseScore': st.property.get('r_psp_MetadynamicsBinding_PoseScore'),
                    'Persistence': st.property.get('r_psp_MetadynamicsBinding_Persistence'),
                    'CompScore': st.property.get('r_psp_MetadynamicsBinding_CompScore'),
                    'GlideScore': st.property.get('r_i_glide_gscore'),
                    'Source_File': os.path.basename(f)
                }
                bpmd_data.append(row)
        except Exception as e:
            print(f"Error reading {f}: {e}")
            
    if not bpmd_data: return
    df = pd.DataFrame(bpmd_data)
    df = df.sort_values(by=['Ligand', 'PoseScore'], ascending=[True, False])
    output_csv = os.path.join(output_dir, "bpmd_validation_results.csv")
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df.to_csv(output_csv, index=False)
    print(f"BPMD Summary generated in {output_dir}")

if __name__ == "__main__":
    main()
