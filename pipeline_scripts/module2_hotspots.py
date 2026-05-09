import os
import sys
import pandas as pd

def main():
    input_csv = "02_Interaction_Analysis/interaction_frequency_summary.csv"
    if not os.path.exists(input_csv):
        print(f"File not found: {input_csv}")
        return

    df = pd.read_csv(input_csv)
    
    # Identify Hotspots (Highest frequency interactions per ligand)
    # We'll label them based on their frequency
    def get_confidence(freq):
        if freq == 100: return "Absolute (100%)"
        if freq >= 75: return "Very High (>75%)"
        if freq >= 50: return "High (>50%)"
        if freq >= 25: return "Medium (>25%)"
        return "Low (<25%)"

    df['Confidence'] = df['Frequency (%)'].apply(get_confidence)
    
    # Create a simplified report of "Key Conserved Residues"
    # We want to show which residues are the 'anchor' points
    hotspots = df.sort_values(by=['Ligand', 'Frequency (%)'], ascending=[True, False])
    
    output_csv = "02_Interaction_Analysis/key_hotspot_residues.csv"
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    hotspots.to_csv(output_csv, index=False)
    
    print("\n" + "="*60)
    print("SCIENTIFIC SUMMARY: KEY BINDING HOTSPOTS (Structural Anchors)")
    print("="*60)
    print("In computational docking theory, residues that interact with a ligand")
    print("across multiple poses (Consensus Interactions) are identified as ")
    print("'Hotspots' critical for binding stability. These residues act as")
    print("structural anchors that maintain the ligand's orientation.")
    print("-" * 60)

    for ligand in hotspots['Ligand'].unique():
        lig_df = hotspots[hotspots['Ligand'] == ligand].head(3) # Show top 3
        print(f"\nLigand: {ligand}")
        for _, row in lig_df.iterrows():
            print(f"  - {row['Residue']} ({row['Interaction_Type']}): {row['Frequency (%)']}% frequency -> [{row['Confidence']}]")
            
    print("\nFull report saved to: 02_Interaction_Analysis/key_hotspot_residues.csv")

if __name__ == "__main__":
    main()
