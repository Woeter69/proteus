import pandas as pd
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import Descriptors

def generate_dummy_data(n_samples=500, output_path="data/training_data.csv"):
    """
    Generates synthetic polymer data for testing the ML pipeline.
    Hypothesis: Rg ~ sqrt(Chain Length) * Monomer Size + Noise
    """
    print(f"[*] Generating {n_samples} synthetic samples...")
    
    # Common Monomers
    monomers = [
        "CC", "C=C", "C=CC1=CC=CC=C1", "C(=O)O", "C(C)O", "COC", "C1CCCCC1", 
        "C(=O)NC", "C(F)(F)C(F)(F)", "C=CC#N"
    ]
    
    data = []
    
    for _ in range(n_samples):
        smiles = np.random.choice(monomers)
        chain_length = np.random.randint(10, 200)
        
        # Calculate simplistic 'Ground Truth' based on physics
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.MolWt(mol)
        
        # Flory theory: Rg ~ N^0.6 (approx for good solvent)
        # Adding some randomness
        rg = (chain_length ** 0.6) * (mw / 50.0) + np.random.normal(0, 0.5)
        
        data.append({
            "smiles": smiles,
            "chain_length": chain_length,
            "rg": abs(rg)
        })
        
    df = pd.DataFrame(data)
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"[*] Saved to {output_path}")
    print(df.head())

if __name__ == "__main__":
    generate_dummy_data()
