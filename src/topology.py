"""
Module I: Topology Architect
Responsible for converting chemical information (SMILES) into physical topology (LAMMPS Data File).
"""

import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def generate_topology(smiles: str, output_path: Path, padding: float = 20.0):
    """
    Generates a LAMMPS data file. 
    Embeds each molecule individually and places them randomly.
    """
    print(f"[*] Generating topology for SMILES: {smiles}")
    
    individual_smiles = smiles.split('.')
    all_mols = []
    
    for s in individual_smiles:
        mol = Chem.MolFromSmiles(s)
        if not mol: continue
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=100)
        except:
            pass
        all_mols.append(mol)

    if not all_mols:
        print("Error: No valid molecules generated.")
        sys.exit(1)

    # Element/Bond/Angle Mappings
    element_map = {'C': 1, 'H': 2, 'O': 3, 'N': 4, 'S': 5}
    mass_map = {1: 12.011, 2: 1.008, 3: 15.999, 4: 14.007, 5: 32.06}
    
    def get_bond_type(bond):
        bt = bond.GetBondType()
        if bt == Chem.rdchem.BondType.SINGLE: return 1
        if bt == Chem.rdchem.BondType.DOUBLE: return 2
        if bt == Chem.rdchem.BondType.TRIPLE: return 3
        return 4

    total_atoms = []
    total_bonds = []
    total_angles = []
    
    # Box Calculation
    # Estimate box size based on number of molecules
    box_size = (len(all_mols) ** (1/3)) * 15.0 + padding
    
    atom_offset = 0
    for mol_idx, mol in enumerate(all_mols):
        AllChem.ComputeGasteigerCharges(mol)
        conf = mol.GetConformer()
        
        # Random Translation
        offset = (np.random.rand(3) - 0.5) * box_size * 0.8
        
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            total_atoms.append({
                'id': atom_offset + i + 1,
                'mol': mol_idx + 1,
                'type': element_map.get(atom.GetSymbol(), 1),
                'q': float(atom.GetProp('_GasteigerCharge')) if atom.HasProp('_GasteigerCharge') else 0.0,
                'x': pos.x + offset[0],
                'y': pos.y + offset[1],
                'z': pos.z + offset[2]
            })
            
        for bond in mol.GetBonds():
            total_bonds.append({
                'type': get_bond_type(bond),
                'a1': bond.GetBeginAtomIdx() + atom_offset + 1,
                'a2': bond.GetEndAtomIdx() + atom_offset + 1
            })
            
        # Detect Angles
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
            if len(neighbors) >= 2:
                import itertools
                for a1, a2 in itertools.combinations(neighbors, 2):
                    total_angles.append({
                        'type': 1,
                        'a1': a1 + atom_offset + 1,
                        'a2': idx + atom_offset + 1,
                        'a3': a2 + atom_offset + 1
                    })
        
        atom_offset += mol.GetNumAtoms()

    # Write File
    half_box = box_size / 2.0
    with open(output_path, 'w') as f:
        f.write(f"LAMMPS data file: {smiles}\n\n")
        f.write(f"{len(total_atoms)} atoms\n")
        f.write(f"{len(total_bonds)} bonds\n")
        f.write(f"{len(total_angles)} angles\n\n")
        f.write(f"5 atom types\n4 bond types\n1 angle types\n\n")
        f.write(f"{-half_box:.4f} {half_box:.4f} xlo xhi\n")
        f.write(f"{-half_box:.4f} {half_box:.4f} ylo yhi\n")
        f.write(f"{-half_box:.4f} {half_box:.4f} zlo zhi\n\n")
        
        f.write("Masses\n\n")
        for t, m in sorted(mass_map.items()):
            f.write(f"{t} {m}\n")
        f.write("\n")
        
        f.write("Atoms # full\n\n")
        for a in total_atoms:
            f.write(f"{a['id']} {a['mol']} {a['type']} {a['q']:.5f} {a['x']:.5f} {a['y']:.5f} {a['z']:.5f}\n")
        f.write("\n")
        
        f.write("Bonds\n\n")
        for i, b in enumerate(total_bonds):
            f.write(f"{i+1} {b['type']} {b['a1']} {b['a2']}\n")
        f.write("\n")
        
        f.write("Angles\n\n")
        for i, ang in enumerate(total_angles):
            f.write(f"{i+1} {ang['type']} {ang['a1']} {ang['a2']} {ang['a3']}\n")
            
    print(f"[*] Topology written to {output_path}")
    return element_map, mass_map

            
    print(f"[*] Topology written to {output_path}")
    return element_map, mass_map

            
    print(f"[*] Topology written to {output_path}")
    return element_map, mass_map

if __name__ == "__main__":
    # Test block
    generate_topology("CCO", Path("test.data"))
