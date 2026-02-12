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

    # Bond and Angle definitions
    unique_bonds = {} # (type_id) -> (k, r0)
    bond_params_map = {} # (rounded_r0, type) -> type_id
    
    unique_angles = {} # (type_id) -> (k, theta0)
    angle_params_map = {} # (rounded_theta0) -> type_id

    unique_dihedrals = {} # (type_id) -> (k, d, n)
    dihedral_params_map = {} # (k, d, n) -> type_id

    total_atoms = []
    total_bonds = []
    total_angles = []
    total_dihedrals = []
    
    # Box Calculation
    box_size = (len(all_mols) ** (1/3)) * 15.0 + padding
    
    atom_offset = 0
    for mol_idx, mol in enumerate(all_mols):
        AllChem.ComputeGasteigerCharges(mol)
        conf = mol.GetConformer()
        
        # Random Translation
        offset = (np.random.rand(3) - 0.5) * box_size * 0.8
        
        for i, atom in enumerate(mol.GetAtoms()):
            symbol = atom.GetSymbol()
            if symbol not in element_map:
                raise ValueError(f"Error: Element '{symbol}' is not supported by the current force field (CHONS only).")
            
            pos = conf.GetAtomPosition(i)
            total_atoms.append({
                'id': atom_offset + i + 1,
                'mol': mol_idx + 1,
                'type': element_map[symbol],
                'q': float(atom.GetProp('_GasteigerCharge')) if atom.HasProp('_GasteigerCharge') else 0.0,
                'x': pos.x + offset[0],
                'y': pos.y + offset[1],
                'z': pos.z + offset[2]
            })
            
        for bond in mol.GetBonds():
            # Calculate actual length in optimized structure
            a1_idx = bond.GetBeginAtomIdx()
            a2_idx = bond.GetEndAtomIdx()
            p1 = conf.GetAtomPosition(a1_idx)
            p2 = conf.GetAtomPosition(a2_idx)
            r0 = np.sqrt((p1.x-p2.x)**2 + (p1.y-p2.y)**2 + (p1.z-p2.z)**2)
            r0_rounded = round(r0, 2)
            
            bt = bond.GetBondType()
            order = 1
            if bt == Chem.rdchem.BondType.DOUBLE: order = 2
            elif bt == Chem.rdchem.BondType.TRIPLE: order = 3
            
            key = (r0_rounded, order)
            if key not in bond_params_map:
                type_id = len(bond_params_map) + 1
                bond_params_map[key] = type_id
                # Stiffness k based on order
                k = 300.0 * order
                unique_bonds[type_id] = (k, r0_rounded)
            
            total_bonds.append({
                'type': bond_params_map[key],
                'a1': a1_idx + atom_offset + 1,
                'a2': a2_idx + atom_offset + 1
            })

            # Detect Dihedrals around this bond
            neigh1 = [n.GetIdx() for n in mol.GetAtomWithIdx(a1_idx).GetNeighbors() if n.GetIdx() != a2_idx]
            neigh2 = [n.GetIdx() for n in mol.GetAtomWithIdx(a2_idx).GetNeighbors() if n.GetIdx() != a1_idx]
            
            for n1 in neigh1:
                for n2 in neigh2:
                    if n1 == n2: continue
                    
                    # Heuristic parameters for CHONS
                    # Double bond -> planar (n=2, d=-1), Single bond -> staggered (n=3, d=1)
                    if bt == Chem.rdchem.BondType.DOUBLE:
                        k_di = 10.0
                        d_di = -1
                        n_di = 2
                    else:
                        k_di = 1.0
                        d_di = 1
                        n_di = 3
                    
                    di_key = (k_di, d_di, n_di)
                    if di_key not in dihedral_params_map:
                        type_id = len(dihedral_params_map) + 1
                        dihedral_params_map[di_key] = type_id
                        unique_dihedrals[type_id] = di_key
                    
                    total_dihedrals.append({
                        'type': dihedral_params_map[di_key],
                        'a1': n1 + atom_offset + 1,
                        'a2': a1_idx + atom_offset + 1,
                        'a3': a2_idx + atom_offset + 1,
                        'a4': n2 + atom_offset + 1
                    })
            
        # Detect Angles
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
            if len(neighbors) >= 2:
                import itertools
                for a1, a2 in itertools.combinations(neighbors, 2):
                    # Calculate angle
                    p_c = conf.GetAtomPosition(idx)
                    p1 = conf.GetAtomPosition(a1)
                    p2 = conf.GetAtomPosition(a2)
                    
                    v1 = np.array([p1.x-p_c.x, p1.y-p_c.y, p1.z-p_c.z])
                    v2 = np.array([p2.x-p_c.x, p2.y-p_c.y, p2.z-p_c.z])
                    
                    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                    theta0 = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
                    theta0_rounded = round(theta0, 1)
                    
                    if theta0_rounded not in angle_params_map:
                        type_id = len(angle_params_map) + 1
                        angle_params_map[theta0_rounded] = type_id
                        unique_angles[type_id] = (60.0, theta0_rounded)
                        
                    total_angles.append({
                        'type': angle_params_map[theta0_rounded],
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
        f.write(f"{len(total_angles)} angles\n")
        f.write(f"{len(total_dihedrals)} dihedrals\n\n")
        f.write(f"5 atom types\n{len(unique_bonds)} bond types\n{len(unique_angles)} angle types\n{len(unique_dihedrals)} dihedral types\n\n")
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
        f.write("\n")

        if total_dihedrals:
            f.write("Dihedrals\n\n")
            for i, di in enumerate(total_dihedrals):
                f.write(f"{i+1} {di['type']} {di['a1']} {di['a2']} {di['a3']} {di['a4']}\n")
            
    print(f"[*] Topology written to {output_path}")
    return unique_bonds, unique_angles, unique_dihedrals

if __name__ == "__main__":
    # Test block
    generate_topology("CCO", Path("test.data"))
