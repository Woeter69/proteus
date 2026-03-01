"""
Proteus: Automated Polymer Nanoprecipitation Pipeline
Main Controller
"""

import argparse
import sys
from pathlib import Path
from rdkit import Chem

# Add src to python path if needed
sys.path.append(str(Path(__file__).parent / "src"))

from src import topology, simulation, analysis, visualization

def validate_smiles(smiles: str, label: str = "Molecule"):
    """Validates a SMILES string using RDKit."""
    if not smiles: return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Error: Invalid {label} SMILES: {smiles}")
        sys.exit(1)
    return smiles

def main():
    parser = argparse.ArgumentParser(description="Proteus: Polymer Nanoprecipitation Simulator")
    parser.add_argument("--smiles", type=str, required=False, help="SMILES string of the polymer")
    parser.add_argument("--name", type=str, default="simulation", help="Name of the simulation run")
    parser.add_argument("--steps", type=int, default=10000, help="Number of simulation steps")
    parser.add_argument("--count", type=int, default=1, help="Number of copies of the molecule to simulate")
    parser.add_argument("--payload", type=str, default=None, help="SMILES of the drug/payload molecule (Optional)")
    parser.add_argument("--payload_count", type=int, default=1, help="Number of payload molecules")
    parser.add_argument("--render", action="store_true", help="Render a GIF of the simulation (requires Ovito)")
    parser.add_argument("--plot", action="store_true", default=True, help="Generate a stability plot (default: True)")
    parser.add_argument("--no-plot", action="store_false", dest="plot", help="Disable stability plot generation")
    
    # Advanced Physics Parameters
    parser.add_argument("--temp", type=float, default=300.0, help="Temperature (K)")
    parser.add_argument("--damp", type=float, default=20.0, help="Langevin damping parameter (fs). Lower = higher viscosity.")
    parser.add_argument("--epsilon", type=float, default=None, help="LJ Epsilon (interaction strength). Defaults to OPLS-AA.")
    parser.add_argument("--sigma", type=float, default=None, help="LJ Sigma (particle size). Defaults to OPLS-AA.")
    parser.add_argument("--timestep", type=float, default=1.0, help="Simulation timestep (fs)")
    parser.add_argument("--padding", type=float, default=20.0, help="Simulation box padding (Angstroms)")
    parser.add_argument("--gpus", type=int, default=1, help="Number of GPUs to use for simulation")

    parser.add_argument("--version", action="version", version="%(prog)s v1.0.1")
    
    args = parser.parse_args()

    if not args.smiles:
        parser.print_help()
        sys.exit(1)
    
    # 0. Validation & Input Prep
    validate_smiles(args.smiles, "Polymer")
    validate_smiles(args.payload, "Payload")
    
    # Construct final system SMILES
    polymer_list = [args.smiles] * args.count
    payload_list = [args.payload] * args.payload_count if args.payload else []
    system_smiles = ".".join(polymer_list + payload_list)
    
    # Setup Paths
    output_dir = Path(__file__).parent / "output" / args.name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    paths = {
        "data": output_dir / "polymer.data",
        "input": output_dir / "simulation.in",
        "log": output_dir / "simulation.log",
        "dump": output_dir / "trajectory.dump",
        "gif": output_dir / "animation.gif",
        "plot": output_dir / "stability.png" if args.plot else None
    }
    
    print("=" * 40)
    print(f"Proteus Pipeline: {args.name}")
    print(f"System: {args.count}x Polymer, {args.payload_count if args.payload else 0}x Payload")
    print("=" * 40)
    
    try:
        # 1. Topology
        print("[*] Phase 1: Topology Generation")
        bond_p, angle_p, dihedral_p = topology.generate_topology(system_smiles, paths["data"], padding=args.padding)
        
        # 2. Simulation
        print("[*] Phase 2: LAMMPS Simulation")
        simulation.generate_input_file(
            paths["data"], paths["input"], paths["dump"], 
            steps=args.steps, temp=args.temp, damp=args.damp,
            epsilon=args.epsilon, sigma=args.sigma, timestep=args.timestep,
            bond_params=bond_p, angle_params=angle_p, dihedral_params=dihedral_p
        )
        simulation.run_simulation(paths["input"], paths["log"], gpus=args.gpus)
        
        # 3. Analysis
        print("[*] Phase 3: Analytics")
        analysis.analyze_results(
            paths["log"], 
            output_plot=paths["plot"],
            polymer_count=args.count,
            payload_count=args.payload_count if args.payload else 0,
            dump_path=paths["dump"]
        )

        # 4. Visualization (Optional)
        if args.render:
            print("[*] Phase 4: Visualization")
            visualization.render_trajectory(dump_path=paths["dump"], output_gif=paths["gif"])

        print("=" * 40)
        print("Pipeline Finished Successfully.")
        print(f"Results located in: {output_dir}")

    except Exception as e:
        print(f"\n[!] Pipeline Failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
