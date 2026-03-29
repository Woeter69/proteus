"""
Proteus: Automated Polymer Nanoprecipitation Pipeline
Main Controller
"""

import argparse
import sys
import os
from pathlib import Path
from rdkit import Chem

# Add src to python path if needed
sys.path.append(str(Path(__file__).parent / "src"))

from src import topology, simulation, analysis, visualization, report

def validate_smiles(smiles: str, label: str = "Molecule"):
    """Validates a SMILES string using RDKit."""
    if not smiles: return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Error: Invalid {label} SMILES: {smiles}")
        return None
    return smiles

def run_pipeline(args):
    """
    Orchestrates the Proteus pipeline for a single simulation run.
    Returns a dictionary of results.
    """
    # 0. Validation & Input Prep
    if not validate_smiles(args.smiles, "Polymer"):
        raise ValueError(f"Invalid Polymer SMILES: {args.smiles}")
    if args.payload and not validate_smiles(args.payload, "Payload"):
        raise ValueError(f"Invalid Payload SMILES: {args.payload}")
    
    # Construct final system SMILES
    polymer_list = [args.smiles] * args.count
    payload_list = [args.payload] * args.payload_count if args.payload else []
    system_smiles = ".".join(polymer_list + payload_list)
    
    # Setup Paths
    output_dir = Path(os.getcwd()) / "output" / args.name
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
    results = analysis.analyze_results(
        paths["log"], 
        output_plot=paths["plot"],
        polymer_count=args.count,
        payload_count=args.payload_count if args.payload else 0,
        dump_path=paths["dump"]
    )

    # 4. Visualization (Optional)
    if args.render:
        print("[*] Phase 4: Visualization")
        visualization.render_trajectory(dump_path=paths["dump"], output_gif=paths["gif"], cut=args.cut)

    # 5. Automated Lab Notebook (Optional)
    if args.report:
        print("[*] Phase 5: Automated Lab Notebook")
        report.generate_report(
            output_dir,
            args.name,
            args.smiles,
            args.steps,
            args.temp,
            results["rg"],
            efficiency=results.get("efficiency"),
            plot_path=paths["plot"]
        )
    
    if args.json:
        print("[*] Generating JSON Report")
        import json
        json_path = output_dir / "report.json"
        report_data = {
            "name": args.name,
            "smiles": args.smiles,
            "steps": args.steps,
            "temp": args.temp,
            "rg": results["rg"],
            "efficiency": results.get("efficiency"),
            "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        with open(json_path, "w") as f:
            json.dump(report_data, f, indent=4)
        print(f"[*] JSON report saved to: {json_path}")

    print("=" * 40)
    print(f"Pipeline Finished Successfully for: {args.name}")
    return {
        "name": args.name,
        "smiles": args.smiles,
        "rg": results["rg"],
        "efficiency": results.get("efficiency"),
        "output_dir": str(output_dir)
    }

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
    parser.add_argument("--report", action="store_true", help="Generate a professional PDF lab report")
    parser.add_argument("--json", action="store_true", help="Generate a machine-readable JSON report")
    parser.add_argument("--hts", type=str, default=None, help="Path to a CSV file for High-Throughput Screening")
    parser.add_argument("--rank-by", type=str, choices=["rg", "efficiency"], default="efficiency", help="HTS ranking priority (default: efficiency)")
    parser.add_argument("--solvent", type=str, choices=["water", "ethanol", "dmso", "custom"], default="custom", help="Predefined solvent profiles (sets temp/damp)")
    
    # Advanced Physics Parameters
    parser.add_argument("--temp", type=float, default=None, help="Temperature (K)")
    parser.add_argument("--damp", type=float, default=None, help="Langevin damping parameter (fs). Lower = higher viscosity.")
    parser.add_argument("--epsilon", type=float, default=None, help="LJ Epsilon (interaction strength). Defaults to OPLS-AA.")
    parser.add_argument("--sigma", type=float, default=None, help="LJ Sigma (particle size). Defaults to OPLS-AA.")
    parser.add_argument("--timestep", type=float, default=1.0, help="Simulation timestep (fs)")
    parser.add_argument("--padding", type=float, default=20.0, help="Simulation box padding (Angstroms)")
    parser.add_argument("--gpus", type=int, default=1, help="Number of GPUs to use for simulation")
    parser.add_argument("--cut", action="store_true", help="Render a cross-section (clipping plane) of the nanoparticle")

    parser.add_argument("--version", action="version", version="%(prog)s v1.2.1")
    
    args = parser.parse_args()

    # Predefined Solvent Profiles
    solvent_profiles = {
        "water": {"temp": 298.15, "damp": 50.0},
        "ethanol": {"temp": 298.15, "damp": 20.0},
        "dmso": {"temp": 310.15, "damp": 100.0}
    }

    if args.solvent != "custom":
        profile = solvent_profiles[args.solvent]
        if args.temp is None: args.temp = profile["temp"]
        if args.damp is None: args.damp = profile["damp"]
    else:
        # Default custom values if not provided
        if args.temp is None: args.temp = 300.0
        if args.damp is None: args.damp = 20.0

    if args.hts:
        from src import hts
        hts.run_hts(args.hts, args)
        return

    if not args.smiles:
        parser.print_help()
        sys.exit(1)
    
    try:
        run_pipeline(args)
    except Exception as e:
        print(f"\n[!] Pipeline Failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
