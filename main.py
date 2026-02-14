"""
Proteus: Automated Polymer Nanoprecipitation Pipeline
Main Controller
"""

import argparse
import sys
from pathlib import Path

# Add src to python path if needed, though local import usually works if running from root
sys.path.append(str(Path(__file__).parent / "src"))

from src import topology, simulation, analysis, visualization
from src.ml import predict

def main():
    parser = argparse.ArgumentParser(description="Proteus: Polymer Nanoprecipitation Simulator")
    parser.add_argument("--smiles", type=str, required=False, help="SMILES string of the polymer")
    parser.add_argument("--predict", action="store_true", help="Use AI to predict Rg instead of simulating")
    parser.add_argument("--model", type=str, default="default", help="Name of the AI model to use (e.g. 'v3')")
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

    parser.add_argument("--version", action="version", version="%(prog)s v1.0.0")
    
    args = parser.parse_args()

    # 0. AI Prediction Mode
    if args.predict:
        if not args.smiles:
            print("Error: --smiles is required for prediction.")
            sys.exit(1)
        print("=========================================")
        print(f"Proteus AI Inference")
        print(f"Model: {args.model}")
        print("=========================================")
        rg = predict.predict_rg(args.smiles, args.model)
        if rg:
            print(f"[*] Molecule: {args.smiles}")
            print(f"[*] Predicted Radius of Gyration (Rg): {rg:.4f} Å")
        sys.exit(0)
    
    if not args.smiles:
        parser.print_help()
        sys.exit(1)
    
    # Handle Molecule Count
    monomer_smiles = ".".join([args.smiles] * args.count)
    
    # Handle Payload
    if args.payload:
        print(f"[*] Adding Payload: {args.payload} (x{args.payload_count})")
        payload_smiles = ".".join([args.payload] * args.payload_count)
        args.smiles = f"{monomer_smiles}.{payload_smiles}"
    else:
        args.smiles = monomer_smiles
    
    # Setup Paths
    base_dir = Path(__file__).parent
    output_dir = base_dir / "output" / args.name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    data_file = output_dir / "polymer.data"
    input_file = output_dir / "simulation.in"
    log_file = output_dir / "simulation.log"
    dump_file = output_dir / "trajectory.dump"
    gif_file = output_dir / "animation.gif"
    plot_file = output_dir / "stability.png"
    
    print("=========================================")
    print(f"Proteus Pipeline: {args.name}")
    print("=========================================")
    
    # 1. Topology
    try:
        bond_params, angle_params, dihedral_params = topology.generate_topology(args.smiles, data_file, padding=args.padding)
    except Exception as e:
        print(f"Topology Generation Failed: {e}")
        sys.exit(1)
        
    # 2. Simulation Setup & Run
    try:
        simulation.generate_input_file(
            data_file, 
            input_file, 
            dump_file, 
            steps=args.steps,
            temp=args.temp,
            damp=args.damp,
            epsilon=args.epsilon,
            sigma=args.sigma,
            timestep=args.timestep,
            bond_params=bond_params,
            angle_params=angle_params,
            dihedral_params=dihedral_params
        )
        simulation.run_simulation(input_file, log_file)
    except Exception as e:
        print(f"Simulation Failed: {e}")
        sys.exit(1)
        
    # 3. Analysis
    try:
        plot_path = plot_file if args.plot else None
        analysis.analyze_results(log_file, output_plot=plot_path)
    except Exception as e:
        print(f"Analysis Failed: {e}")
        sys.exit(1)

    # 4. Visualization (Optional)
    if args.render:
        try:
            visualization.render_trajectory(dump_path=dump_file, output_gif=gif_file)
        except Exception as e:
            print(f"Visualization Failed: {e}")

    print("=========================================")
    print("Pipeline Finished Successfully.")
    print(f"Results located in: {output_dir}")

if __name__ == "__main__":
    main()
