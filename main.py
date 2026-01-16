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

def main():
    parser = argparse.ArgumentParser(description="Proteus: Polymer Nanoprecipitation Simulator")
    parser.add_argument("--smiles", type=str, required=True, help="SMILES string of the polymer")
    parser.add_argument("--name", type=str, default="simulation", help="Name of the simulation run")
    parser.add_argument("--steps", type=int, default=10000, help="Number of simulation steps")
    parser.add_argument("--count", type=int, default=1, help="Number of copies of the molecule to simulate")
    parser.add_argument("--render", action="store_true", help="Render a GIF of the simulation (requires Ovito)")
    parser.add_argument("--version", action="version", version="%(prog)s v1.0.0")
    
    args = parser.parse_args()
    
    # Handle Molecule Count
    if args.count > 1:
        print(f"[*] Replicating molecule {args.count} times...")
        # "CCO" * 3 -> "CCO.CCO.CCO"
        args.smiles = ".".join([args.smiles] * args.count)
    
    # Setup Paths
    base_dir = Path(__file__).parent
    output_dir = base_dir / "output" / args.name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    data_file = output_dir / "polymer.data"
    input_file = output_dir / "simulation.in"
    log_file = output_dir / "simulation.log"
    dump_file = output_dir / "trajectory.dump"
    gif_file = output_dir / "animation.gif"
    
    print("=========================================")
    print(f"Proteus Pipeline: {args.name}")
    print("=========================================")
    
    # 1. Topology
    try:
        topology.generate_topology(args.smiles, data_file)
    except Exception as e:
        print(f"Topology Generation Failed: {e}")
        sys.exit(1)
        
    # 2. Simulation Setup & Run
    try:
        simulation.generate_input_file(data_file, input_file, dump_file, steps=args.steps)
        simulation.run_simulation(input_file, log_file)
    except Exception as e:
        print(f"Simulation Failed: {e}")
        sys.exit(1)
        
    # 3. Analysis
    try:
        analysis.analyze_results(log_file)
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
