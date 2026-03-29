"""
Module VII: High-Throughput Screening (HTS)
Batch processes multiple simulations from a CSV file and ranks them.
"""

import csv
import sys
from pathlib import Path
from copy import deepcopy

def run_hts(csv_path: str, global_args):
    """
    Orchestrates a High-Throughput Screening (HTS) campaign from a CSV dataset.

    This module:
    1. Reads a CSV file where each row defines a simulation (name, smiles, count, etc.).
    2. Sequentially runs the full Proteus pipeline for each row.
    3. Aggregates results into a list of physical metrics.
    4. Automatically generates an 'hts_errors.log' if any simulations fail.
    5. Passes results to 'rank_and_summarize' for final leaderboard generation.

    Args:
        csv_path (str): Path to the input CSV file.
        global_args (argparse.Namespace): Global configuration from main.py.
    """
    from main import run_pipeline
    
    csv_file = Path(csv_path)
    if not csv_file.exists():
        print(f"Error: Batch CSV file not found: {csv_path}")
        sys.exit(1)
        
    print(f"[*] Starting High-Throughput Screening (HTS): {csv_path}")
    
    results = []
    errors = []
    
    try:
        with open(csv_file, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Create a fresh args object for this run
                run_args = deepcopy(global_args)
                
                # Update with CSV values (handle defaults/types)
                run_args.name = row.get('name', f"HTS_{len(results)+1}")
                run_args.smiles = row.get('smiles')
                run_args.count = int(row.get('count', run_args.count))
                run_args.payload = row.get('payload')
                if run_args.payload == "": run_args.payload = None
                run_args.payload_count = int(row.get('payload_count', run_args.payload_count))
                run_args.steps = int(row.get('steps', run_args.steps))
                
                if not run_args.smiles:
                    print(f"Skipping row: missing SMILES in {row}")
                    errors.append(f"{run_args.name}: Missing SMILES")
                    continue
                
                try:
                    result = run_pipeline(run_args)
                    results.append(result)
                except Exception as e:
                    print(f"[!] Run failed for {run_args.name}: {e}")
                    errors.append(f"{run_args.name}: {e}")
                    continue
                    
        # Log errors to file
        if errors:
            error_log = csv_file.parent / "hts_errors.log"
            with open(error_log, "w") as f:
                f.write("\n".join(errors))
            print(f"[!] HTS Errors logged to: {error_log.name}")

        if not results:
            print("No simulations completed successfully in HTS run.")
            return

        # Rank Results
        rank_and_summarize(results, csv_file.parent, global_args)
        
    except Exception as e:
        print(f"Error during HTS screening: {e}")

def rank_and_summarize(results, output_dir: Path, global_args):
    """
    Sorts and ranks simulation results based on user-defined priority.

    The ranking uses a two-tier sort:
    - If priority is 'efficiency': Primary = Encapsulation % (Desc), Secondary = Rg (Asc).
    - If priority is 'rg': Primary = Rg (Asc), Secondary = Efficiency % (Desc).

    Args:
        results (list): List of dictionaries containing individual simulation results.
        output_dir (Path): Directory to save the summary CSV.
        global_args (argparse.Namespace): Contains 'rank_by' choice and other metadata.
    """
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_path = output_dir / f"hts_summary_{timestamp}.csv"
    
    print(f"[*] Ranking {len(results)} results...")
    
    # Ranking logic
    if global_args.rank_by == "efficiency":
        # Primary: Efficiency (Desc), Secondary: Rg (Asc)
        ranked = sorted(
            results, 
            key=lambda x: (-(x.get('efficiency') or 0.0), x['rg'])
        )
    else:
        # Primary: Rg (Asc), Secondary: Efficiency (Desc)
        ranked = sorted(
            results, 
            key=lambda x: (x['rg'], -(x.get('efficiency') or 0.0))
        )
    
    # Save to CSV
    keys = ["name", "rg", "efficiency", "smiles", "output_dir"]
    with open(summary_path, mode='w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for r in ranked:
            writer.writerow({k: r.get(k) for k in keys})
            
    print(f"[*] HTS Ranking Summary saved to: {summary_path}")
    
    print("\n" + "="*40)
    print("HTS LEADERBOARD (Top 5)")
    print("="*40)
    for i, r in enumerate(ranked[:5]):
        eff_str = f"{r['efficiency']:.2f}%" if r['efficiency'] is not None else "N/A"
        print(f"{i+1}. {r['name']} | Eff: {eff_str} | Rg: {r['rg']:.4f} \u00c5")
    print("="*40 + "\n")

if __name__ == "__main__":
    pass
