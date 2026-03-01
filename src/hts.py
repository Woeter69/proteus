"""
Module VII: High-Throughput Screening (HTS)
Batch processes multiple simulations from a CSV file and ranks them.
"""

import csv
import sys
from pathlib import Path
from copy import deepcopy

def run_batch(csv_path: str, global_args):
    """
    Reads a CSV and runs the pipeline for each row.
    CSV Format: name,smiles,count,payload,payload_count,steps
    """
    from main import run_pipeline
    
    csv_file = Path(csv_path)
    if not csv_file.exists():
        print(f"Error: Batch CSV file not found: {csv_path}")
        sys.exit(1)
        
    print(f"[*] Starting High-Throughput Screening (HTS): {csv_path}")
    
    results = []
    
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
                    continue
                
                try:
                    result = run_pipeline(run_args)
                    results.append(result)
                except Exception as e:
                    print(f"[!] Run failed for {run_args.name}: {e}")
                    continue
                    
        if not results:
            print("No simulations completed successfully in batch.")
            return

        # Rank Results
        rank_and_summarize(results, csv_file.parent)
        
    except Exception as e:
        print(f"Error during HTS batch run: {e}")

def rank_and_summarize(results, output_dir: Path):
    """
    Sorts results by Rg (folding stability) and encapsulation efficiency.
    Generates a summary CSV.
    """
    summary_path = output_dir / "hts_ranking_summary.csv"
    
    print(f"[*] Ranking {len(results)} results...")
    
    # Ranking logic: 
    # 1. Primary: Encapsulation Efficiency (Descending)
    # 2. Secondary: Radius of Gyration (Ascending - tighter folding is often better)
    
    # Sort: Efficiency (highest first), then Rg (smallest first)
    # Note: Efficiency might be None, so we treat None as 0.0
    ranked = sorted(
        results, 
        key=lambda x: (-(x.get('efficiency') or 0.0), x['rg'])
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
