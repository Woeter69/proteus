"""
Module IV: Analytics
Parses simulation logs to calculate physical metrics (Radius of Gyration).
"""

from pathlib import Path
import sys

def analyze_results(log_path: Path):
    """
    Parses the LAMMPS log file to find the final Radius of Gyration (Rg).

    Args:
        log_path (Path): Path to the LAMMPS log file.
    
    Returns:
        float: The final Rg value.
    """
    print(f"[*] Analyzing results from {log_path}")
    
    if not log_path.exists():
        print(f"Error: Log file {log_path} not found.")
        return None
        
    print(f"[*] Parsing log file: {log_path.name}")

    rg_values = []
    
    try:
        with open(log_path, 'r') as f:
            lines = f.readlines()
            
        # Find the thermo data section
        # Look for the line starting with "Step"
        start_idx = -1
        for i, line in enumerate(lines):
            if line.strip().startswith("Step"):
                # verify it has c_rg
                if "c_rg" in line:
                    start_idx = i
                    # Don't break immediately, there might be multiple runs (minimization + run)
                    # We want the LAST run.
        
        if start_idx == -1:
            print("Error: Could not find thermo output with c_rg in log.")
            return None
            
        # Parse from start_idx + 1 until "Loop time" or end
        # Since we might have found an early "Step" line, let's look for the data chunks.
        # A robust way is to read all lines that look like numbers after the last "Step" header.
        
        # Re-scan for the LAST "Step" header
        last_step_idx = -1
        col_map = {}
        
        for i, line in enumerate(lines):
            if line.strip().startswith("Step") and "c_rg" in line:
                last_step_idx = i
                headers = line.strip().split()
                for col_i, h in enumerate(headers):
                    col_map[h] = col_i
        
        if last_step_idx == -1:
             print("Error: Could not find thermo output with c_rg in log.")
             return None

        rg_col = col_map.get("c_rg")
        
        for line in lines[last_step_idx+1:]:
            tokens = line.strip().split()
            if not tokens: continue
            
            # Check if line is numeric (start of data)
            if tokens[0].replace('.', '', 1).isdigit():
                try:
                    rg = float(tokens[rg_col])
                    rg_values.append(rg)
                except (ValueError, IndexError):
                    pass
            elif "Loop time" in line:
                break
                
        if not rg_values:
            print("Warning: No Rg data extraction.")
            return 0.0
            
        final_rg = rg_values[-1]
        print(f"[*] Final Radius of Gyration (Rg): {final_rg:.4f} Angstroms")
        return final_rg

    except Exception as e:
        print(f"Error analyzing log file: {e}")
        return None

if __name__ == "__main__":
    pass
