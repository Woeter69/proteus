"""
Module IV: Analytics
Parses simulation logs to calculate physical metrics (Radius of Gyration) and generates stability plots.
"""

from pathlib import Path
import sys

try:
    import matplotlib
    matplotlib.use('Agg') # Use non-interactive backend
    import matplotlib.pyplot as plt
    import numpy as np
    PLOT_AVAILABLE = True
except ImportError:
    PLOT_AVAILABLE = False

def moving_average(data, window_size):
    if window_size <= 1:
        return data
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

def analyze_results(log_path: Path, output_plot: Path = None):
    """
    Parses the LAMMPS log file to find physical metrics and optionally generates a stability plot.

    Args:
        log_path (Path): Path to the LAMMPS log file.
        output_plot (Path, optional): Path to save the stability plot.
    
    Returns:
        float: The final Rg value.
    """
    print(f"[*] Analyzing results from {log_path}")
    
    if not log_path.exists():
        print(f"Error: Log file {log_path} not found.")
        return None
        
    print(f"[*] Parsing log file: {log_path.name}")

    steps = []
    rg_values = []
    temp_values = []
    energy_values = []
    
    try:
        with open(log_path, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
            
        # Find the LAST run's thermo data section
        last_step_idx = -1
        col_map = {}
        
        for i, line in enumerate(lines):
            if line.strip().startswith("Step") and "c_rg" in line:
                last_step_idx = i
                headers = line.strip().split()
                col_map = {h: idx for idx, h in enumerate(headers)}
        
        if last_step_idx == -1:
             print("Error: Could not find thermo output with c_rg in log.")
             return None

        rg_col = col_map.get("c_rg")
        temp_col = col_map.get("Temp")
        energy_col = col_map.get("PotEng") or col_map.get("E_pair") # Proteus uses epair in thermo_style

        for line in lines[last_step_idx+1:]:
            tokens = line.strip().split()
            if not tokens: continue
            
            if tokens[0].replace('.', '', 1).isdigit():
                try:
                    step = int(tokens[0])
                    rg = float(tokens[rg_col])
                    temp = float(tokens[temp_col])
                    energy = float(tokens[energy_col])
                    
                    steps.append(step)
                    rg_values.append(rg)
                    temp_values.append(temp)
                    energy_values.append(energy)
                except (ValueError, IndexError):
                    pass
            elif "Loop time" in line:
                break
                
        if not rg_values:
            print("Warning: No data extraction.")
            return 0.0
            
        final_rg = rg_values[-1]
        print(f"[*] Final Radius of Gyration (Rg): {final_rg:.4f} Angstroms")

        # Generate stability plot if requested and possible
        if output_plot and steps:
            if not PLOT_AVAILABLE:
                print("[!] Matplotlib/Numpy not found. Skipping stability plot.")
            else:
                print(f"[*] Generating stability plot: {output_plot.name}")
                
                # Smooth data if we have many points
                window = max(1, len(steps) // 50)
                s_steps = steps[window-1:]
                s_temp = moving_average(temp_values, window)
                s_energy = moving_average(energy_values, window)

                fig, ax1 = plt.subplots(figsize=(10, 6))

                color = 'tab:red'
                ax1.set_xlabel('Step')
                ax1.set_ylabel('Temperature (K)', color=color)
                ax1.plot(steps, temp_values, color=color, alpha=0.2) # Raw data faint
                ax1.plot(s_steps, s_temp, color=color, linewidth=2, label='Temp (MA)')
                ax1.tick_params(axis='y', labelcolor=color)

                ax2 = ax1.twinx()
                color = 'tab:blue'
                ax2.set_ylabel('Potential Energy (kcal/mol)', color=color)
                ax2.plot(steps, energy_values, color=color, alpha=0.2) # Raw data faint
                ax2.plot(s_steps, s_energy, color=color, linewidth=2, label='Energy (MA)')
                ax2.tick_params(axis='y', labelcolor=color)

                plt.title(f'System Equilibrium Stability (Smoothing Window: {window})')
                fig.tight_layout()
                plt.savefig(output_plot)
                plt.close()

        return final_rg

    except Exception as e:
        print(f"Error analyzing log file: {e}")
        return None

if __name__ == "__main__":
    pass