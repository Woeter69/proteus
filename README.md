# Proteus: Automated Polymer Nanoprecipitation Pipeline

**Proteus** is a computational research pipeline designed to simulate the self-assembly and nanoprecipitation of polymer chains. By automating the transition from chemical text (SMILES) to physical simulation (Molecular Dynamics), it allows researchers to rapidly screen and analyze polymer folding behaviors.

## Features

*   **Text-to-Structure**: Instantly converts SMILES strings into 3D molecular geometries with explicit hydrogens.
*   **Automated Topology**: Generates LAMMPS-compliant data files with generic force field parameters (Lennard-Jones).
*   **Physics Engine**: Runs implicit solvent simulations using Langevin dynamics with optimized viscosity for realistic molecular drifting.
*   **Built-in Analytics**: Automatically parses simulation logs to calculate the final Radius of Gyration ($R_g$).
*   **Automated Visualization**: Uses **Ovito** to render high-quality, color-coded GIF animations of the simulation.

## Installation

Proteus uses **Conda** to manage dependencies (Python, RDKit, LAMMPS, NumPy).

1.  **Clone the repository:**
    ```bash
    git clone <repository-url>
    cd proteus
    ```

2.  **Create the environment:**
    ```bash
    conda env create -f environment.yml
    ```

3.  **Install Visualization Engine (Optional but Recommended):**
    ```bash
    pip install ovito
    ```

## Usage

The easiest way to run the pipeline is using the `run.sh` helper script.

### Using the Helper Script
```bash
chmod +x run.sh
./run.sh "<SMILES>" "<NAME>" [STEPS]
```
*   **Example (3-unit PEO)**: `./run.sh "CCOCCOCCO" "PEO_3" 10000`
*   **Example (Multi-chain)**: `./run.sh "CCO.CCO.CCO" "Triple_Chain" 10000`

### Arguments
*   `SMILES`: The chemical structure. Use a dot `.` to separate independent molecules.
*   `NAME`: The folder name for your results in `output/`.
*   `STEPS`: (Optional) Simulation duration. Default is 10,000.

## Output

Results are saved in `output/<NAME>/`:

*   `animation.gif`: **High-quality video** of the simulation (color-coded by molecule).
*   `polymer.data`: LAMMPS topology file (3D structure).
*   `simulation.in`: The generated LAMMPS input script.
*   `simulation.log`: Raw simulation data (energies, temperatures).
*   `trajectory.dump`: Raw atom positions (viewable in external tools like VMD/Ovito).

## Architecture

1.  **Topology Architect (`topology.py`)**: SMILES $\to$ 3D Topology via RDKit & UFF.
2.  **Simulation Engine (`simulation.py`)**: Langevin dynamics (300K) with tuned viscosity for smooth motion.
3.  **Analytics (`analysis.py`)**: Log parsing and $R_g$ calculation.
4.  **Visualization (`visualization.py`)**: Headless rendering of trajectories via Ovito.

## Simulation Physics

The simulation uses a **Generic Hydrophobic Interaction** model:
*   **Force Field**: Lennard-Jones (epsilon=0.105, sigma=2.5).
*   **Solvent**: Implicit solvent via Langevin thermostat.
*   **Viscosity**: Tuned damping for realistic drifting rather than high-frequency vibrations.
