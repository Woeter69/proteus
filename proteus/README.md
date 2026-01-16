# Proteus: Automated Polymer Nanoprecipitation Pipeline

**Proteus** is a computational research pipeline designed to simulate the self-assembly and nanoprecipitation of polymer chains. By automating the transition from chemical text (SMILES) to physical simulation (Molecular Dynamics), it allows researchers to rapidly screen and analyze polymer folding behaviors.

## 🚀 Features

*   **Text-to-Structure**: Instantly converts SMILES strings into 3D molecular geometries with explicit hydrogens.
*   **Automated Topology**: Generates LAMMPS-compliant data files with generic force field parameters (Lennard-Jones).
*   **Physics Engine**: Runs implicit solvent simulations using Langevin dynamics to model polymer collapse.
*   **Built-in Analytics**: Automatically parses simulation logs to calculate the final Radius of Gyration ($R_g$).
*   **Modular Design**: Clean separation of topology generation, simulation, and analysis.

## 🛠️ Installation

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

3.  **Activate the environment:**
    ```bash
    conda activate proteus_env
    ```

## 💻 Usage

Run the pipeline using the central `main.py` controller. You can specify the polymer structure, simulation name, and duration.

### Basic Command
```bash
python main.py --smiles "CCOCCOCCO" --name "PEO_Oligomer" --steps 10000
```

### Arguments
*   `--smiles`: (Required) The SMILES string of the polymer (e.g., `CCO` for PEG monomer).
*   `--name`: (Optional) Name of the simulation run. Output will be saved to `output/<name>/`. Default: `simulation`.
*   `--steps`: (Optional) Number of simulation steps to run. Default: `10000`.

## 📂 Output

Results are saved in the `output/` directory, organized by simulation name.

```text
output/
└── PEO_Oligomer/
    ├── polymer.data       # LAMMPS topology file (3D structure)
    ├── simulation.in      # LAMMPS input script
    ├── simulation.log     # Raw simulation logs
    └── trajectory.dump    # Atom trajectories (view in VMD/Ovito)
```

## 🏗️ Architecture

The codebase is split into four specialized modules within `src/`:

1.  **Topology Architect (`topology.py`)**:
    *   Uses **RDKit** to embed SMILES strings into 3D space.
    *   Applies UFF optimization to relax geometry.
    *   Assigns generic masses and atom types.

2.  **Simulation Engine (`simulation.py`)**:
    *   Constructs the LAMMPS input file.
    *   Configures Langevin thermostat (300K) and generic hydrophobic interactions.
    *   Executes `lmp` (LAMMPS) as a subprocess.

3.  **Analytics (`analysis.py`)**:
    *   Parses thermodynamic output from logs.
    *   Extracts the Radius of Gyration ($R_g$) to quantify compactness.

4.  **Controller (`main.py`)**:
    *   CLI entry point that orchestrates the data flow between modules.

## ⚠️ Limitations

*   **Force Field**: Uses a generic "one-size-fits-all" Lennard-Jones potential. For high-accuracy chemical properties, a specific force field (like OPLS-AA or CHARMM) implementation is required.
*   **Solvent**: Uses an implicit solvent model (Langevin dynamics), which approximates solvent effects via friction and random force, rather than explicit water molecules.
