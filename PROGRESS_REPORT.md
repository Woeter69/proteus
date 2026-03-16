# Proteus: A High-Performance Molecular Dynamics Pipeline

**Version**: 1.5.0 (Post-Optimization Sprint)
**Date**: March 11, 2026
**Lead Maintainer**: [Woeter69](https://github.com/Woeter69)

---

## 1. Project Mission & Philosophy

Proteus is a command-line-driven scientific computing pipeline designed to automate the process of running Molecular Dynamics (MD) simulations, starting from a simple SMILES string. Its core mission is to provide researchers with a powerful, flexible, and efficient tool to explore the conformational space of polymers and other molecules without the manual overhead typically associated with preparing, running, and analyzing simulations.

The project adheres to a **CLI-First** philosophy. This strategic choice prioritizes performance, scriptability, and seamless integration into larger computational research workflows, such as those found on High-Performance Computing (HPC) clusters. By focusing on a robust and powerful command-line interface, Proteus ensures maximum utility for its core scientific audience.

---

## 2. System Architecture

The architecture of Proteus is modular and streamlined, reflecting its focused mission. Key decisions have been made to concentrate development efforts on the core simulation engine.

### 2.1. Architectural Evolution: The CLI Focus

In a recent major strategic pivot, all web platform components (including a Next.js frontend and FastAPI backend) were removed from the repository. This decision was made to:
-   **Streamline Maintenance**: Reduce the complexity and overhead of maintaining a web stack.
-   **Enhance Focus**: Concentrate all development resources on the core scientific computing engine.
-   **Promote Integration**: A CLI is far more straightforward to integrate into automated, multi-step research scripts and batch-job schedulers (e.g., Slurm, PBS).

The project's previous exploration into AI/ML-based property prediction has also been formally archived, solidifying the commitment to high-fidelity, physics-based simulation.

### 2.2. Codebase Structure

The heart of the project resides in the `src/` directory. Each module has a distinct and clear responsibility:

-   `main.py`: The main entry point for the CLI. It uses `argparse` to handle user input, validate arguments, and orchestrate the overall simulation pipeline.

-   `topology.py`: Responsible for converting a SMILES string into a 3D molecular structure with the correct topology for simulation. It uses RDKit for initial structure generation and contains all the logic for assigning atom types, bonds, angles, and dihedrals.

-   `simulation.py`: The core simulation engine. This module takes the generated topology, creates input files for the LAMMPS simulator, executes the simulation as a subprocess, and monitors its progress. It has been heavily optimized for performance.

-   `analysis.py`: Handles the post-processing of simulation output. It can parse massive trajectory (`.dump`) and log (`.log`) files efficiently to calculate key physical properties like the Radius of Gyration (Rg).

-   `visualization.py`: Creates high-quality animated GIFs of the simulation trajectory using Ovito. This module is designed to produce clear, informative, and publication-ready visuals of the molecular dynamics.

-   `report.py`: A new module that automates the creation of comprehensive PDF lab reports, summarizing the parameters and results of a simulation.

-   `hts.py`: The High-Throughput Screening engine. This module manages the execution of multiple simulations in a batch, processes the results, and ranks them according to user-defined criteria.

---

## 3. Core Engine Features & Optimizations

Significant effort has been invested in making the simulation engine not only physically accurate but also highly performant.

### 3.1. Performance & GPU Acceleration

A "micro-optimization sprint" yielded substantial performance gains:

-   **30% Simulation Speedup**: Achieved by intelligently tuning the LAMMPS neighbor list update frequency (`neighbor 2.0 bin`, `every 10 delay 0 check yes`), reducing redundant calculations in the main simulation loop.
-   **GPU Support**: The CLI now supports offloading the simulation workload to a GPU, which can provide an order-of-magnitude speedup for large systems. This is activated via a command-line flag.
-   **Efficient I/O**: The simulation module (`src/simulation.py`) now uses `subprocess.Popen` to stream `stdout` in real-time. This prevents memory bloat that could previously occur from buffering gigabytes of log data.
-   **Memory-Conscious Analysis**: The analysis module (`src/analysis.py`) was refactored to read large trajectory files line-by-line, allowing it to process terabyte-scale data without high memory consumption.

### 3.2. Advanced Physics Modeling

-   **Dihedral Support**: The automatic detection and application of harmonic dihedrals provide the necessary rotational constraints to prevent molecules from adopting unrealistic, high-energy conformations. This was a critical step in improving simulation stability.
-   **CHONS Forcefield**: Full support for Carbon, Hydrogen, Oxygen, Nitrogen, and Sulfur using OPLS-AA parameters allows for the simulation of a wide range of organic and biological molecules.
-   **The Payload Feature**: The engine can simulate a system containing a primary polymer "solvent" and a secondary "payload" molecule (e.g., a drug). This is essential for studying drug delivery, encapsulation, and other multi-component systems.

---

## 4. Research Automation Features

Two powerful new features have been added to automate key parts of the research workflow: HTS and PDF reporting.

### 4.1. High-Throughput Screening (HTS) Workflow

The HTS module transforms Proteus from a single-molecule tool into a materials discovery engine. The workflow is simple yet powerful:

1.  **Create an Input File**: Prepare a simple CSV file (e.g., `molecules.csv`) listing the molecules to be screened.

    ```csv
    name,smiles
    PEO,CCOCCOCCO
    Polystyrene,CC(C1=CC=CC=C1)C
    PVA,C(C(O))C
    ```

2.  **Execute the HTS Command**: Run the main script with the `hts` subcommand, pointing to your input file and specifying a ranking metric.

    ```bash
    # Run simulations for all molecules and rank by final Radius of Gyration (Rg)
    python main.py hts --input molecules.csv --rank-by rg
    ```

3.  **Analyze the Results**: The `hts` module creates a summary file, `hts_ranking_summary.csv`, with the calculated properties for each molecule, sorted by your chosen metric.

    ```csv
    name,rg,efficiency,smiles,output_dir
    PEO,21.01,80.0,CCOCCOCCO,/path/to/output/PEO
    Polystyrene,18.54,75.0,CC(C1=CC=CC=C1)C,/path/to/output/Polystyrene
    PVA,15.32,90.0,C(C(O))C,/path/to/output/PVA
    ```

### 4.2. Automated PDF Lab Reports

To streamline documentation and reporting, Proteus can now automatically generate a PDF summary of any simulation.

-   **Content**: The report includes a title page, a table of all simulation parameters (temperature, pressure, etc.), plots of energy and temperature over time to show system stability, and a 2D image of the molecule's final structure.
-   **Trigger**: This feature is integrated directly into the main simulation pipeline and can be enabled with a flag.

---

## 5. Usage Examples

### 5.1. Running a Single Simulation with GPU

```bash
# Run a simulation for a PEO oligomer at 300K, with 10 chains,
# enabling GPU acceleration and generating a final report.

python main.py --smiles "CCOCCOCCO" \
               --name "PEO_Oligomer_300K" \
               --chains 10 \
               --temp 300 \
               --gpu 1 \
               --report
```

### 5.2. Running a Batch Screening (HTS)

```bash
# Execute a high-throughput screening run from a list of SMILES
# in `candidate_polymers.csv` and rank the output by the lowest Rg.

python main.py hts --input candidate_polymers.csv --rank-by rg --sort-order asc
```

---

## 6. Conclusion & Future Outlook

Proteus has successfully evolved into a mature, high-performance, command-line-first tool for molecular dynamics simulation. By strategically focusing on its core engine and shedding non-essential components, the project has delivered powerful new features for research automation. The addition of High-Throughput Screening and automated PDF reporting transforms Proteus from a simple simulator to a genuine research pipeline.

Future development will focus on expanding the capabilities of the core engine, with a particular emphasis on enhancing the **"Payload"** feature to support more complex multi-component systems and interactions.
