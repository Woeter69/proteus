# Proteus: Project History & Overview

## What is Proteus?
**Proteus** is an automated computational research pipeline designed for Materials Science and Pharmacology. It simulates the self-assembly and nanoprecipitation of polymer chains. By automating the transition from simple chemical text codes (SMILES) to complex physical Molecular Dynamics (MD) simulations, it allows researchers to rapidly screen and analyze how different polymers fold and interact with other molecules.

A major focus of the project is **"The Payload" (Drug Encapsulation)**—simulating how a polymer "suitcase" (nanoparticle) wraps around and traps smaller drug molecules for targeted medical delivery.

## What Help Does It Bring?
Proteus acts as a **virtual crash test for medicine and materials**. 
In the real world, synthesizing and testing hundreds of polymer variations in a lab is incredibly expensive and time-consuming. Proteus solves this by allowing researchers to:
1.  **Screen rapidly:** Test thousands of polymer/drug combinations on a computer for pennies.
2.  **Measure Efficiency:** Automatically calculate the "Encapsulation Efficiency" (what percentage of the drug actually stayed inside the polymer).
3.  **Visualize:** Generate high-quality 3D animations and cross-sections of the molecular folding process to see exactly *why* a design succeeded or failed.
4.  **Save Resources:** Identify the top 2 or 3 candidates virtually, ensuring physical lab resources are only spent on the most promising designs.

## How We Achieved That (System Architecture)
Proteus is built on a highly optimized, modular pipeline:
*   **Module I - Topology Architect (`src/topology.py`):** Uses **RDKit** and the Universal Force Field (UFF) to convert 1D SMILES strings into 3D geometries. It automatically detects bonds, angles, and proper dihedrals, assigning physical properties using a CHONS (Carbon, Hydrogen, Oxygen, Nitrogen, Sulfur) OPLS-AA force field.
*   **Module II - Simulation Engine (`src/simulation.py`):** Generates and executes **LAMMPS** simulation scripts. It uses Langevin dynamics to simulate an implicit solvent (like water or ethanol) with specifically tuned viscosity/damping to create realistic molecular drifting and folding.
*   **Module III - Analytics (`src/analysis.py`):** Features a highly memory-efficient parser that reads gigabyte-sized simulation logs line-by-line. It calculates thermodynamic stability, the final size of the nanoparticle (Radius of Gyration, $R_g$), and the exact geometric Encapsulation Efficiency of payload molecules.
*   **Module IV - Visualization (`src/visualization.py`):** Uses **Ovito** (Open Visualization Tool) via Python to automatically render color-coded, GPU-accelerated 3D GIFs and cross-sectional views of the simulation.

---

## Comprehensive Project History (Start to Finish)

### Phase 1: Inception & Core Pipeline
*   **Added:** Project scaffolding, basic SMILES-to-LAMMPS pipeline.
*   **Added:** Core modules: `topology.py`, `simulation.py`, `analysis.py`, and `main.py` controller.
*   **Added:** Basic visualization via Ovito for rendering animations.
*   **Added:** `run.sh` helper script for easy CLI usage.

### Phase 2: The Web Era (Later Removed)
*   *Note: Proteus was briefly developed as a full-stack web application before being split.*
*   **Added:** A Next.js frontend featuring a highly complex "Triple-Nested Tesseract" 3D glowing hero component using Three.js and WebGL.
*   **Added:** FastAPI backend, Celery workers, Redis, and PostgreSQL via Docker to handle asynchronous simulation requests from the web.
*   **Added:** Email notification system for completed simulations.
*   **Removed:** **All web and backend components were entirely migrated out of this repository** to a dedicated platform repo to keep this CLI tool lightweight and focused strictly on the physics engine.

### Phase 3: The AI / Machine Learning Era (Archived)
*   **Added:** An AI inference engine (`src/ml/`) designed to instantly predict polymer properties without running the full physics simulation.
*   **Added:** A PyTorch-based ChemBERTa (Transformer) model.
*   **Achieved:** Trained the model on synthetic datasets to predict $R_g$ with 90% accuracy.
*   **Removed/Archived:** The entire AI suite (`src/ml/`, weights, and datasets) was moved to `archive/ai/`. The project's philosophy pivoted back to prioritizing high-fidelity physical Molecular Dynamics over AI estimations.

### Phase 4: Advanced Physics & "The Payload"
*   **Added:** **"The Payload" (Drug Encapsulation)** feature. The CLI now accepts `--payload` and `--payload_count` to inject secondary drug molecules into the polymer solvent.
*   **Added:** Comprehensive CHONS OPLS-AA forcefield support.
*   **Added:** Dynamic bond order detection (Single, Double, Triple) and automatic dihedral (torsion) detection using RDKit for realistic polymer stiffness.
*   **Added:** `VARIABLES.md` established as the Single Source of Truth for all physical parameters.
*   **UX Philosophy Enforced:** All complex features, especially the Payload, were made strictly optional to maintain a clean "UX First" user experience.

### Phase 5: Automation & High-Throughput Screening (HTS)
*   **Added:** `--hts` flag for processing CSV spreadsheets containing hundreds of polymer/payload combinations.
*   **Added:** Automatic ranking systems that sort HTS results by structural stability ($R_g$) or Encapsulation Efficiency.
*   **Added:** Automated PDF lab report generation including 2D chemical structures, stability plots, and analysis results.
*   **Added:** `--cut` flag to render 3D cross-sections of the nanoparticles to see the drug trapped inside.
*   **Added:** Solvent profiles (`--solvent water/ethanol`) that automatically configure physical simulation parameters.

### Phase 6: Micro-Optimization Sprint (The Final Polish)
*   **Optimized (`topology.py`):** Replaced manual math with C++ optimized RDKit `rdMolTransforms` routines.
*   **Optimized (`simulation.py`):** Implemented streaming `stdout` via `subprocess.Popen` to prevent memory bloat during long LAMMPS runs. Added GPU scaling.
*   **Optimized (`analysis.py`):** Rewrote the trajectory parser to seek from the end of the file and read line-by-line, preventing memory crashes when analyzing multi-gigabyte files.
*   **Optimized (`visualization.py`):** Added element-specific atomic radii and dynamic frame sampling to drastically reduce GIF rendering times while improving visual quality.
*   **Added:** A complete unit testing suite using `pytest`.