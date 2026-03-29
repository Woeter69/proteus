# Proteus Documentation

Welcome to the developer documentation for **Proteus**, the automated polymer nanoprecipitation pipeline.

## Core Modules

- **Topology Architect**: Converts SMILES to LAMMPS data files.
- **Simulation Engine**: Configures and runs Molecular Dynamics via LAMMPS.
- **Analytics**: Parses logs and computes physical metrics like $R_g$.
- **HTS**: Orchestrates high-throughput screening from CSV datasets.
- **Visualization**: Headless rendering of trajectories via Ovito.
- **Report**: Generates professional PDF lab notebooks.

## Getting Started

To install dependencies and set up the environment:

```bash
make setup
```

To run tests:

```bash
make test
```
