# Proteus Simulation Variables

This document serves as the single source of truth for all configurable variables in the Proteus pipeline.
**Rule:** If a new variable is added to the codebase, it MUST be documented here immediately.

## Input & Composition
| Flag | Variable | Default | Description |
| :--- | :--- | :--- | :--- |
| `--smiles` | `SMILES` | **Required** | The SMILES string of the primary polymer monomer. |
| `--count` | `count` | `1` | Number of polymer monomers to simulate (Aggregation mode). |
| `--payload` | `payload` | `None` | SMILES string of the secondary "drug" or payload molecule. |
| `--payload_count` | `payload_count` | `1` | Number of payload molecules to inject. |
| `--name` | `name` | `simulation` | Name of the job (creates `output/<name>` directory). |

## Physics & Environment
| Flag | Variable | Default | Description |
| :--- | :--- | :--- | :--- |
| `--temp` | `temp` | `300.0` | **Temperature (K)**. Controls thermal energy in the Langevin thermostat. |
| `--damp` | `damp` | `20.0` | **Damping (fs)**. Inverse of viscosity. Lower value = Higher Friction (Thicker Solvent). |
| `--timestep` | `timestep` | `1.0` | **Time Step (fs)**. Resolution of the simulation integration. |
| `--padding` | `padding` | `20.0` | **Padding (Å)**. Extra space around molecules to determine Simulation Box size. |
| `--steps` | `steps` | `10000` | **Duration**. Total number of simulation steps to run. |

## Force Field (Lennard-Jones)
*Note: CHONS atoms (C, H, O, N, S) use specific OPLS-AA parameters by default. These flags act as global overrides or defaults for unknown types.*

| Flag | Variable | Default | Description |
| :--- | :--- | :--- | :--- |
| `--epsilon` | `epsilon` | `0.105` | **Interaction Strength**. Depth of the potential well. |
| `--sigma` | `sigma` | `2.5` | **Particle Size (Å)**. Distance where potential is zero. |

## Visualization
| Flag | Variable | Default | Description |
| :--- | :--- | :--- | :--- |
| `--render` | `render` | `False` | If set, generates an animated GIF of the trajectory using Ovito. |
