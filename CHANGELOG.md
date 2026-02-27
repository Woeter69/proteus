# Proteus Changelog

All notable changes to the Proteus project will be documented in this file.

## [Unreleased]
### Changed
- **Archived AI Module**: The entire machine learning and AI property prediction suite has been moved to `archive/ai/`.
    - Removed `src/ml/` (Dataset, Model, Training, Prediction logic).
    - Removed `models/` (Pre-trained ChemBERTa model weights).
    - Removed `data/train_supplement/` (Training datasets).
    - Removed `train.sh` helper script.
- **CLI Clean-up**: 
    - Removed `--predict` and `--model` flags from `main.py`.
    - Removed `predict` command from `run.sh`.
- **Dependency Optimization**: 
    - Removed `transformers`, `torch`, `pandas`, and `scikit-learn` from `requirements.txt` and `environment.yml` to reduce environment footprint and resolve potential conflicts.
- **Roadmap Update**: Shifted focus back to physical Molecular Dynamics simulations and "The Payload" (Drug Encapsulation) features.
- **Micro-optimization**: 
    - Tuned neighbor list update frequency (`every 10`) achieving a **30% performance speedup** in simulation loop times.
    - Increased minimization iterations to `1000` to ensure better structural stability for complex chains before dynamics.

## [1.4.0] - 2026-02-17
### Changed
- **Project Restructuring**: 
    - Moved all web-related components (Frontend, Backend, Docker, Nginx) into a dedicated `platform/` directory for better separation of concerns between the CLI core and the Web UI.
    - Renamed `web/` to `platform/frontend/` and `backend/` to `platform/backend/`.
    - Consolidated deployment configurations into `platform/` (Dockerfile, docker-compose.yml, nginx.conf).
- **Backend Refactoring**:
    - Updated `sys.path` and base directory calculations in FastAPI and Celery worker to support the new directory depth.
    - Implemented robust root directory detection in `platform/backend/main.py`, `worker.py`, and `database.py`.
- **Deployment**:
    - Updated `Dockerfile` and `vercel.json` to reflect new paths.
    - Updated root `Makefile` and root `.gitignore`.

## [1.3.2] - 2026-02-14
### Added
- **Dihedral Support**:
    - Implemented automatic dihedral detection in `src/topology.py` using RDKit neighbor analysis.
    - Added `dihedral_style harmonic` to LAMMPS input generation in `src/simulation.py`.
    - Differentiates between Single (n=3) and Double (n=2, planar) bonds for improved structural stability.
    - Updated `main.py` to pass dihedral parameters through the pipeline.
### Fixed
- **System Stability**: Resolved excessive gyration and temperature fluctuations in small simulation boxes by constraining bond rotations with proper dihedrals.

## [Unreleased] - 2026-02-15
### Changed
- **Git History**: Restructured recent feature implementation into logical commits across Feb 12 and Feb 14 for better traceability.

## [1.3.1] - 2026-01-27
### Added
- **AI Inference Engine**:
    - Trained a high-accuracy Transformer model (`models/v3`) achieving **90.0% Accuracy** (MAE 1.59 Å) on the NeurIPS 2025 dataset.
    - Updated `predict.py` and `main.py` to support loading custom model versions via `--model`.
    - Enhanced `run.sh` to support `predict <smiles> <model_name>` syntax.
    - Added `train.sh` helper script for parameterized training runs.
- **Training Infrastructure**:
    - Implemented metadata saving (JSON) for training runs to track hyperparameters and metrics.
    - Fixed JSON serialization bug for NumPy types in training logs.

## [1.3.0] - 2026-01-27
### Added
- **Machine Learning Module (Phase 2)**:
    - Implemented `src/ml/` for AI-driven property prediction.
    - **Architecture**: `ChemBERTa` (Transformer) model adapted for regression ($R_g$ prediction).
    - **Infrastructure**: Set up PyTorch with CUDA 12.1 support for RTX 40-series GPUs.
    - **Training Pipeline**: Created `dataset.py`, `model.py`, and `train.py` with SafeTensors support.
    - **Synthetic Data**: Added `synthetic_data.py` to generate dummy polymer datasets for immediate pipeline verification.
- **Environment**: 
    - Resolved dependency conflicts between `lammps` and `pytorch` by implementing a hybrid Conda/Pip installation strategy in `environment.yml`.
    - Added `transformers`, `pandas`, `scikit-learn` to the environment.

## [1.2.1] - 2026-01-27
### Added
- **Web Payload Support**: Integrated the "Payload" (Drug Encapsulation) feature into the Web Interface.
    - Updated `SimulationRequest` API to accept `payload` and `payload_count`.
    - Added optional payload fields to the simulation submission form.
    - Database schema migrated to store payload configuration.
    - Celery worker updated to handle multi-molecule payload injections.

## [1.2.0] - 2026-01-21
### Added
- **Email Notifications**: Implemented automated HTML email alerts for completed or failed simulations using `smtplib` (Backend) and a new user input field (Frontend).
- **PostgreSQL Support**: Added production-grade database support. Users can now choose between SQLite (default) and PostgreSQL via `DATABASE_URL`.
- **Docker Infrastructure**: Added `docker-compose.yml` to spin up PostgreSQL and Redis services.
- **Documentation Site**: Created a comprehensive `/docs` page in the Web Interface, aggregating info from all markdown files into a navigable UI.
- **Environment Config**: Added `.env.example` to document all required environment variables for Email, Database, and App settings.
- **Dynamic Geometry**: Refactored the pipeline to automatically detect and map unique bond lengths and angles from RDKit structures. No more "exploding" molecules due to uniform C-H bond lengths.
- **Physics Refinement**: Integrated `special_bonds` scaling for 1-4 interactions and dynamic LJ cutoffs based on atom size.

### Changed
- **Database Schema**: Added `user_email` column to the `Simulation` table.
- **Environment**: Added `psycopg2-binary` and `fastapi-mail` (placeholder) dependencies to `environment.yml`.
- **Frontend**: Updated Landing Page to link to the new Documentation section.

## [1.1.0] - 2026-01-20
### Added
- **CHONS Forcefield**: Comprehensive support for Carbon, Hydrogen, Oxygen, Nitrogen, and Sulfur using OPLS-AA parameters.
- **Angle Generation**: Automatic detection and harmonic angle coefficient application for more realistic polymer structures.
- **Bond Order Support**: Differentiates between Single, Double, and Triple bonds with specific harmonic constants.
- **The Payload (Drug Encapsulation)**: New `--payload` and `--payload_count` flags to simulate drug molecules within the polymer solvent.
- **Improved Topology**: Switched to individual molecule embedding and random spatial distribution to prevent high-energy overlaps and LAMMPS "missing atom" crashes.
- **Advanced Physics Control**: Exposed `--temp`, `--damp`, `--timestep`, and `--padding` to the CLI for fine-tuned simulation control.

### Changed
- Refactored `src/topology.py` to handle mixed-molecule systems robustly.
- Updated `src/simulation.py` to handle multi-type bond and pair coefficients with arithmetic mixing.
- Enhanced `run.sh` to support both legacy and advanced pass-through command styles.

### Fixed
- Fixed LAMMPS "Missing Atoms" crash caused by high-energy initial overlaps in multi-molecule simulations.
- Fixed string interpolation bug in `simulation.py` for the timestep variable.


### Changed
- **Visuals**: Removed orbiting rings from the main 3D Hero visualization (`IridescentReactor`) for a cleaner look.
- **Visuals**: Removed the inner structure (icosahedron and sphere) from the main 3D Hero cube.
- **Visuals**: Enhanced Iridescence effect to use true thin-film interference physics (soap bubble spectrum) by adjusting `iridescenceIOR` to 1.33 and adding volume/thickness.
- **Visuals**: Removed `Lightformer` elements to prevent white washout and reflection blowouts, ensuring clear rainbow dispersion.
- **Visuals**: Fixed "milky" appearance of the iridescent cube by setting base color to black.
- **Visuals**: Replaced `Stars` component in reflection environment with a custom **Manual Star Field** (scattered bright spheres) to guarantee visibility in reflections.
- **Visuals**: Optimized iridescence (`iridescenceIOR: 1.8`, `range: [100, 400]`) to create a **prismatic chromatic aberration effect** on reflected star light, causing the star reflections to split into rainbow spectrums.
- **Visuals**: Removed refraction (opaque black base) for ultra-sharp star reflections.
- **Visuals**: Added glowy white `<Edges />` to the Hero Cube to enhance visibility and aesthetics.

## [Unreleased]



### Added

- Automated stability graph generation (`stability.png`) after simulations to visualize Temperature and Potential Energy equilibrium.

- Matplotlib dependency for scientific plotting.



### Changed

- Recreated `proteus_env` with updated dependencies.

- Enhanced `src/analysis.py` to parse thermodynamic data for stability reporting.

 - 2026-01-19

### Added
- **Visual Overhaul**: Implemented a "Triple-Nested Tesseract" hero component with:
    - Outer refractive glass shell with rapid rainbow dispersion.
    - Middle wireframe cube with glowing HDR outlines.
    - Inner rotating core cube.
- **Rendering**: Enhanced 3D scene with high-intensity HDR outlines (`Edges`) and optimized `Bloom` settings for a "dark atmospheric" look.

### Added (Previous)
- **Visual Enhancements**:
    - Redesigned 3D Cube Hero and decoration components with **true filled glowy trapeziums** and **filled central squares** on each face.
    - Added **Twin Blade Connectors** (long, glowy diamond blades) that link opposite trapeziums across each face, creating an intersecting cross-hair effect.
    - Preserved high-quality glass `MeshTransmissionMaterial` while adding a structured, high-intensity HUD aesthetic.
- **Web Interface (Phase 1)**:
    - Re-implemented and enhanced the **3D Cube Hero** component with a glass-like `MeshTransmissionMaterial` and glowy white outlines.
    - Integrated `Bloom` post-processing for a futuristic "Black & Glowy" aesthetic.
    - Confirmed landing page introduces the Proteus project with interactive 3D elements.
- **Backend Infrastructure**:
    - Added `celery` and `redis` to `environment.yml` for asynchronous task processing.
    - Configured Celery app in `backend/celery_app.py`.
    - Created `backend/worker.py` to handle simulation pipeline tasks in the background.
    - Exposed `POST /api/simulate` endpoint to trigger simulations.
- **Automation**: Added `conda env update --file environment.yml --prune` to `dev.sh` and `Makefile` to automate dependency management.
- **Monorepo Architecture**: Established a unified repository structure for CLI, Backend, and Frontend.
- **Backend**:
    - Initialized `backend/` directory using **FastAPI** (Python).
    - Configured `uvicorn` as the ASGI server.
    - Added `backend/main.py` with basic health check endpoints.
    - Configured CORS to allow requests from `localhost:3000`.
- **Frontend**:
    - Initialized `web/` directory using **Next.js 14+** (App Router, TypeScript, Tailwind CSS).
    - Pre-installed **Three.js** dependencies (`three`, `@react-three/fiber`, `@react-three/drei`) for future 3D visualization features.
- **Environment Management**:
    - Updated `environment.yml` to include `fastapi`, `uvicorn`, and `python-multipart`.
    - Integrated Conda environment usage into development scripts.
- **Development Tooling**:
    - Created `dev.sh` to launch both Backend and Frontend simultaneously with a single command.
    - Created `Makefile` with a `up` target (`make up`) for easy startup.
    - Updated `.gitignore` to handle `backend/` and `web/` specific artifacts.

### Changed
- **Dependency Management**: Moved backend dependencies from `backend/requirements.txt` to the central `environment.yml` to leverage the existing `proteus_env` Conda environment.
- **Run Scripts**: Updated `dev.sh` to explicitly run the backend using `conda run -n proteus_env` to ensure access to scientific libraries (`rdkit`, `lammps`).
- **Deployment**: Updated `Dockerfile` and `start.sh` to include `redis-server` and `celery worker`, making the Docker image fully functional for simulations out-of-the-box.
- **Dependencies**: Added `libgomp1` to the system dependencies for LAMMPS support in the container.

### Fixed
- **Docker Build**: Fixed a syntax error in `Dockerfile` line 50 where a multi-line `echo` command was incorrectly formatted, causing a parse error.
- **Frontend Build**: Resolved a TypeScript error in `web/app/page.tsx` where `disableNormalPass` was used instead of the correct `enableNormalPass={false}` prop for `EffectComposer`.
- **Frontend Integration**: Updated `web/app/simulation/[id]/page.tsx` to use relative URLs for result files, ensuring compatibility with Nginx proxying.
- **Frontend Logic**: Fixed a `ReferenceError` in `web/app/simulation/page.tsx` where `API_BASE` was not defined within the scope of the `SimulationPage` component.
- **Nginx Config**: Added a `/files` proxy location to `nginx.conf.template` to allow serving simulation results (logs, GIFs) from the backend.
- **Backend**: Fixed undefined references to `Base` and `Simulation` in `backend/main.py` by correctly namespacing them under `models`.