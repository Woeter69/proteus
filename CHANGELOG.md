# Proteus Changelog

All notable changes to the Proteus project will be documented in this file.

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

## [Unreleased] - 2026-01-19

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
