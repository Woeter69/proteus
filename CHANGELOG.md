# Proteus Changelog

All notable changes to the Proteus project will be documented in this file.

## [Unreleased] - 2026-01-17

### Added
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
- **Frontend Stability**: Wrapped 3D components (`FancyCube`, `Environment`, `DecorationCube`) in `Suspense` boundaries to prevent the 3D Canvas from disappearing during navigation or resource loading.
- **Nginx Config**: Added a `/files` proxy location to `nginx.conf.template` to allow serving simulation results (logs, GIFs) from the backend.
- **Backend**: Fixed undefined references to `Base` and `Simulation` in `backend/main.py` by correctly namespacing them under `models`.
