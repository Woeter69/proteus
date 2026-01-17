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

### Fixed
- N/A
