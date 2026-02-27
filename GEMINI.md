# Proteus: Gemini Project Memories

## Core UX Principles
- **Optionality is Key**: New features, especially specialized ones like "The Payload" (Drug Encapsulation), must be strictly optional. 
- **User Control**: Never force complex simulation parameters or extra molecules on the user unless they explicitly request them via flags or configuration.
- **UX First**: Always prioritize a clean, understandable interface and output over feature bloating.

## Project Evolution
- Initial implementation focused on automated SMILES-to-MD pipeline.
- Added Ovito integration for high-quality, color-coded GIF visualization.
- **Archived AI**: The ML-based property prediction (CamemBERT/ChemBERTa) has been moved to `archive/ai/`. The project now focuses on simulation accuracy and optional "Payload" features.

## Documentation Rules
- **Variables**: `VARIABLES.md` is the Single Source of Truth for simulation parameters. If a new CLI argument or physics variable is added, it **MUST** be added to `VARIABLES.md` in the same turn.
