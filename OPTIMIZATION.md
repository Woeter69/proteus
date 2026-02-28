# Micro-Optimization Sprint

## Objective
Thoroughly scan the codebase one file at a time, analyzing every function individually. 
For each function, we define:
1.  **Functionality:** What does the function do?
2.  **Analysis:** Is there anything redundant, slow, or incorrect?
3.  **Optimization:** Apply surgical improvements to maximize performance and code quality.

## Progress Tracking

### [X] `src/topology.py`
- **Function:** `generate_topology`
    - **Status:** Optimized.
    - **Changes:** Moved imports to top-level, removed dead code (`get_bond_type`), and replaced manual NumPy geometric calculations with optimized RDKit `rdMolTransforms` C++ routines.

### [ ] `src/simulation.py`
- **Status:** Pending.

### [ ] `src/analysis.py`
- **Status:** Pending.

### [ ] `src/visualization.py`
- **Status:** Pending.
