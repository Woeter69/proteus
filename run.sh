#!/bin/bash
# Helper script to run Proteus
# Usage: ./run.sh <SMILES> <NAME> [STEPS] [COUNT]

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <SMILES> <NAME> [STEPS] [COUNT]"
    exit 1
fi

SMILES=$1
NAME=$2
STEPS=${3:-10000} # Default 10k steps
COUNT=${4:-1}     # Default 1 molecule

echo "Running Proteus:"
echo "  Polymer: $SMILES"
echo "  Copies:  $COUNT"
echo "  Job Name: $NAME"
echo "  Duration: $STEPS steps"

# Use 'conda run' which works in scripts without needing 'conda activate' or 'conda init'
conda run --no-capture-output -n proteus_env python main.py --smiles "$SMILES" --name "$NAME" --steps "$STEPS" --count "$COUNT" --render
