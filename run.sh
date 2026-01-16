#!/bin/bash
# Helper script to run Proteus
# Usage: ./run.sh <SMILES> <NAME> [STEPS]

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <SMILES> <NAME> [STEPS]"
    exit 1
fi

SMILES=$1
NAME=$2
STEPS=${3:-10000} # Default to 10000 if not provided

echo "Running Proteus:"
echo "  Polymer: $SMILES"
echo "  Job Name: $NAME"
echo "  Duration: $STEPS steps"

python main.py --smiles "$SMILES" --name "$NAME" --steps "$STEPS"