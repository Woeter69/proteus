#!/bin/bash
# Helper script to run Proteus
# Usage: ./run.sh <SMILES> <NAME>

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <SMILES> <NAME>"
    exit 1
fi

SMILES=$1
NAME=$2

echo "Running Proteus for $NAME with SMILES: $SMILES"
python main.py --smiles "$SMILES" --name "$NAME"
