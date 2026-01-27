#!/bin/bash
# Helper script to run Proteus
# Usage Legacy: ./run.sh <SMILES> <NAME> [STEPS] [COUNT] [EXTRA_ARGS...]
# Usage Advanced: ./run.sh --smiles <SMILES> --name <NAME> ...

if [ "$#" -eq 0 ]; then
    echo "Usage Legacy: $0 <SMILES> <NAME> [STEPS] [COUNT]"
    echo "Usage Advanced: $0 --smiles <SMILES> ... (Pass any main.py flags)"
    echo "New Variables available: --temp, --damp, --epsilon, --sigma, --timestep, --padding"
    echo "Payload Variables: --payload <SMILES>, --payload_count <N>"
    exit 1
fi

# Check if first argument is a flag (Advanced Mode)
if [[ "$1" == "predict" ]]; then
    echo "Running Proteus AI Inference..."
    shift
    # Usage: ./run.sh predict <SMILES> [MODEL_NAME]
    if [ -z "$2" ]; then
         conda run --no-capture-output -n proteus_env python main.py --predict --smiles "$1"
    else
         conda run --no-capture-output -n proteus_env python main.py --predict --smiles "$1" --model "$2"
    fi
elif [[ "$1" == --* ]]; then
    echo "Running in Advanced Mode (Pass-through)..."
    conda run --no-capture-output -n proteus_env python main.py "$@"
else
    # Legacy Mode with optional extras
    SMILES=$1
    NAME=$2
    
    # Check if $3 is a number (STEPS)
    if [[ "$3" =~ ^[0-9]+$ ]]; then
        STEPS=$3
        shift 3
        # Check if $1 (originally $4) is a number (COUNT)
        if [[ "$1" =~ ^[0-9]+$ ]]; then
            COUNT=$1
            shift
        else
            COUNT=1
        fi
    else
        STEPS=10000
        COUNT=1
        shift 2
    fi
    
    # Remaining arguments in $@ are passed as extra flags
    echo "Running Proteus (Legacy Mode):"
    echo "  Polymer: $SMILES"
    echo "  Job Name: $NAME"
    echo "  Duration: $STEPS steps"
    echo "  Copies:  $COUNT"
    echo "  Extra Args: $@"

    conda run --no-capture-output -n proteus_env python main.py \
        --smiles "$SMILES" \
        --name "$NAME" \
        --steps "$STEPS" \
        --count "$COUNT" \
        --render \
        "$@"
fi
