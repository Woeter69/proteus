#!/bin/bash
# Proteus Training Helper Script

# Default values
EPOCHS=100
BATCH_SIZE=16
LR=2e-5
NAME="default"

# Help message
show_help() {
    echo "Usage: ./train.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --epochs N       Number of training epochs (default: $EPOCHS)"
    echo "  --batch N        Batch size (default: $BATCH_SIZE)"
    echo "  --lr N           Learning rate (default: $LR)"
    echo "  --name STR       Name of the run (creates models/NAME/) (default: $NAME)"
    echo "  --help           Show this help message"
    echo ""
    echo "Example:"
    echo "  ./train.sh --name proteus_v1 --epochs 200 --batch 32"
}

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --epochs) EPOCHS="$2"; shift ;;
        --batch) BATCH_SIZE="$2"; shift ;;
        --lr) LR="$2"; shift ;;
        --name) NAME="$2"; shift ;;
        --help) show_help; exit 0 ;;
        *) echo "Unknown parameter: $1"; show_help; exit 1 ;;
    esac
    shift
done

# Run training using the conda environment
echo "[*] Starting Proteus Training Pipeline..."
echo "[*] Run Name: $NAME"
echo "[*] Configuration: Epochs=$EPOCHS, Batch=$BATCH_SIZE, LR=$LR"

conda run --no-capture-output -n proteus_env python -m src.ml.train \
    --epochs "$EPOCHS" \
    --batch_size "$BATCH_SIZE" \
    --lr "$LR" \
    --name "$NAME"
