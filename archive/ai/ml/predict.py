import torch
from transformers import AutoTokenizer
from .model import PolymerPredictor
import argparse
import sys
import os

def predict_rg(smiles: str, model_name: str = "default"):
    """
    Predicts the Radius of Gyration (Rg) for a given SMILES string using the trained Transformer.
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Check for direct path or named model in 'models/' dir
    if os.path.isfile(model_name):
        model_path = model_name
    else:
        model_path = os.path.join("models", model_name, "model.pt")
        
    # Fallback to root if 'default' fails
    if not os.path.exists(model_path) and model_name == "default":
        model_path = "proteus_model.pt"

    # 1. Initialize Model Architecture
    try:
        model = PolymerPredictor().to(device)
    except Exception as e:
        print(f"Error initializing model: {e}")
        return None

    # 2. Load Trained Weights
    if not os.path.exists(model_path):
        print(f"Error: Model file '{model_path}' not found. Have you trained it yet?")
        return None
        
    try:
        # Load state dict (map_location handles CPU/GPU mismatch automatically)
        state_dict = torch.load(model_path, map_location=device, weights_only=True)
        model.load_state_dict(state_dict)
        model.eval()
    except Exception as e:
        print(f"Error loading model weights: {e}")
        return None

    # 3. Tokenize Input
    tokenizer = AutoTokenizer.from_pretrained("seyonec/ChemBERTa-zinc-base-v1")
    
    # Preprocessing (same as training)
    clean_smiles = smiles.replace("*", "C") 
    
    inputs = tokenizer(
        clean_smiles,
        return_tensors="pt",
        padding="max_length",
        truncation=True,
        max_length=128
    )
    
    input_ids = inputs['input_ids'].to(device)
    attention_mask = inputs['attention_mask'].to(device)

    # 4. Inference
    with torch.no_grad():
        prediction = model(input_ids, attention_mask)
        rg_value = prediction.item()
        
    return rg_value

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict Polymer Rg using AI")
    parser.add_argument("smiles", type=str, help="SMILES string of the polymer")
    parser.add_argument("--model", type=str, default="proteus_model.pt", help="Path to trained model")
    
    args = parser.parse_args()
    
    print(f"[*] Analyzing: {args.smiles}")
    rg = predict_rg(args.smiles, args.model)
    
    if rg is not None:
        print(f"[*] Predicted Rg: {rg:.4f} Angstroms")
    else:
        sys.exit(1)
