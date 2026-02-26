import torch
from torch.utils.data import DataLoader
from .dataset import PolymerDataset
from .model import PolymerPredictor
from torch.optim import AdamW
import torch.nn as nn
from tqdm import tqdm
import os
import argparse
import numpy as np
from sklearn.metrics import r2_score

import json
from datetime import datetime

def train_model(csv_path, epochs=100, batch_size=16, lr=2e-5, model_name="default"):
    # 1. Setup Paths
    output_dir = os.path.join("models", model_name)
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, "model.pt")
    log_path = os.path.join(output_dir, "training_log.txt")
    
    print(f"[*] Output Directory: {output_dir}")

    # 1. Setup Device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"[*] Training on: {device}")

    # 2. Load Data
    print(f"[*] Loading Dataset: {csv_path}")
    dataset = PolymerDataset(csv_path)
    train_size = int(0.8 * len(dataset))
    val_size = len(dataset) - train_size
    train_dataset, val_dataset = torch.utils.data.random_split(dataset, [train_size, val_size])

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size)

    # 3. Initialize Model
    model = PolymerPredictor().to(device)
    optimizer = AdamW(model.parameters(), lr=lr)
    criterion = nn.MSELoss() 

    # 4. Training Loop
    best_loss = float('inf')
    best_metrics = {}
    
    for epoch in range(epochs):
        model.train()
        total_loss = 0
        
        loop = tqdm(train_loader, desc=f"Epoch {epoch+1}/{epochs}")
        for batch in loop:
            input_ids = batch['input_ids'].to(device)
            mask = batch['attention_mask'].to(device)
            labels = batch['labels'].to(device)

            optimizer.zero_grad()
            
            outputs = model(input_ids, mask)
            loss = criterion(outputs.squeeze(), labels)
            
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
            loop.set_postfix(loss=loss.item())

        avg_loss = total_loss / len(train_loader)
        
        # Validation with Accuracy Metrics
        val_loss, mae, mape, r2 = evaluate(model, val_loader, criterion, device)
        
        print(f"Epoch {epoch+1} | Loss: {val_loss:.2f} | MAE: {mae:.2f} | Accuracy: {100-mape:.1f}% | R2: {r2:.2f}")

        if val_loss < best_loss:
            best_loss = val_loss
            torch.save(model.state_dict(), save_path)
            best_metrics = {
                "epoch": int(epoch + 1),
                "loss": float(val_loss),
                "mae": float(mae),
                "accuracy_percent": float(100 - mape),
                "r2": float(r2)
            }
            print(f" -> Best Model Saved to {save_path}!")

    # Save Metadata
    metadata = {
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "config": {
            "epochs": int(epochs),
            "batch_size": int(batch_size),
            "lr": float(lr)
        },
        "best_metrics": best_metrics
    }
    
    with open(os.path.join(output_dir, "metadata.json"), "w") as f:
        json.dump(metadata, f, indent=4)
    
    print("[*] Training Complete. Metadata saved.")

def evaluate(model, loader, criterion, device):
    model.eval()
    total_loss = 0
    all_preds = []
    all_labels = []
    
    with torch.no_grad():
        for batch in loader:
            input_ids = batch['input_ids'].to(device)
            mask = batch['attention_mask'].to(device)
            labels = batch['labels'].to(device)

            outputs = model(input_ids, mask).squeeze()
            loss = criterion(outputs, labels)
            total_loss += loss.item()
            
            all_preds.extend(outputs.cpu().numpy())
            all_labels.extend(labels.cpu().numpy())
            
    # Metrics
    all_preds = np.array(all_preds)
    all_labels = np.array(all_labels)
    
    mae = np.mean(np.abs(all_preds - all_labels))
    
    # Avoid division by zero for MAPE
    with np.errstate(divide='ignore', invalid='ignore'):
        mape = np.mean(np.abs((all_labels - all_preds) / all_labels)) * 100
        if np.isnan(mape): mape = 100.0
        
    r2 = r2_score(all_labels, all_preds)

    return total_loss / len(loader), mae, mape, r2

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=100, help="Number of training epochs")
    parser.add_argument("--batch_size", type=int, default=16, help="Batch size")
    parser.add_argument("--lr", type=float, default=2e-5, help="Learning rate")
    parser.add_argument("--name", type=str, default="default", help="Name of the model run (creates folder)")
    args = parser.parse_args()

    # If run directly, look for data
    if os.path.exists("data/train.csv"):
        train_model("data/train.csv", epochs=args.epochs, batch_size=args.batch_size, lr=args.lr, model_name=args.name)
    elif os.path.exists("data/training_data.csv"):
        print("Warning: Using synthetic data.")
        train_model("data/training_data.csv", epochs=args.epochs, batch_size=args.batch_size, lr=args.lr, model_name=args.name)
    else:
        print("Please run 'python src/ml/synthetic_data.py' first to generate dummy data.")