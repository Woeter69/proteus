import torch
import torch.nn as nn
from transformers import AutoModel

class PolymerPredictor(nn.Module):
    def __init__(self, model_name="seyonec/ChemBERTa-zinc-base-v1", hidden_dim=64):
        super(PolymerPredictor, self).__init__()
        
        # 1. The Pre-trained Chemical Brain (Transformer)
        try:
            self.bert = AutoModel.from_pretrained(model_name, use_safetensors=True)
        except Exception:
             # Fallback if safetensors unavailable
             print("[!] SafeTensors failed, falling back to standard load.")
             self.bert = AutoModel.from_pretrained(model_name)

        embedding_dim = self.bert.config.hidden_size # Usually 768

        # 2. The Prediction Head
        # Just use the CLS embedding
        self.regressor = nn.Sequential(
            nn.Linear(embedding_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim, 1) # Output: Predicted Rg
        )

    def forward(self, input_ids, attention_mask):
        # Pass through Transformer
        outputs = self.bert(input_ids=input_ids, attention_mask=attention_mask)
        
        # Get the "CLS" token embedding (representation of the whole molecule)
        cls_embedding = outputs.last_hidden_state[:, 0, :]
        
        # Predict Rg
        return self.regressor(cls_embedding)
