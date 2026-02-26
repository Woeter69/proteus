import torch
from torch.utils.data import Dataset
from transformers import AutoTokenizer
import pandas as pd
import numpy as np

class PolymerDataset(Dataset):
    def __init__(self, csv_file, tokenizer_name="seyonec/ChemBERTa-zinc-base-v1", max_length=128):
        """
        Args:
            csv_file (str): Path to Kaggle CSV.
            tokenizer_name (str): HuggingFace tokenizer model.
            max_length (int): Max token length.
        """
        raw_data = pd.read_csv(csv_file)
        
        # Filter: Only rows with valid Rg values
        if 'Rg' in raw_data.columns:
            self.data = raw_data.dropna(subset=['Rg']).reset_index(drop=True)
            self.target_col = 'Rg'
        else:
            # Fallback for synthetic data
            self.data = raw_data
            self.target_col = 'rg'

        self.tokenizer = AutoTokenizer.from_pretrained(tokenizer_name)
        self.max_length = max_length

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        row = self.data.iloc[idx]
        
        # Handle Column Names (Kaggle 'SMILES' vs Synthetic 'smiles')
        smiles = row.get('SMILES') or row.get('smiles')
        
        # Preprocessing: Replace Wildcards '*' with Carbon 'C' for valid tokenization
        # (Polymer end-capping)
        if isinstance(smiles, str):
            smiles = smiles.replace("*", "C")
        
        target = row[self.target_col]

        encoding = self.tokenizer(
            smiles,
            add_special_tokens=True,
            max_length=self.max_length,
            return_token_type_ids=False,
            padding='max_length',
            truncation=True,
            return_attention_mask=True,
            return_tensors='pt',
        )

        return {
            'input_ids': encoding['input_ids'].flatten(),
            'attention_mask': encoding['attention_mask'].flatten(),
            'labels': torch.tensor(target, dtype=torch.float)
        }
