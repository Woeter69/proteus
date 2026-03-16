#!/bin/bash
# Recreation script for the virtual environment using uv
echo "Recreating virtual environment with uv..."
rm -rf venv
uv venv venv
uv pip install -r requirements.txt
echo "Environment recreated in venv"
