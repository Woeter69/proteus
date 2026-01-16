from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
import sys
from pathlib import Path

# Add the parent directory to sys.path to access 'src'
sys.path.append(str(Path(__file__).parent.parent))

from src import topology  # Example import to ensure connection works

app = FastAPI(title="Proteus API")

# Allow CORS for the frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],  # Next.js default port
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
def read_root():
    return {"message": "Proteus API is running"}

@app.get("/health")
def health_check():
    return {"status": "ok"}
