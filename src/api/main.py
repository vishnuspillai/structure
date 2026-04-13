import os
import yaml
import json
from fastapi import FastAPI, WebSocket, WebSocketDisconnect, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from typing import List, Optional
import asyncio
from src.api.orchestrator import PipelineOrchestrator

app = FastAPI(title="Structural Prioritization API")

# Setup CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CONFIG_PATH = os.path.join(ROOT_DIR, "config", "parameters.yaml")
DATA_DIR = os.path.join(ROOT_DIR, "data", "processed")

orchestrator = PipelineOrchestrator(ROOT_DIR)

class ConfigUpdate(BaseModel):
    af_threshold: float
    gene_symbol: str
    consequence_filter: str
    structure_id: Optional[str] = "7kox"
    species: Optional[str] = "homo_sapiens"

class PipelineRequest(BaseModel):
    gene_symbol: str
    af_threshold: float
    structure_id: str
    species: str


@app.get("/config")
def get_config():
    with open(CONFIG_PATH, 'r') as f:
        return yaml.safe_load(f)

@app.post("/config")
def update_config(config: ConfigUpdate):
    with open(CONFIG_PATH, 'r') as f:
        data = yaml.safe_load(f)
    
    data.update(config.dict())
    
    with open(CONFIG_PATH, 'w') as f:
        yaml.safe_dump(data, f)
    return {"status": "success"}

@app.post("/run_pipeline")
def trigger_pipeline(req: PipelineRequest):
    with open(CONFIG_PATH, 'r') as f:
        data = yaml.safe_load(f)
    
    data.update({
        "gene_symbol": req.gene_symbol,
        "af_threshold": req.af_threshold,
        "structure_id": req.structure_id,
        "species": req.species
    })
    
    with open(CONFIG_PATH, 'w') as f:
        yaml.safe_dump(data, f)
    
    # We return success here; the frontend will use the websocket to actually stream the run
    return {"status": "ready"}

@app.get("/steps")
def get_steps():
    return orchestrator.get_steps()

@app.get("/data/{filename}")
def get_data(filename: str):
    file_path = os.path.join(DATA_DIR, filename)
    if not os.path.exists(file_path):
        for f in os.listdir(DATA_DIR):
            if f.lower() == filename.lower():
                file_path = os.path.join(DATA_DIR, f)
                break
        else:
            raise HTTPException(status_code=404, detail="File not found")
    
    if filename.endswith(".csv"):
        import pandas as pd
        import numpy as np
        df = pd.read_csv(file_path)
        # Handle NaN/Inf for JSON compliance
        df = df.replace({np.nan: None, np.inf: None, -np.inf: None})
        return df.to_dict(orient="records")
    elif filename.endswith(".json"):
        with open(file_path, 'r') as f:
            return json.load(f)
    return {"error": "Unsupported file type"}

@app.websocket("/ws/run")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()
    try:
        while True:
            data = await websocket.receive_text()
            cmd = json.loads(data)
            
            if cmd.get("action") == "run_all":
                for i in range(len(orchestrator.steps)):
                    await websocket.send_json({"type": "step_start", "index": i})
                    success, desc = await orchestrator.run_step(i)
                    
                    while not orchestrator.output_queue.empty():
                        line = await orchestrator.output_queue.get()
                        await websocket.send_json({"type": "log", "message": line})
                        
                    await websocket.send_json({
                        "type": "step_end", 
                        "index": i, 
                        "success": success
                    })
                    if not success:
                        break
                await websocket.send_json({"type": "pipeline_complete"})
                
    except WebSocketDisconnect:
        pass
    except Exception as e:
        await websocket.send_json({"type": "error", "message": str(e)})

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
