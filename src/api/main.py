import os
import yaml
import json
from fastapi import FastAPI, WebSocket, WebSocketDisconnect, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional, List
from dotenv import load_dotenv
from src.api.orchestrator import PipelineOrchestrator

load_dotenv()

app = FastAPI(title="Structural Prioritization API")

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

class ChatMessage(BaseModel):
    role: str
    content: str

class ChatRequest(BaseModel):
    message: str
    history: Optional[List[ChatMessage]] = []
    context: Optional[dict] = {}


def build_system_prompt(context: dict) -> str:
    gene = context.get("gene_symbol", "unknown gene").upper()
    structure = context.get("structure_id", "unknown structure").upper()
    coverage = context.get("mapping_coverage", "N/A")
    top_variants = context.get("top_variants", [])
    enrichment = context.get("enrichment", [])

    variants_text = ""
    for v in top_variants[:5]:
        rsid = v.get("rsid", "?")
        aa = v.get("amino_acid_change", "?")
        domain = v.get("domain_region", "?")
        score = v.get("priority_score", "?")
        cadd = v.get("cadd_phred", "?")
        variants_text += f"  - {rsid} ({aa}): domain={domain}, priority_score={score}, CADD={cadd}\n"

    enrichment_text = ""
    for e in enrichment:
        feat = e.get("feature", "?")
        status = e.get("status", "?")
        if status == "success":
            OR = e.get("odds_ratio", "?")
            p = e.get("p_value", "?")
            enrichment_text += f"  - {feat}: OR={OR:.2f}, p={p:.2e}\n"
        else:
            reason = e.get("reason", "")
            enrichment_text += f"  - {feat}: {status} ({reason})\n"

    return f"""You are a structural genomics AI assistant embedded in RAREMISS, a research pipeline for prioritizing rare missense variants in ion channel and disease-associated proteins.

Current analysis context:
- Gene: {gene}
- Structure: PDB {structure}
- Structural mapping coverage: {coverage}%
- Pipeline: Ensembl rare missense variants → coordinate correction → domain annotation → spatial mapping → priority scoring → Fisher enrichment

Top prioritized variants:
{variants_text if variants_text else "  Not yet available."}

Structural enrichment results:
{enrichment_text if enrichment_text else "  Not yet computed or data unavailable."}

Your role:
- Answer questions about this specific analysis, the variants, the biological significance of the enrichment results, and what the scores mean
- Explain concepts like CADD scores, odds ratios, Fisher's exact test, structural domains, binding sites, and interface residues in plain language
- Suggest follow-up experiments or literature searches when relevant
- Be concise and scientific but accessible to researchers who may not be bioinformaticians
- If you don't have enough data to answer confidently, say so clearly
- Do NOT make up variant data, clinical claims, or drug interactions not supported by the context above

Keep responses focused and to the point. Use markdown formatting where helpful."""


@app.post("/chat")
async def chat(req: ChatRequest):
    api_key = os.getenv("GEMINI_API_KEY")
    if not api_key or api_key == "your_api_key_here":
        raise HTTPException(
            status_code=503,
            detail="Gemini API key not configured. Add GEMINI_API_KEY to your .env file."
        )

    try:
        from google import genai
        from google.genai import types

        client = genai.Client(api_key=api_key)

        system_prompt = build_system_prompt(req.context)

        history_contents = []
        for msg in (req.history or []):
            role = "user" if msg.role == "user" else "model"
            history_contents.append(
                types.Content(role=role, parts=[types.Part(text=msg.content)])
            )

        response = client.models.generate_content(
            model="gemini-2.5-flash",
            contents=history_contents + [
                types.Content(role="user", parts=[types.Part(text=req.message)])
            ],
            config=types.GenerateContentConfig(
                system_instruction=system_prompt,
                temperature=0.7,
                max_output_tokens=1024,
            ),
        )

        return {"response": response.text}

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/config")
def get_config():
    with open(CONFIG_PATH, 'r') as f:
        return yaml.safe_load(f)

@app.post("/config")
def update_config(config: ConfigUpdate):
    with open(CONFIG_PATH, 'r') as f:
        data = yaml.safe_load(f)
    data.update(config.model_dump())
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
                pipeline_succeeded = True
                for i in range(len(orchestrator.steps)):
                    await websocket.send_json({"type": "step_start", "index": i})
                    success, desc = await orchestrator.run_step(i)
                    while not orchestrator.output_queue.empty():
                        line = await orchestrator.output_queue.get()
                        await websocket.send_json({"type": "log", "message": line})
                    await websocket.send_json({"type": "step_end", "index": i, "success": success})
                    if not success:
                        pipeline_succeeded = False
                        await websocket.send_json({
                            "type": "pipeline_failed",
                            "message": f"Step '{desc}' failed. Check logs."
                        })
                        break
                if pipeline_succeeded:
                    await websocket.send_json({"type": "pipeline_complete"})
    except WebSocketDisconnect:
        pass
    except Exception as e:
        await websocket.send_json({"type": "error", "message": str(e)})

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
