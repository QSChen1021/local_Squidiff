from __future__ import annotations

from pathlib import Path

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from ..runtime import job_queue, store

router = APIRouter(prefix="/api/jobs", tags=["jobs"])


class TrainJobPayload(BaseModel):
    dataset_id: str
    gene_size: int = Field(gt=0)
    output_dim: int = Field(gt=0)
    use_drug_structure: bool = False
    control_dataset_id: str | None = None
    batch_size: int = 64
    lr: float = 1e-4


class PredictJobPayload(BaseModel):
    dataset_id: str
    model_id: str | None = None
    model_path: str | None = None
    gene_size: int = Field(gt=0)
    output_dim: int = Field(gt=0)
    use_drug_structure: bool = False


@router.get("")
def list_jobs() -> dict[str, object]:
    return {"items": store.list_jobs()}


@router.get("/{job_id}")
def get_job(job_id: str) -> dict[str, object]:
    job = store.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return {"job": job}


@router.get("/{job_id}/log")
def get_job_log(job_id: str) -> dict[str, str]:
    job = store.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")

    log_path = job.get("log_path")
    if not isinstance(log_path, str):
        raise HTTPException(status_code=404, detail="Log not found")

    path = Path(log_path)
    if not path.exists():
        raise HTTPException(status_code=404, detail="Log file not found")

    text = path.read_text(encoding="utf-8", errors="replace")
    if len(text) > 20000:
        text = text[-20000:]
    return {"log": text}


@router.post("/train")
def submit_train_job(payload: TrainJobPayload) -> dict[str, object]:
    dataset = store.get_dataset(payload.dataset_id)
    if dataset is None:
        raise HTTPException(status_code=404, detail="Dataset not found")

    job = store.create_job(
        {
            "type": "train",
            "status": "queued",
            "dataset_id": payload.dataset_id,
            "params": payload.model_dump(),
            "error_msg": None,
            "started_at": None,
            "ended_at": None,
        }
    )
    job_queue.enqueue(job["id"])
    return {"job": job}


@router.post("/predict")
def submit_predict_job(payload: PredictJobPayload) -> dict[str, object]:
    dataset = store.get_dataset(payload.dataset_id)
    if dataset is None:
        raise HTTPException(status_code=404, detail="Dataset not found")

    if payload.model_id is None and payload.model_path is None:
        raise HTTPException(
            status_code=400, detail="model_id or model_path is required"
        )

    job = store.create_job(
        {
            "type": "predict",
            "status": "queued",
            "dataset_id": payload.dataset_id,
            "params": payload.model_dump(),
            "error_msg": None,
            "started_at": None,
            "ended_at": None,
        }
    )
    job_queue.enqueue(job["id"])
    return {"job": job}
