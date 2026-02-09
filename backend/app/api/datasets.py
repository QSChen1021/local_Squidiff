from __future__ import annotations

import shutil
import re
from pathlib import Path

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from pydantic import BaseModel

from ..core.config import settings
from ..runtime import store
from ..services.data_validator import validate_h5ad
from ..services.seurat_converter import SUPPORTED_INPUTS, convert_to_h5ad

router = APIRouter(prefix="/api/datasets", tags=["datasets"])


class ValidatePayload(BaseModel):
    use_drug_structure: bool = False
    auto_convert: bool = True
    r_exec_mode: str | None = None
    r_conda_env: str | None = None
    r_conda_bat: str | None = None
    rscript_bin: str | None = None


def _save_upload(file: UploadFile, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / file.filename
    with destination.open("wb") as out:
        shutil.copyfileobj(file.file, out)
    return destination


def _safe_dataset_folder_name(name: str) -> str:
    cleaned = re.sub(r"[^a-zA-Z0-9._-]+", "_", name).strip("._")
    return cleaned or "dataset"


@router.get("")
def list_datasets() -> dict[str, object]:
    return {"items": store.list_datasets()}


@router.post("/upload")
def upload_dataset(
    file: UploadFile = File(...),
    smiles_file: UploadFile | None = File(default=None),
    dataset_name: str | None = Form(default=None),
) -> dict[str, object]:
    suffix = Path(file.filename).suffix.lower()
    if suffix not in SUPPORTED_INPUTS:
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported dataset format {suffix}. Supported: {SUPPORTED_INPUTS}",
        )

    folder_name = _safe_dataset_folder_name(dataset_name) if dataset_name else None
    upload_dir = (
        settings.upload_dir / folder_name if folder_name else settings.upload_dir
    )
    raw_path = _save_upload(file, upload_dir)
    smiles_path = None
    if smiles_file:
        smiles_path = str(_save_upload(smiles_file, upload_dir))

    record = store.create_dataset(
        {
            "name": dataset_name or raw_path.stem,
            "status": "uploaded",
            "input_type": suffix.replace(".", ""),
            "path_raw": str(raw_path),
            "path_h5ad": str(raw_path) if suffix == ".h5ad" else None,
            "smiles_path": smiles_path,
            "validation": None,
        }
    )
    return {"dataset": record}


@router.post("/{dataset_id}/validate")
def validate_dataset(dataset_id: str, payload: ValidatePayload) -> dict[str, object]:
    dataset = store.get_dataset(dataset_id)
    if dataset is None:
        raise HTTPException(status_code=404, detail="Dataset not found")

    raw_path = Path(dataset["path_raw"])
    h5ad_path = Path(dataset["path_h5ad"]) if dataset.get("path_h5ad") else None

    if payload.auto_convert and (h5ad_path is None or not h5ad_path.exists()):
        try:
            h5ad_path = convert_to_h5ad(
                input_path=raw_path,
                output_dir=settings.upload_dir / "converted",
                r_exec_mode=payload.r_exec_mode,
                r_conda_env=payload.r_conda_env,
                r_conda_bat=payload.r_conda_bat,
                rscript_bin=payload.rscript_bin,
            )
        except Exception as exc:  # noqa: BLE001
            store.update_dataset(
                dataset_id,
                {
                    "status": "invalid",
                    "validation": {
                        "valid": False,
                        "errors": [str(exc)],
                        "warnings": [],
                    },
                },
            )
            raise HTTPException(status_code=400, detail=str(exc)) from exc

    target = h5ad_path if h5ad_path else raw_path
    report = validate_h5ad(
        data_path=target,
        use_drug_structure=payload.use_drug_structure,
    )

    updated = store.update_dataset(
        dataset_id,
        {
            "path_h5ad": str(target),
            "status": "validated" if report["valid"] else "invalid",
            "validation": report,
        },
    )
    return {"dataset": updated, "validation": report}
