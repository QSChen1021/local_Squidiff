# LabFlow Backend (MVP)

## Run locally

```bash
pip install -r requirements.txt -r backend/requirements.txt
uvicorn backend.app.main:app --reload --host 0.0.0.0 --port 8000
```

## Health check

`GET /api/health`

## Core endpoints

- `POST /api/datasets/upload`
- `POST /api/datasets/{id}/validate`
- `POST /api/jobs/train`
- `POST /api/jobs/predict`
- `GET /api/jobs/{id}`
- `GET /api/results/{id}`

## Notes

- State is file-based JSON under `backend/state/`.
- Set `LABFLOW_DRY_RUN=true` to run without heavy model execution.
- R conversion runtime options:
- `LABFLOW_R_EXEC_MODE=direct|cmd_conda`
- `LABFLOW_RSCRIPT_BIN=Rscript`
- `LABFLOW_R_CONDA_ENV=<your_conda_r_env_name>`
- `LABFLOW_R_CONDA_BAT=<path-to-conda.bat>`
