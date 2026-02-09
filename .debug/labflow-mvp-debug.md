# LabFlow MVP Debug Record

## Metadata
- Module name: labflow-mvp
- Created at: 2026-02-09
- Last updated: 2026-02-10
- Related files:
- `backend/app/main.py`
- `backend/app/api/datasets.py`
- `backend/app/api/jobs.py`
- `backend/app/api/results.py`
- `backend/app/storage/state_manager.py`
- `backend/app/services/job_queue.py`
- `backend/app/services/squidiff_runner.py`
- `backend/app/services/seurat_converter.py`
- `frontend/src/App.tsx`
- `frontend/src/services/api.ts`
- `infra/docker-compose.yml`
- Dependency modules:
- `train_squidiff.py`
- `sample_squidiff.py`

## Runtime Context and Test Rules
- Runtime environment: Local Windows (PowerShell workspace `E:\Development\Squidiff`)
- SSH mode (if remote): Not used in this implementation round
- Remote project path (if remote): N/A
- Validation/Checkfix execution mode: Run commands directly in local shell
- R execution constraint (confirmed by user): R conda env must be activated via `cmd` (not PowerShell) before running `Rscript`.
- R config strategy: support both `.env` defaults (`LABFLOW_R_*`) and per-request frontend overrides (`r_exec_mode`, `r_conda_env`, `r_conda_bat`, `rscript_bin`).

## Context Network
- File layout
- New MVP modules added under `backend/`, `frontend/`, and `infra/`.
- Existing training/predict scripts are kept unchanged and invoked via service wrapper.
- Function call chain
- API (`datasets/jobs/results`) -> runtime singleton (`store`, `job_queue`) -> service layer (`validator/converter/runner`) -> JSON state + scripts.
- Variable/data dependencies
- Dataset upload writes raw paths into `datasets.json`.
- Job payload references dataset IDs and optional model IDs.
- Worker resolves paths then writes model/result entries into `models.json`/`results.json`.
- Data flow
- User upload -> validate/convert -> queue job -> run training/predict -> generate result assets -> query job/result APIs.
- Frontend full flow
- Upload file -> validate (with R runtime config) -> submit train job -> poll job status -> display train outputs/result assets/logs.

## Debug History
### [2026-02-09 23:xx] Bootstrap no-SQL MVP
- Problem
- Need to start implementation from PRD with no SQL requirement and minimal intranet scope.
- Root cause
- Repository only had research scripts, no service architecture.
- Solution
- Added minimal backend API, JSON state store, background worker queue, script wrappers, and frontend shell.
- Code changes (files/functions)
- Added backend runtime and API endpoints, plus file-based persistence and Docker compose deployment.
- Verification results
- `ruff check backend`: passed.
- `ruff format --check backend`: passed.
- `python -m compileall backend/app`: passed.
- `python -c "from fastapi.testclient import TestClient; ... /api/health"`: passed (`200 {"status":"ok"}`).
- `uv run pytest -q backend/tests`: blocked by dependency resolution conflict in current environment (`scanpy` vs broad `requires-python >=3.8` resolution path). Not a code syntax failure.
- `python -m pytest -q backend/tests`: blocked because local Python environment does not include `pytest`.
- Impact assessment
- Existing model scripts untouched; new code isolated under `backend/`, `frontend/`, `infra/`.

### [2026-02-09 23:xx] Startup dependency decoupling
- Problem
- Backend import path required `scanpy` at process startup, causing health check startup failure on lean environments.
- Root cause
- `scanpy` was imported at module import time in validator/runner service modules.
- Solution
- Converted `scanpy` imports to lazy runtime imports inside function scope.
- Code changes (files/functions)
- `backend/app/services/data_validator.py` (`validate_h5ad`)
- `backend/app/services/squidiff_runner.py` (`run_predict`)
- Verification results
- `python -c "from backend.app.main import app; print(app.title)"` prints app title successfully.
- `ruff`/`compileall` still pass.
- Impact assessment
- Service can boot earlier; validation/predict endpoints now return explicit dependency guidance if `scanpy` missing.

### [2026-02-10 00:xx] Frontend full workflow + Windows cmd_conda R support
- Problem
- Need to complete end-to-end UI flow and support user-specified Conda R environment on Windows CMD for Seurat conversion.
- Root cause
- Initial frontend was only a shell; backend R conversion only supported direct `Rscript`.
- Solution
- Implemented complete frontend workflow (upload/validate/train/poll/result).
- Added backend support for `cmd_conda` execution mode and per-request R runtime overrides.
- Added result asset URL mapping and job log API for result page display.
- Code changes (files/functions)
- `backend/app/core/config.py` (new `LABFLOW_R_*` settings)
- `backend/app/services/seurat_converter.py` (`_build_r_command`, `convert_to_h5ad`)
- `backend/app/api/datasets.py` (`ValidatePayload` extended with R runtime fields)
- `backend/app/api/jobs.py` (`GET /api/jobs/{job_id}/log`)
- `backend/app/api/results.py` (model detail API + asset URL + asset file serving)
- `frontend/src/services/api.ts` (full API client for flow)
- `frontend/src/App.tsx` (step-by-step workflow implementation)
- `frontend/src/styles/tokens.css` (form/status/result styles)
- `infra/.env.example` and `infra/docker-compose.yml` (runtime config exposure)
- Verification results
- Backend:
- `ruff check backend`: passed.
- `ruff format --check backend`: passed.
- `python -m compileall backend/app`: passed.
- Frontend:
- `npm install`: passed.
- `npm run lint`: passed.
- `npm run build`: passed.
- Impact assessment
- Lab users can now complete the requested workflow from a single frontend page.
- Windows R/Conda activation requirement is now configurable and executable via CMD.

### [2026-02-10 00:xx] V2 PRD drafted: Seurat interactive selection + 500x500 pipeline
- Problem
- Users still struggle with manual Seurat filtering and metadata preparation before training.
- Root cause
- MVP assumes preprocessed h5ad-style input and lacks in-UI cluster/group selection + bounded preprocessing.
- Solution
- Added a new PRD describing:
- 1) WebUI Seurat inspection and UMAP interaction,
- 2) user-selected metadata mapping (`group_column`, `cluster_column`),
- 3) cell stratified downsampling to max 500,
- 4) DEG-based gene selection to max 500,
- 5) final 500x500 training matrix contract.
- Code changes (files/functions)
- `docs/PRD_Seurat交互筛选与500x500训练管线.md` (new)
- Verification results
- Documentation update only (no runtime behavior changed in this step).
- Impact assessment
- V2 scope is now explicit and implementation-ready for phased development.

## Open Issues
- Real-world Seurat conversion relies on local R/SeuratDisk availability.
- Production auth is intentionally simplified for MVP.
- V2 features in new PRD are not yet implemented in code.

## Technical Debt
- JSON file storage has limited concurrency compared with database-backed approach.
- Current UI is single-page workflow; multi-page routing and better UX states can be added later.
