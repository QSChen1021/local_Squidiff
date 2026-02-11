from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


def _as_bool(value: str, default: bool = False) -> bool:
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "y", "on"}


@dataclass(frozen=True)
class Settings:
    repo_root: Path
    backend_root: Path
    state_dir: Path
    auth_db_path: Path
    upload_dir: Path
    artifact_dir: Path
    max_upload_mb: int
    dry_run: bool
    auth_required: bool
    auth_session_ttl_hours: int
    rscript_bin: str
    r_exec_mode: str
    r_conda_env: str | None
    r_conda_bat: str | None


def load_settings() -> Settings:
    backend_root = Path(__file__).resolve().parents[2]
    repo_root = Path(__file__).resolve().parents[3]

    state_dir = Path(os.getenv("LABFLOW_STATE_DIR", backend_root / "state"))
    auth_db_path = Path(os.getenv("LABFLOW_AUTH_DB_PATH", state_dir / "auth.db"))
    upload_dir = Path(os.getenv("LABFLOW_UPLOAD_DIR", backend_root / "uploads"))
    artifact_dir = Path(os.getenv("LABFLOW_ARTIFACT_DIR", backend_root / "artifacts"))
    max_upload_mb = int(os.getenv("LABFLOW_MAX_UPLOAD_MB", "5120"))
    dry_run = _as_bool(os.getenv("LABFLOW_DRY_RUN"), default=False)
    auth_required = _as_bool(os.getenv("LABFLOW_AUTH_REQUIRED"), default=True)
    auth_session_ttl_hours = int(os.getenv("LABFLOW_AUTH_SESSION_TTL_HOURS", "12"))
    rscript_bin = os.getenv("LABFLOW_RSCRIPT_BIN", "Rscript")
    r_exec_mode = os.getenv("LABFLOW_R_EXEC_MODE", "direct")
    r_conda_env = os.getenv("LABFLOW_R_CONDA_ENV")
    r_conda_bat = os.getenv("LABFLOW_R_CONDA_BAT")

    state_dir.mkdir(parents=True, exist_ok=True)
    upload_dir.mkdir(parents=True, exist_ok=True)
    artifact_dir.mkdir(parents=True, exist_ok=True)

    return Settings(
        repo_root=repo_root,
        backend_root=backend_root,
        state_dir=state_dir,
        auth_db_path=auth_db_path,
        upload_dir=upload_dir,
        artifact_dir=artifact_dir,
        max_upload_mb=max_upload_mb,
        dry_run=dry_run,
        auth_required=auth_required,
        auth_session_ttl_hours=auth_session_ttl_hours,
        rscript_bin=rscript_bin,
        r_exec_mode=r_exec_mode,
        r_conda_env=r_conda_env,
        r_conda_bat=r_conda_bat,
    )


settings = load_settings()
