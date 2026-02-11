"""Runtime discovery API: conda envs and conda.bat candidates for R execution."""

from __future__ import annotations

import csv
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path

from fastapi import APIRouter, Query
from fastapi.responses import JSONResponse

router = APIRouter(prefix="/api/runtime", tags=["runtime"])


def _conda_bat_candidates_windows() -> list[str]:
    """Return list of conda.bat paths that exist (Windows)."""
    candidates: list[str] = []
    # PATH
    for base in os.environ.get("PATH", "").split(os.pathsep):
        p = Path(base) / "conda.bat"
        if p.exists():
            s = str(p.resolve())
            if s not in candidates:
                candidates.append(s)
    # Common install locations
    for base in [
        Path(os.environ.get("PROGRAMDATA", "C:\\ProgramData")) / "Miniconda3",
        Path(os.environ.get("PROGRAMDATA", "C:\\ProgramData")) / "miniconda3",
        Path(os.environ.get("LOCALAPPDATA", "")) / "Programs" / "Miniconda3",
        Path(os.environ.get("LOCALAPPDATA", "")) / "Programs" / "miniconda3",
        Path("F:/software/Miniconda3"),
        Path("C:/ProgramData/Miniconda3"),
    ]:
        if not base:
            continue
        bat = base / "condabin" / "conda.bat"
        if bat.exists():
            s = str(bat.resolve())
            if s not in candidates:
                candidates.append(s)
    return candidates


def _run_conda_env_list(conda_bat: str) -> list[str]:
    """Run conda env list for given conda.bat (Windows cmd), return env names."""
    env_names: list[str] = []
    try:
        if sys.platform == "win32":
            out = subprocess.run(
                ["cmd", "/c", f'"{conda_bat}"', "env", "list", "--json"],
                capture_output=True,
                text=True,
                timeout=15,
            )
        else:
            out = subprocess.run(
                ["conda", "env", "list", "--json"],
                capture_output=True,
                text=True,
                timeout=15,
            )
        if out.returncode != 0 or not out.stdout:
            # Fallback: parse text output
            if sys.platform == "win32":
                out2 = subprocess.run(
                    ["cmd", "/c", f'"{conda_bat}"', "env", "list"],
                    capture_output=True,
                    text=True,
                    timeout=15,
                )
            else:
                out2 = subprocess.run(
                    ["conda", "env", "list"],
                    capture_output=True,
                    text=True,
                    timeout=15,
                )
            text = (out2.stdout or "").strip()
            for line in text.splitlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                # First column is env name (may have * for current)
                parts = re.split(r"\s+", line, 1)
                name = parts[0].rstrip("*").strip()
                if name:
                    env_names.append(name)
            return env_names
        data = json.loads(out.stdout)
        # conda env list --json: envs can be list of paths or list of {name, prefix}
        for env in data.get("envs", []):
            if isinstance(env, dict) and "name" in env:
                env_names.append(str(env["name"]))
            elif isinstance(env, str):
                env_names.append(Path(env).name)
    except (json.JSONDecodeError, subprocess.TimeoutExpired, FileNotFoundError):
        pass
    return env_names


def _parse_int(value: str) -> int | None:
    try:
        return int(value.strip())
    except (ValueError, AttributeError):
        return None


def _parse_float(value: str) -> float | None:
    try:
        return float(value.strip())
    except (ValueError, AttributeError):
        return None


def _read_gpu_stats() -> dict[str, object]:
    query_fields = [
        "index",
        "name",
        "utilization.gpu",
        "utilization.memory",
        "memory.used",
        "memory.total",
        "temperature.gpu",
        "power.draw",
        "power.limit",
    ]
    cmd = [
        "nvidia-smi",
        f"--query-gpu={','.join(query_fields)}",
        "--format=csv,noheader,nounits",
    ]

    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=4)
    except FileNotFoundError:
        return {
            "available": False,
            "reason": "nvidia-smi is not installed or not in PATH",
            "gpus": [],
            "updated_at": datetime.now(timezone.utc).isoformat(),
        }
    except subprocess.TimeoutExpired:
        return {
            "available": False,
            "reason": "nvidia-smi timed out",
            "gpus": [],
            "updated_at": datetime.now(timezone.utc).isoformat(),
        }

    if out.returncode != 0:
        detail = (out.stderr or out.stdout or "").strip()
        return {
            "available": False,
            "reason": detail or f"nvidia-smi exited with code {out.returncode}",
            "gpus": [],
            "updated_at": datetime.now(timezone.utc).isoformat(),
        }

    reader = csv.reader(StringIO(out.stdout.strip()))
    gpu_rows: list[dict[str, object]] = []
    for row in reader:
        if len(row) < len(query_fields):
            continue
        gpu_rows.append(
            {
                "index": _parse_int(row[0]),
                "name": row[1].strip(),
                "utilization_gpu": _parse_float(row[2]),
                "utilization_memory": _parse_float(row[3]),
                "memory_used_mb": _parse_float(row[4]),
                "memory_total_mb": _parse_float(row[5]),
                "temperature_c": _parse_float(row[6]),
                "power_draw_w": _parse_float(row[7]),
                "power_limit_w": _parse_float(row[8]),
            }
        )

    return {
        "available": True,
        "reason": None,
        "gpus": gpu_rows,
        "updated_at": datetime.now(timezone.utc).isoformat(),
    }


@router.get("/conda-envs")
async def get_conda_envs(
    conda_bat: str | None = Query(
        None, description="Optional conda.bat path to list envs for"
    ),
) -> JSONResponse:
    """Return conda.bat candidates and conda env names for dropdowns.
    If conda_bat is provided, return envs for that path; else envs for first candidate.
    """
    if sys.platform == "win32":
        candidates = _conda_bat_candidates_windows()
    else:
        candidates = []
    conda_envs: list[str] = []
    if conda_bat:
        conda_envs = _run_conda_env_list(conda_bat)
    elif candidates:
        conda_envs = _run_conda_env_list(candidates[0])
    return JSONResponse(
        content={
            "conda_bat_candidates": candidates,
            "conda_envs": conda_envs,
        }
    )


@router.get("/gpu-stats")
async def get_gpu_stats() -> JSONResponse:
    """Return GPU utilization summary from nvidia-smi for UI polling."""
    return JSONResponse(content=_read_gpu_stats())
