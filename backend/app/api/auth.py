from __future__ import annotations

import re
import sqlite3
from pathlib import Path
from typing import Any

from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field

from ..auth import get_current_token, get_current_user
from ..core.config import settings
from ..runtime import auth_service

router = APIRouter(prefix="/api/auth", tags=["auth"])

USERNAME_PATTERN = re.compile(r"^[A-Za-z0-9_.-]{3,32}$")


class AuthPayload(BaseModel):
    username: str = Field(min_length=3, max_length=32)
    password: str = Field(min_length=8, max_length=128)


def _validate_username(username: str) -> str:
    clean = username.strip()
    if not USERNAME_PATTERN.fullmatch(clean):
        raise HTTPException(
            status_code=400,
            detail=(
                "Username must be 3-32 chars and only contain letters, numbers, "
                "underscore, dot, or hyphen."
            ),
        )
    return clean


def _build_auth_response(user: dict[str, Any]) -> dict[str, Any]:
    token = auth_service.create_session(int(user["id"]))
    return {
        "access_token": token,
        "token_type": "bearer",
        "user": user,
    }


@router.post("/register")
async def register(payload: AuthPayload) -> dict[str, Any]:
    username = _validate_username(payload.username)
    try:
        user = auth_service.register(username, payload.password)
        return _build_auth_response(user)
    except sqlite3.IntegrityError:
        raise HTTPException(status_code=409, detail="Username already exists") from None
    except sqlite3.Error as exc:
        raise HTTPException(
            status_code=503,
            detail=f"Auth database unavailable: {exc}",
        ) from exc


@router.post("/login")
async def login(payload: AuthPayload) -> dict[str, Any]:
    try:
        user = auth_service.authenticate(payload.username, payload.password)
        if user is None:
            raise HTTPException(status_code=401, detail="Invalid username or password")
        return _build_auth_response(user)
    except sqlite3.Error as exc:
        raise HTTPException(
            status_code=503,
            detail=f"Auth database unavailable: {exc}",
        ) from exc


@router.post("/logout")
async def logout(token: str = Depends(get_current_token)) -> dict[str, str]:
    try:
        auth_service.revoke_session(token)
    except sqlite3.Error as exc:
        raise HTTPException(
            status_code=503,
            detail=f"Auth database unavailable: {exc}",
        ) from exc
    return {"status": "ok"}


@router.get("/me")
async def me(user: dict[str, Any] = Depends(get_current_user)) -> dict[str, Any]:
    return {"user": user}


@router.get("/user-guide")
async def user_guide() -> FileResponse:
    guide_path = settings.repo_root / "docs" / "LabFlow前端用户操作说明.md"
    if not Path(guide_path).exists():
        raise HTTPException(status_code=404, detail="User guide not found")
    return FileResponse(path=guide_path, media_type="text/markdown; charset=utf-8")
