from __future__ import annotations

import sqlite3
from typing import Any

from fastapi import Header, HTTPException

from .core.config import settings
from .runtime import auth_service


def _extract_bearer_token(authorization: str | None) -> str:
    if not authorization:
        raise HTTPException(status_code=401, detail="Missing Authorization header")
    parts = authorization.split(" ", 1)
    if len(parts) != 2 or parts[0].lower() != "bearer" or not parts[1].strip():
        raise HTTPException(status_code=401, detail="Invalid Authorization header")
    return parts[1].strip()


def get_current_token(authorization: str | None = Header(None)) -> str:
    return _extract_bearer_token(authorization)


def get_current_user(authorization: str | None = Header(None)) -> dict[str, Any]:
    token = _extract_bearer_token(authorization)
    try:
        user = auth_service.get_user_by_token(token)
    except sqlite3.Error as exc:
        raise HTTPException(
            status_code=503, detail=f"Auth database unavailable: {exc}"
        ) from exc
    if user is None:
        raise HTTPException(status_code=401, detail="Invalid or expired token")
    return user


def require_auth(authorization: str | None = Header(None)) -> dict[str, Any] | None:
    if not settings.auth_required:
        return None
    token = _extract_bearer_token(authorization)
    try:
        user = auth_service.get_user_by_token(token)
    except sqlite3.Error as exc:
        raise HTTPException(
            status_code=503, detail=f"Auth database unavailable: {exc}"
        ) from exc
    if user is None:
        raise HTTPException(status_code=401, detail="Invalid or expired token")
    return user
