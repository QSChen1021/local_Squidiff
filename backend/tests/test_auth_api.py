from __future__ import annotations

from uuid import uuid4

from fastapi.testclient import TestClient

from backend.app.main import app


def test_register_login_me_logout_flow() -> None:
    client = TestClient(app)
    username = f"user_{uuid4().hex[:8]}"
    password = "labflow-pass-123"

    register_resp = client.post(
        "/api/auth/register",
        json={"username": username, "password": password},
    )
    assert register_resp.status_code == 200
    register_payload = register_resp.json()
    token = register_payload["access_token"]
    assert register_payload["user"]["username"] == username

    me_resp = client.get("/api/auth/me", headers={"Authorization": f"Bearer {token}"})
    assert me_resp.status_code == 200
    assert me_resp.json()["user"]["username"] == username

    logout_resp = client.post(
        "/api/auth/logout", headers={"Authorization": f"Bearer {token}"}
    )
    assert logout_resp.status_code == 200

    me_after_logout = client.get(
        "/api/auth/me", headers={"Authorization": f"Bearer {token}"}
    )
    assert me_after_logout.status_code == 401


def test_login_with_wrong_password_returns_401() -> None:
    client = TestClient(app)
    username = f"user_{uuid4().hex[:8]}"
    password = "labflow-pass-123"

    client.post(
        "/api/auth/register",
        json={"username": username, "password": password},
    )
    bad_login = client.post(
        "/api/auth/login",
        json={"username": username, "password": "wrong-password"},
    )
    assert bad_login.status_code == 401
