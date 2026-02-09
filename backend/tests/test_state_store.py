from __future__ import annotations

from backend.app.storage.state_manager import JsonStateStore


def test_create_and_update_dataset(tmp_path):
    store = JsonStateStore(tmp_path)
    ds = store.create_dataset({"name": "demo", "status": "uploaded"})
    assert ds["id"]
    assert ds["name"] == "demo"

    updated = store.update_dataset(ds["id"], {"status": "validated"})
    assert updated is not None
    assert updated["status"] == "validated"
