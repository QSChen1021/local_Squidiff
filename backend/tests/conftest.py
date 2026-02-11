from __future__ import annotations

import os

# Keep existing API tests unchanged; auth can be enabled explicitly in runtime env.
os.environ.setdefault("LABFLOW_AUTH_REQUIRED", "false")
