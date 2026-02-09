from __future__ import annotations

from .core.config import settings
from .services.job_queue import JobQueue
from .services.squidiff_runner import SquidiffRunner
from .storage.state_manager import JsonStateStore

store = JsonStateStore(settings.state_dir)
runner = SquidiffRunner(settings.repo_root)
job_queue = JobQueue(store=store, runner=runner)
