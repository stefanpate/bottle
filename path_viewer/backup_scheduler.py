import logging
import os
import subprocess
import threading
import time
from pathlib import Path

import streamlit as st

log = logging.getLogger(__name__)

BACKUP_SCRIPT = Path(os.environ.get("BACKUP_SCRIPT", "/app/scripts/backup_feedback.sh"))
FEEDBACK_ROOT = Path(os.environ.get("FEEDBACK_ROOT", "/feedback"))
LAST_BACKUP_MARKER = FEEDBACK_ROOT / ".last_backup"
INTERVAL_HOURS = float(os.environ.get("BACKUP_INTERVAL_HOURS", "24"))
CHECK_SECONDS = float(os.environ.get("BACKUP_CHECK_SECONDS", "3600"))


def _due() -> bool:
    if not LAST_BACKUP_MARKER.exists():
        return True
    age_h = (time.time() - LAST_BACKUP_MARKER.stat().st_mtime) / 3600.0
    return age_h >= INTERVAL_HOURS


def _run_backup() -> None:
    if not BACKUP_SCRIPT.exists():
        log.warning("backup script not found at %s", BACKUP_SCRIPT)
        return
    if not os.environ.get("BACKUP_REMOTE"):
        log.info("BACKUP_REMOTE not set; skipping backup")
        return
    try:
        result = subprocess.run(
            ["sh", str(BACKUP_SCRIPT)],
            capture_output=True,
            text=True,
            timeout=300,
        )
        if result.returncode == 0:
            LAST_BACKUP_MARKER.parent.mkdir(parents=True, exist_ok=True)
            LAST_BACKUP_MARKER.touch()
            log.info("backup ok: %s", result.stdout.strip())
        else:
            log.error("backup failed (rc=%s): %s", result.returncode, result.stderr.strip())
    except Exception:
        log.exception("backup script raised")


def _loop() -> None:
    while True:
        if _due():
            _run_backup()
        time.sleep(CHECK_SECONDS)


@st.cache_resource
def start_backup_thread_once() -> threading.Thread:
    t = threading.Thread(target=_loop, name="feedback-backup", daemon=True)
    t.start()
    log.info("feedback backup thread started")
    return t
