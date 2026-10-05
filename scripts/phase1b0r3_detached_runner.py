#!/usr/bin/env python3
"""Run the frozen Phase 1B-0R3 scorer while ignoring terminal interrupts."""
from __future__ import annotations

import runpy
import signal
import os
import subprocess
import sys
from pathlib import Path


PIPELINE = Path(__file__).with_name("phase1b0r3_score_pipeline.py").resolve()

# The scoring process is checkpointed and can run for many hours. Ignore the
# SIGINT sent when its launching terminal/session closes; failures still exit
# normally and can resume from the validated checkpoint.
signal.signal(signal.SIGINT, signal.SIG_IGN)
signal.signal(signal.SIGHUP, signal.SIG_IGN)
subprocess.Popen(
    ["/usr/bin/caffeinate", "-i", "-m", "-w", str(os.getpid())],
    stdin=subprocess.DEVNULL,
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
    start_new_session=True,
)
sys.argv = [str(PIPELINE), "--phase", "inference", "--stage", "full", "--cache-mode", "readwrite"]
runpy.run_path(str(PIPELINE), run_name="__main__")
