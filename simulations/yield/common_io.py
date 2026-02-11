from __future__ import annotations

import json
import logging
import os
import platform
import subprocess
import sys
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _json_default(obj: Any) -> Any:
    if is_dataclass(obj):
        return asdict(obj)
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, default=_json_default))


def stable_seed(base_seed: int, *indices: int) -> int:
    """Deterministic 32-bit seed from a base seed + integer indices.

    This avoids Python's randomized hash and stays stable across runs.
    """
    value = base_seed & 0xFFFFFFFF
    for idx in indices:
        value = (value * 1664525 + (idx & 0xFFFFFFFF) + 1013904223) & 0xFFFFFFFF
    return int(value)


def get_git_hash(cwd: Path) -> str | None:
    try:
        output = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=str(cwd),
            stderr=subprocess.STDOUT,
        )
    except (subprocess.SubprocessError, FileNotFoundError):
        return None
    return output.decode("utf-8").strip() or None


def _cpu_model_darwin() -> str | None:
    try:
        output = subprocess.check_output(
            ["sysctl", "-n", "machdep.cpu.brand_string"],
            stderr=subprocess.STDOUT,
        )
    except (subprocess.SubprocessError, FileNotFoundError):
        return None
    value = output.decode("utf-8").strip()
    return value or None


def _cpu_model_linux() -> str | None:
    cpuinfo = Path("/proc/cpuinfo")
    if not cpuinfo.exists():
        return None
    for line in cpuinfo.read_text().splitlines():
        if "model name" in line:
            _, value = line.split(":", 1)
            return value.strip() or None
    return None


def get_cpu_model() -> str | None:
    system = platform.system()
    if system == "Darwin":
        return _cpu_model_darwin()
    if system == "Linux":
        return _cpu_model_linux()
    value = platform.processor()
    return value or None


def build_run_manifest(cwd: Path, cli_args: Iterable[str]) -> dict:
    timestamp = datetime.now(timezone.utc).isoformat()
    manifest = {
        "timestamp_utc": timestamp,
        "os": platform.platform(),
        "python_version": sys.version.split()[0],
        "cpu_model": get_cpu_model(),
        "git_hash": get_git_hash(cwd),
        "cwd": str(cwd),
        "cli_args": list(cli_args),
    }
    return manifest


def write_run_manifest(output_dir: Path, cwd: Path, cli_args: Iterable[str]) -> Path:
    manifest = build_run_manifest(cwd, cli_args)
    path = output_dir / "run_manifest.json"
    write_json(path, manifest)
    return path


def setup_logging(log_dir: Path, stage: str) -> tuple[logging.Logger, Path]:
    ensure_dir(log_dir)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    log_path = log_dir / f"{stage}_{timestamp}.log"

    logger = logging.getLogger(stage)
    logger.setLevel(logging.INFO)
    logger.handlers = []
    logger.propagate = False

    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")

    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    return logger, log_path


def validate_required_keys(payload: dict, required_keys: Iterable[str]) -> None:
    missing = [key for key in required_keys if key not in payload]
    if missing:
        raise ValueError(f"Missing required keys: {missing}")
