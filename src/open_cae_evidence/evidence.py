"""Common evidence persistence and integrity checks, independent of the solver."""

import hashlib
import json
from pathlib import PurePosixPath
import platform
import subprocess
import sys
import os


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest(directory):
    return {
        str(p.relative_to(directory)): sha256(p)
        for p in sorted(directory.rglob("*"))
        if p.is_file() and p != directory / "manifest.json"
    }


def verify_manifest(directory):
    expected = json.loads((directory / "manifest.json").read_text())
    for name in expected:
        path = PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts or str(path) != name:
            raise ValueError("Unsafe evidence manifest path")
    if any(p.is_symlink() for p in directory.rglob("*")):
        raise ValueError("Evidence symlinks are unsupported")
    actual = manifest(directory)
    if set(actual) != set(expected):
        raise ValueError("Evidence file set differs from manifest")
    for name, digest in actual.items():
        if expected[name] != digest:
            raise ValueError(f"Evidence hash mismatch: {name}")


def provenance(executable, adapter):
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True
    )
    dirty = subprocess.run(
        ["git", "status", "--porcelain"], capture_output=True, text=True
    )
    from pathlib import Path

    return {
        "python": sys.version,
        "platform": platform.platform(),
        "solver_path": executable,
        "solver_binary_sha256": sha256(Path(executable)),
        "solver_name": adapter.name,
        "adapter_version": adapter.version,
        "git_revision": revision.stdout.strip() or None,
        "working_tree_dirty": bool(dirty.stdout.strip()),
        "github_run_id": os.environ.get("GITHUB_RUN_ID"),
        "github_sha": os.environ.get("GITHUB_SHA"),
        "threads": 1,
    }
