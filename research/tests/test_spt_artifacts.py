"""Verify shared SPT publication and execution behavior without research datasets."""

import os
from pathlib import Path

import pytest

from research.spt.artifacts import (
    atomic_json,
    atomic_output,
    atomic_text,
    fingerprint,
    publish,
)
from research.spt.execution import execute_jobs


def _square(value: int) -> tuple[int, int]:
    if value == 0:
        raise ValueError("Invalid job")
    return value * value, os.getpid()


@pytest.mark.parametrize("workers", [1, 2])
def test_jobs_isolate_failures_and_complete_remaining_work(workers: int) -> None:
    """Return every job once and convert failures in the coordinator in both modes."""
    failures: list[tuple[int, str]] = []

    def failed(job: int, error: Exception) -> tuple[int, int]:
        failures.append((job, str(error)))
        return -1, os.getpid()

    results = dict(execute_jobs([2, 0, 3], _square, workers=workers, failed=failed))
    assert failures == [(0, "Invalid job")]
    assert {job: result[0] for job, result in results.items()} == {2: 4, 0: -1, 3: 9}
    assert all((results[job][1] == os.getpid()) == (workers == 1) for job in (2, 3))
    assert results[0][1] == os.getpid()


@pytest.mark.parametrize("jobs", [[], [3]])
def test_small_batches_never_start_workers(
    jobs: list[int], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Avoid worker startup for unchanged selections and single pending histories."""

    def unexpected(*args: object, **kwargs: object) -> None:
        pytest.fail("Small batches must not start worker processes")

    monkeypatch.setattr("research.spt.execution.ProcessPoolExecutor", unexpected)
    results = list(
        execute_jobs(jobs, _square, workers=8, failed=lambda _job, _error: (-1, 0))
    )
    assert results == [(job, (job * job, os.getpid())) for job in jobs]
    with pytest.raises(ValueError, match="workers must be positive"):
        list(execute_jobs(jobs, _square, workers=0, failed=lambda *_: (-1, 0)))


@pytest.mark.parametrize("failure", ["changed source", "manifest replacement"])
def test_failed_publication_retains_previous_generation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    """Keep the published pointer and files when a replacement cannot complete."""
    source = tmp_path / "source.txt"
    source.write_text("input")
    directory = tmp_path / "output"
    staging = directory / "staging"
    staging.mkdir(parents=True)
    (staging / "data.txt").write_text("previous result")
    records = [fingerprint(source)]
    previous = publish(
        directory,
        staging,
        {},
        (source,),
        records,
        catalog=lambda generation: {"data": str(generation / "data.txt")},
    )
    pointer = directory / "manifest.json"
    before = pointer.read_bytes()
    staging.mkdir()
    (staging / "data.txt").write_text("new result")
    if failure == "changed source":
        source.write_text("modified input")
    else:
        replace = Path.replace

        def fail_pointer(path: Path, target: Path) -> Path:
            if target == pointer:
                raise OSError("Cannot publish pointer")
            return replace(path, target)

        monkeypatch.setattr(Path, "replace", fail_pointer)
    with pytest.raises((RuntimeError, OSError)):
        publish(directory, staging, {}, (source,), records)
    assert pointer.read_bytes() == before
    assert (previous / "data.txt").read_text() == "previous result"
    assert list(directory.glob("generation-*")) == [previous]
    assert not list(directory.glob(".partial-*"))


def test_atomic_writes_retain_previous_content_on_failure(tmp_path: Path) -> None:
    """Never replace a report with partial text or invalid JSON."""
    path = tmp_path / "report.json"
    atomic_json(path, {"value": 1})
    before = path.read_bytes()
    with (
        pytest.raises(OSError, match="Interrupted"),
        atomic_output(path) as staging,
    ):
        staging.write_text("partial")
        raise OSError("Interrupted")
    assert path.read_bytes() == before
    with pytest.raises(ValueError):
        atomic_json(path, {"value": float("nan")})
    assert path.read_bytes() == before
    assert not list(tmp_path.glob(".partial-*"))
    atomic_text(path, "complete")
    assert path.read_text() == "complete"


@pytest.mark.parametrize("kind", ["file", "directory"])
def test_atomic_writes_reject_symbolic_links(tmp_path: Path, kind: str) -> None:
    """Protect unrelated outputs from both final and intermediate symlink paths."""
    protected = tmp_path / "protected"
    protected.mkdir()
    (protected / "result.txt").write_text("protected")
    destination = tmp_path / "output"
    if kind == "file":
        destination.symlink_to(protected / "result.txt")
    else:
        destination.symlink_to(protected, target_is_directory=True)
        destination /= "result.txt"
    with pytest.raises(ValueError, match="symbolic links"):
        atomic_text(destination, "unrelated result")
    assert (protected / "result.txt").read_text() == "protected"
