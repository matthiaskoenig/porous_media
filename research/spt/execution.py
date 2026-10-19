"""Execute independent SPT jobs with the same serial and spawned-worker behavior."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from research.spt.cache import ProcessingResult


def execute_jobs[Job, Result](
    jobs: Sequence[Job],
    action: Callable[[Job], Result],
    *,
    workers: int,
    failed: Callable[[Job, Exception], Result],
) -> Iterator[tuple[Job, Result]]:
    """Yield completed jobs, isolating failures and avoiding unnecessary processes.

    The action must be picklable for spawned workers. Failure conversion happens
    in the coordinator. Empty batches never start a pool; single jobs run in the
    caller. Consumers sort their results independently of completion order.
    """
    if workers < 1:
        raise ValueError("workers must be positive")
    if workers == 1 or len(jobs) <= 1:
        for job in jobs:
            try:
                result = action(job)
            except Exception as error:
                result = failed(job, error)
            yield job, result
        return
    with ProcessPoolExecutor(
        max_workers=min(workers, len(jobs)), mp_context=get_context("spawn")
    ) as executor:
        futures = {executor.submit(action, job): job for job in jobs}
        for future in as_completed(futures):
            job = futures[future]
            try:
                result = future.result()
            except Exception as error:
                result = failed(job, error)
            yield job, result


def print_result(result: ProcessingResult, *, prefix: str = "") -> None:
    """Print one completed workflow result consistently across coordinators."""
    print(
        f"{prefix}{result.study}/{result.run_id}: {result.status} "
        f"({result.elapsed_s:.2f}s) {result.message}",
        flush=True,
    )
