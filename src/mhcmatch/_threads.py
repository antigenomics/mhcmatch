"""One explicit budget for native searches; serial by default for composable library calls."""
import operator
import os


def resolve_threads(threads, tasks):
    """Validate a native budget, capping it to work; zero opts into allocated CPUs.

    Automatic mode respects process affinity and SLURM's per-task allocation. It never passes
    zero to a native library, where it could mean the entire host despite an outer worker pool.
    Explicit positive budgets are the caller's allocation, shared by all searches in that call.
    """
    threads = operator.index(threads)
    if threads < 0:
        raise ValueError("threads must be non-negative (0 = allocated CPUs)")
    if threads == 0:
        count = getattr(os, "process_cpu_count", os.cpu_count)() or 1
        if hasattr(os, "sched_getaffinity"):
            count = min(count, len(os.sched_getaffinity(0)))
        allocation = os.environ.get("SLURM_CPUS_PER_TASK")
        if allocation:
            try:
                count = min(count, max(1, int(allocation)))
            except ValueError:
                raise ValueError("SLURM_CPUS_PER_TASK must be an integer") from None
        threads = count
    return max(1, min(threads, tasks))
