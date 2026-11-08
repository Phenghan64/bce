"""Run the official OpenFold3 CLI and measure this process on macOS.

Usage: ./openfold3/run.sh python openfold3/benchmark.py REPORT.json predict ...
Memory samples are overlapping measures; do not sum RSS and Metal memory.
"""
import json
import resource
import sys
import threading
import time
from importlib.metadata import version
from pathlib import Path

import psutil
import torch


def main():
    report_path = Path(sys.argv[1]).resolve()
    report_path.parent.mkdir(parents=True, exist_ok=True)
    if not torch.backends.mps.is_available():
        raise RuntimeError("Apple MPS is not available; refusing a silent CPU run")
    torch.set_num_threads(6)
    torch.mps.set_per_process_memory_fraction(0.7)
    start = time.perf_counter()
    process = psutil.Process()
    swap_before = psutil.swap_memory().used
    samples = []
    done = threading.Event()

    def sample():
        while not done.is_set():
            samples.append({
                "seconds": round(time.perf_counter() - start, 3),
                "process_rss_bytes": process.memory_info().rss,
                "mps_tensor_bytes": torch.mps.current_allocated_memory(),
                "mps_driver_bytes": torch.mps.driver_allocated_memory(),
                "system_swap_used_bytes": psutil.swap_memory().used,
            })
            done.wait(1)

    monitor = threading.Thread(target=sample, daemon=True)
    monitor.start()
    success = False
    error = None
    try:
        from openfold3.run_openfold import cli
        cli.main(args=sys.argv[2:], standalone_mode=False)
        torch.mps.synchronize()
        success = True
    except BaseException as exc:
        error = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        done.set()
        monitor.join()
        elapsed = time.perf_counter() - start
        report = {
            "cli_completed": success,
            "error": error,
            "elapsed_seconds": round(elapsed, 3),
            "timing_scope": "CLI import and execution; excludes Python/PyTorch startup and prior MSA/download stages",
            "openfold3_version": version("openfold3"),
            "torch_version": torch.__version__,
            "mps_available": torch.backends.mps.is_available(),
            "mps_memory_fraction": 0.7,
            "mps_recommended_max_bytes": torch.mps.recommended_max_memory(),
            "process_peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "sampled_peak_mps_driver_bytes": max(s["mps_driver_bytes"] for s in samples),
            "system_swap_before_bytes": swap_before,
            "system_swap_after_bytes": psutil.swap_memory().used,
            "memory_note": "RSS and MPS driver allocations overlap on unified memory; system swap includes all apps",
            "samples": samples,
        }
        report_path.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps({k: v for k, v in report.items() if k != "samples"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
