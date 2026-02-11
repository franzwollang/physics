from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

from simulations.vacuum.vacuum_geometry_simulation import SimulationParams, run_simulation_with_graph

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from common_io import ensure_dir, setup_logging, write_json, write_run_manifest  # noqa: E402


def compute_graph_slope(evals: np.ndarray) -> tuple[float, float, float]:
    evals_sorted = np.sort(evals)
    if evals_sorted.size <= 1:
        return float("nan"), float("nan"), float("nan")
    k = np.sqrt(evals_sorted[1:])
    i_lo = int(math.floor(0.2 * len(k)))
    i_hi = int(math.floor(0.8 * len(k)))
    i_hi = max(i_hi, i_lo + 1)
    k_lo = k[i_lo]
    k_hi = k[i_hi]

    band_mask = (k >= k_lo) & (k <= k_hi)
    k_band = k[band_mask]
    if k_band.size < 2:
        return float("nan"), float(k_lo), float(k_hi)

    hist, edges = np.histogram(k_band, bins=30, range=(k_lo, k_hi), density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    mask = (hist > 0) & np.isfinite(hist) & (centers > 0)
    if np.sum(mask) < 2:
        return float("nan"), float(k_lo), float(k_hi)

    x = np.log(centers[mask])
    y = np.log(hist[mask] / (centers[mask] ** 2))
    slope, _ = np.polyfit(x, y, 1)
    n_s_graph = -float(slope)
    return n_s_graph, float(k_lo), float(k_hi)


def _flat_to_indices(flat: int, n_beta: int, n_seeds: int) -> tuple[int, int, int]:
    per_alpha = n_beta * n_seeds
    i_alpha = flat // per_alpha
    rem = flat % per_alpha
    i_beta = rem // n_seeds
    i_seed = rem % n_seeds
    return i_alpha, i_beta, i_seed


def _format_eta(elapsed: float, progress: int, total: int) -> str:
    if progress <= 0:
        return "ETA: unknown"
    rate = elapsed / progress
    remaining = max(total - progress, 0)
    eta = rate * remaining
    return f"ETA: {eta:.1f}s"


class RunTimeout(Exception):
    pass


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage 3: Graph origin check.")
    parser.add_argument("--stage2-json", type=Path, default=Path("stage2_fit.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--checkpoint-dir", type=Path, default=Path("checkpoints/stage3"))
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--max-wall-seconds", type=int, default=1800)
    parser.add_argument("--log-dir", type=Path, default=Path("logs"))
    args = parser.parse_args()

    output_dir = args.output_dir
    ensure_dir(output_dir)
    checkpoint_dir = args.checkpoint_dir
    ensure_dir(checkpoint_dir)

    logger, log_path = setup_logging(args.log_dir, "stage3")
    logger.info("Stage 3 starting. Logs at %s", log_path)

    stage2 = json.loads(args.stage2_json.read_text())
    n_s_ci = stage2["bootstrap"]["n_s_ci_95"]
    n_s_lo, n_s_hi = float(n_s_ci[0]), float(n_s_ci[1])

    alpha_grid = np.array([0.5, 1.0, 2.0], dtype=float)
    beta_grid = np.array([10.0, 30.0, 100.0], dtype=float)
    seeds = np.arange(10, dtype=int)

    shape = (alpha_grid.size, beta_grid.size, seeds.size)
    kept_mask = np.zeros(shape, dtype=bool)
    final_d_s = np.full(shape, np.nan)
    fiedler = np.full(shape, np.nan)
    n_s_graph = np.full(shape, np.nan)
    k_fit_band_lo = np.full(shape, np.nan)
    k_fit_band_hi = np.full(shape, np.nan)

    vacuum_runs = []

    checkpoint_path = checkpoint_dir / "stage3_checkpoint.npz"
    start_flat = 0
    if args.resume and checkpoint_path.exists() and not args.force_recompute:
        data = np.load(checkpoint_path, allow_pickle=True)
        start_flat = int(data["next_flat"])
        kept_mask = data["kept_mask"]
        final_d_s = data["final_d_s"]
        fiedler = data["fiedler_eigenvalue"]
        n_s_graph = data["n_s_graph"]
        k_fit_band_lo = data["k_fit_band_lo"]
        k_fit_band_hi = data["k_fit_band_hi"]
        logger.info("Resuming from checkpoint %s at flat index %d", checkpoint_path, start_flat)
    elif args.resume and not checkpoint_path.exists():
        raise SystemExit(f"Checkpoint not found at {checkpoint_path}")

    total = alpha_grid.size * beta_grid.size * seeds.size
    start_time = time.monotonic()
    last_log = {"time": start_time}
    last_checkpoint = {"time": start_time}
    log_interval = 30.0
    checkpoint_interval = 60.0

    current_flat = {"value": start_flat}

    def write_checkpoint(next_flat: int) -> None:
        np.savez(
            checkpoint_path,
            next_flat=next_flat,
            sweep_alpha_grid=alpha_grid,
            sweep_beta_grid=beta_grid,
            sweep_seeds=seeds,
            kept_mask=kept_mask,
            final_d_s=final_d_s,
            fiedler_eigenvalue=fiedler,
            n_s_graph=n_s_graph,
            k_fit_band_lo=k_fit_band_lo,
            k_fit_band_hi=k_fit_band_hi,
        )

    def check_progress() -> None:
        now = time.monotonic()
        progress = current_flat["value"] - start_flat
        if now - last_log["time"] >= log_interval:
            elapsed = now - start_time
            eta = _format_eta(elapsed, max(progress, 1), total - start_flat)
            i_a, i_b, i_s = _flat_to_indices(current_flat["value"], beta_grid.size, seeds.size)
            logger.info(
                "Progress %d/%d (alpha_idx=%d, beta_idx=%d, seed_idx=%d). Elapsed %.1fs. %s",
                progress,
                total - start_flat,
                i_a,
                i_b,
                i_s,
                elapsed,
                eta,
            )
            last_log["time"] = now

        if now - last_checkpoint["time"] >= checkpoint_interval:
            write_checkpoint(current_flat["value"])
            logger.info("Checkpoint saved at %s (next_flat=%d)", checkpoint_path, current_flat["value"])
            last_checkpoint["time"] = now

        if now - start_time >= args.max_wall_seconds:
            write_checkpoint(current_flat["value"])
            logger.error("RUN_TIMEOUT_CHECKPOINT_WRITTEN")
            print("RUN_TIMEOUT_CHECKPOINT_WRITTEN")
            raise RunTimeout()

    try:
        for flat in range(start_flat, total):
            current_flat["value"] = flat
            i_a, i_b, i_s = _flat_to_indices(flat, beta_grid.size, seeds.size)
            alpha = float(alpha_grid[i_a])
            beta = float(beta_grid[i_b])
            seed = int(seeds[i_s])
            check_progress()

            params = SimulationParams(
                n=200,
                alpha=float(alpha),
                temperature=0.1,
                beta=float(beta),
                d_star=3.0,
                tau_star=1.0,
                steps=2000,
                eta=0.01,
                initial_graph="ErdosRenyi_p0.3",
                seed=int(seed),
            )
            metadata, _coords, _j, evals = run_simulation_with_graph(params)
            outputs = metadata["outputs"]
            ds_final = float(outputs["final_d_s"])
            fiedler_val = float(outputs["fiedler_eigenvalue"])

            final_d_s[i_a, i_b, i_s] = ds_final
            fiedler[i_a, i_b, i_s] = fiedler_val

            keep = (2.9 <= ds_final <= 3.1) and (fiedler_val > 1e-6)
            kept_mask[i_a, i_b, i_s] = keep

            if keep:
                slope, k_lo, k_hi = compute_graph_slope(evals)
                n_s_graph[i_a, i_b, i_s] = slope
                k_fit_band_lo[i_a, i_b, i_s] = k_lo
                k_fit_band_hi[i_a, i_b, i_s] = k_hi

            if flat == total - 1:
                write_checkpoint(flat + 1)
    except RunTimeout:
        sys.exit(1)

    write_checkpoint(total)

    vacuum_runs = []
    for i_a, alpha in enumerate(alpha_grid):
        for i_b, beta in enumerate(beta_grid):
            for i_s, seed in enumerate(seeds):
                if not kept_mask[i_a, i_b, i_s]:
                    continue
                ns_val = float(n_s_graph[i_a, i_b, i_s])
                k_lo = float(k_fit_band_lo[i_a, i_b, i_s])
                k_hi = float(k_fit_band_hi[i_a, i_b, i_s])
                vacuum_runs.append(
                    {
                        "alpha": float(alpha),
                        "T": 0.1,
                        "beta": float(beta),
                        "d_s": float(final_d_s[i_a, i_b, i_s]),
                        "fiedler_eigenvalue": float(fiedler[i_a, i_b, i_s]),
                        "n_s_graph": ns_val,
                        "k_fit_band": [k_lo, k_hi],
                    }
                )

    match = False
    for run in vacuum_runs:
        ns_val = run["n_s_graph"]
        if np.isfinite(ns_val) and n_s_lo <= ns_val <= n_s_hi:
            match = True
            break

    stage3_json = {
        "stage": 3,
        "vacuum_graph_source": "vacuum_geometry_simulation.py (fixed sweep)",
        "graph_phase_model": "OU: dphi/dt = -L phi + sqrt(2) xi, so P_eff(k)=rho(k)/k^2",
        "sweep": {
            "N": 200,
            "steps": 2000,
            "eta": 0.01,
            "tau_star": 1.0,
            "d_star": 3.0,
            "initial_graph": "ErdosRenyi_p0.3",
            "temperature": 0.1,
            "alpha_grid": alpha_grid.tolist(),
            "beta_grid": beta_grid.tolist(),
            "seeds": seeds.tolist(),
            "keep_filters": {
                "final_d_s": [2.9, 3.1],
                "fiedler_eigenvalue_min": 1e-6,
            },
        },
        "vacuum_runs": vacuum_runs,
        "stage2_target": {
            "n_s_ci_95": [n_s_lo, n_s_hi],
        },
        "match": match,
        "seed": int(args.seed),
    }

    write_json(output_dir / "stage3_graph_check.json", stage3_json)

    np.savez(
        output_dir / "raw_aggregates.npz",
        sweep_alpha_grid=alpha_grid,
        sweep_beta_grid=beta_grid,
        sweep_seeds=seeds,
        kept_mask=kept_mask,
        final_d_s=final_d_s,
        fiedler_eigenvalue=fiedler,
        n_s_graph=n_s_graph,
        k_fit_band_lo=k_fit_band_lo,
        k_fit_band_hi=k_fit_band_hi,
    )

    write_run_manifest(output_dir, Path.cwd(), sys.argv)

    if not match:
        (output_dir / "stage3_failure_note.txt").write_text("STAGE3_FAIL_GRAPH_SPECTRUM_MISMATCH")


if __name__ == "__main__":
    main()
