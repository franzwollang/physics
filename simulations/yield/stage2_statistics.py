from __future__ import annotations

import argparse
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from common_grid import generate_grf_field, spectral_gradient  # noqa: E402
from common_io import (  # noqa: E402
    ensure_dir,
    setup_logging,
    stable_seed,
    write_json,
    write_run_manifest,
)


@dataclass
class SelectionMap:
    r_values: np.ndarray
    sigma_values: np.ndarray
    p_grid: np.ndarray


def load_selection_map(path: Path) -> SelectionMap:
    data = json.loads(path.read_text())
    if "grid" in data:
        r_values = np.array(data["grid"]["R_values"], dtype=float)
        sigma_values = np.array(data["grid"]["sigma_phi_values"], dtype=float)
    else:
        r_min, r_max = data["param_R"]
        s_min, s_max = data["param_sigma"]
        r_res, s_res = data["map_resolution"]
        r_values = np.logspace(np.log10(r_min), np.log10(r_max), int(r_res))
        sigma_values = np.linspace(s_min, s_max, int(s_res))
    p_grid = np.array(data["transition_matrix"], dtype=float)
    return SelectionMap(r_values=r_values, sigma_values=sigma_values, p_grid=p_grid)


def interpolate_selection(
    r: np.ndarray,
    sigma: np.ndarray,
    selection: SelectionMap,
) -> np.ndarray:
    r_vals = selection.r_values
    s_vals = selection.sigma_values
    p = selection.p_grid

    r = np.clip(r, r_vals[0], r_vals[-1])
    sigma = np.clip(sigma, s_vals[0], s_vals[-1])

    i = np.searchsorted(r_vals, r, side="right") - 1
    j = np.searchsorted(s_vals, sigma, side="right") - 1
    i = np.clip(i, 0, r_vals.size - 2)
    j = np.clip(j, 0, s_vals.size - 2)

    r0 = r_vals[i]
    r1 = r_vals[i + 1]
    s0 = s_vals[j]
    s1 = s_vals[j + 1]

    t = (r - r0) / (r1 - r0 + 1e-12)
    u = (sigma - s0) / (s1 - s0 + 1e-12)

    p00 = p[i, j]
    p10 = p[i + 1, j]
    p01 = p[i, j + 1]
    p11 = p[i + 1, j + 1]

    p_interp = (1 - t) * (1 - u) * p00 + t * (1 - u) * p10 + (1 - t) * u * p01 + t * u * p11
    return np.clip(p_interp, 0.0, 1.0)


def find_peaks(u: np.ndarray) -> np.ndarray:
    neighbors = [
        np.roll(u, 1, axis=0),
        np.roll(u, -1, axis=0),
        np.roll(u, 1, axis=1),
        np.roll(u, -1, axis=1),
        np.roll(np.roll(u, 1, axis=0), 1, axis=1),
        np.roll(np.roll(u, 1, axis=0), -1, axis=1),
        np.roll(np.roll(u, -1, axis=0), 1, axis=1),
        np.roll(np.roll(u, -1, axis=0), -1, axis=1),
    ]
    is_peak = np.ones_like(u, dtype=bool)
    for nbr in neighbors:
        is_peak &= u > nbr
    peak_indices = np.argwhere(is_peak)
    return peak_indices


def non_max_suppression(
    indices: np.ndarray,
    values: np.ndarray,
    grid_size: int,
    max_keep: int = 2000,
) -> np.ndarray:
    order = np.argsort(values)[::-1]
    candidates = indices[order]
    n = grid_size
    suppressed = np.zeros((n, n), dtype=bool) if n > 0 else None
    kept: list[tuple[int, int]] = []

    for idx in candidates:
        i, j = int(idx[0]), int(idx[1])
        if suppressed is not None and suppressed[i, j]:
            continue
        kept.append((i, j))
        if suppressed is not None:
            i0 = max(0, i - 3)
            i1 = min(n - 1, i + 3)
            j0 = max(0, j - 3)
            j1 = min(n - 1, j + 3)
            suppressed[i0 : i1 + 1, j0 : j1 + 1] = True
        if len(kept) >= max_keep:
            break
    return np.array(kept, dtype=int)


def select_peak_centers(u: np.ndarray) -> np.ndarray:
    peaks = find_peaks(u)
    if peaks.size == 0:
        return peaks
    mu = float(np.mean(u))
    sigma = float(np.std(u))
    values = u[peaks[:, 0], peaks[:, 1]]

    threshold = mu + 2.0 * sigma
    mask = values > threshold
    selected = peaks[mask]

    if selected.shape[0] < 2000:
        threshold = mu + 1.5 * sigma
        mask = values > threshold
        selected = peaks[mask]

    if selected.shape[0] < 2000:
        flat = u.ravel()
        top_idx = np.argpartition(flat, -2000)[-2000:]
        top_idx = top_idx[np.argsort(flat[top_idx])[::-1]]
        candidates = np.column_stack(np.unravel_index(top_idx, u.shape))
        selected = non_max_suppression(candidates, flat[top_idx], u.shape[0], max_keep=2000)

    if selected.shape[0] > 5000:
        selected_values = u[selected[:, 0], selected[:, 1]]
        order = np.argsort(selected_values)[::-1][:5000]
        selected = selected[order]

    return selected


def compute_shell_rms(
    grad_phi_sq: np.ndarray,
    centers: np.ndarray,
    radii: np.ndarray,
    dx: float,
    shell_width: float,
    check_fn=None,
    check_stride: int = 200,
) -> np.ndarray:
    n = grad_phi_sq.shape[0]
    max_radius = 8.0 + shell_width
    max_pix = int(np.ceil(max_radius / dx))
    offsets = np.arange(-max_pix, max_pix + 1)
    offset_x = offsets[None, :] * dx
    offset_y = offsets[:, None] * dx
    offset_r = np.sqrt(offset_x**2 + offset_y**2)

    sigma_phi = np.zeros(centers.shape[0], dtype=float)

    for idx, (i, j) in enumerate(centers):
        if check_fn is not None and idx % check_stride == 0:
            check_fn()
        radius = radii[idx]
        mask = (offset_r >= radius) & (offset_r <= radius + shell_width)
        rows = (i + offsets) % n
        cols = (j + offsets) % n
        patch = grad_phi_sq[np.ix_(rows, cols)]
        values = patch[mask]
        if values.size == 0:
            sigma_phi[idx] = 0.0
        else:
            sigma_phi[idx] = float(np.sqrt(np.mean(values)))
    return sigma_phi


def _flat_to_indices(flat: int, n_sigma: int, n_realizations: int) -> tuple[int, int, int]:
    per_ns = n_sigma * n_realizations
    i_ns = flat // per_ns
    rem = flat % per_ns
    j_s = rem // n_realizations
    r_idx = rem % n_realizations
    return i_ns, j_s, r_idx


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
    parser = argparse.ArgumentParser(description="Stage 2: Rogue Wave fit.")
    parser.add_argument("--stage1-json", type=Path, default=Path("stage1_selection.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--checkpoint-dir", type=Path, default=Path("checkpoints/stage2"))
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--max-wall-seconds", type=int, default=1800)
    parser.add_argument("--log-dir", type=Path, default=Path("logs"))
    args = parser.parse_args()

    output_dir = args.output_dir
    ensure_dir(output_dir)
    checkpoint_dir = args.checkpoint_dir
    ensure_dir(checkpoint_dir)

    logger, log_path = setup_logging(args.log_dir, "stage2")
    logger.info("Stage 2 starting. Logs at %s", log_path)

    selection = load_selection_map(args.stage1_json)

    # Fixed parameters
    l = 80.0
    n = 1024
    dx = l / n
    k_min = 2.0 * np.pi / l
    k_max = np.pi / dx

    n_s_grid = np.round(np.arange(0.5, 4.0 + 1e-9, 0.05), 2)
    s_sigma_grid = np.logspace(np.log10(0.25), np.log10(4.0), 25)

    num_realizations = 200
    base_seed = args.seed

    shape = (n_s_grid.size, s_sigma_grid.size, num_realizations)
    E_B_real = np.zeros(shape, dtype=float)
    E_DM_real = np.zeros(shape, dtype=float)
    N_B_real = np.zeros(shape, dtype=float)
    N_DM_real = np.zeros(shape, dtype=float)

    checkpoint_path = checkpoint_dir / "stage2_checkpoint.npz"
    start_flat = 0
    if args.resume and checkpoint_path.exists() and not args.force_recompute:
        data = np.load(checkpoint_path, allow_pickle=True)
        start_flat = int(data["next_flat"])
        E_B_real = data["E_B_real"]
        E_DM_real = data["E_DM_real"]
        N_B_real = data["N_B_real"]
        N_DM_real = data["N_DM_real"]
        logger.info("Resuming from checkpoint %s at flat index %d", checkpoint_path, start_flat)
    elif args.resume and not checkpoint_path.exists():
        raise SystemExit(f"Checkpoint not found at {checkpoint_path}")

    total = n_s_grid.size * s_sigma_grid.size * num_realizations
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
            n_s_grid=n_s_grid,
            s_sigma_grid=s_sigma_grid,
            E_B_real=E_B_real,
            E_DM_real=E_DM_real,
            N_B_real=N_B_real,
            N_DM_real=N_DM_real,
        )

    def check_progress() -> None:
        now = time.monotonic()
        progress = current_flat["value"] - start_flat
        if now - last_log["time"] >= log_interval:
            elapsed = now - start_time
            eta = _format_eta(elapsed, max(progress, 1), total - start_flat)
            i_ns, j_s, r_idx = _flat_to_indices(
                current_flat["value"], s_sigma_grid.size, num_realizations
            )
            logger.info(
                "Progress %d/%d (n_s_idx=%d, s_sigma_idx=%d, realization_idx=%d). Elapsed %.1fs. %s",
                progress,
                total - start_flat,
                i_ns,
                j_s,
                r_idx,
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
            i_ns, j_s, r_idx = _flat_to_indices(flat, s_sigma_grid.size, num_realizations)
            n_s = float(n_s_grid[i_ns])
            s_sigma = float(s_sigma_grid[j_s])
            check_progress()

            phi_seed = stable_seed(base_seed, i_ns, j_s, r_idx, 0)
            u_seed = stable_seed(base_seed, i_ns, j_s, r_idx, 1)
            rng_phi = np.random.default_rng(phi_seed)
            rng_u = np.random.default_rng(u_seed)

            def power_phi(k: np.ndarray) -> np.ndarray:
                return np.where((k >= k_min) & (k <= k_max), k ** (-n_s), 0.0)

            phi = generate_grf_field(rng_phi, n, l, power_phi)
            grad_phi_x, grad_phi_y = spectral_gradient(phi, l)
            g = float(np.sqrt(np.mean(grad_phi_x**2 + grad_phi_y**2)))
            phi = phi / (g + 1e-12)
            phi = phi * s_sigma
            grad_phi_x, grad_phi_y = spectral_gradient(phi, l)
            grad_phi_sq = grad_phi_x**2 + grad_phi_y**2

            def power_u(k: np.ndarray) -> np.ndarray:
                return np.where((k >= k_min) & (k <= k_max), 1.0, 0.0)

            u = generate_grf_field(rng_u, n, l, power_u)

            centers = select_peak_centers(u)
            n_b = 0.0
            n_dm = 0.0
            e_b = 0.0
            e_dm = 0.0

            if centers.size != 0:
                lap_u = (
                    np.roll(u, 1, axis=0)
                    + np.roll(u, -1, axis=0)
                    + np.roll(u, 1, axis=1)
                    + np.roll(u, -1, axis=1)
                    - 4.0 * u
                ) / (dx * dx)

                u_vals = u[centers[:, 0], centers[:, 1]]
                lap_vals = lap_u[centers[:, 0], centers[:, 1]]
                kappa_vals = np.maximum(1e-6, -(lap_vals) / (u_vals + 1e-12))
                radii = np.clip(kappa_vals ** -0.5, 0.5, 8.0)

                sigma_phi = compute_shell_rms(
                    grad_phi_sq, centers, radii, dx, 2.0, check_fn=check_progress
                )
                p_i = interpolate_selection(radii, sigma_phi, selection)

                n_b = float(np.sum(p_i))
                n_dm = float(np.sum(1.0 - p_i))
                weights = radii**3
                e_b = float(np.sum(weights * p_i))
                e_dm = float(np.sum(weights * (1.0 - p_i)))

            E_B_real[i_ns, j_s, r_idx] = e_b
            E_DM_real[i_ns, j_s, r_idx] = e_dm
            N_B_real[i_ns, j_s, r_idx] = n_b
            N_DM_real[i_ns, j_s, r_idx] = n_dm

            if flat == total - 1:
                write_checkpoint(flat + 1)
    except RunTimeout:
        sys.exit(1)

    write_checkpoint(total)

    E_B = np.sum(E_B_real, axis=2)
    E_DM = np.sum(E_DM_real, axis=2)
    N_B = np.sum(N_B_real, axis=2)
    N_DM = np.sum(N_DM_real, axis=2)

    eps = 1e-12
    ratio_energy = E_DM / (E_B + eps)
    ratio_count = N_DM / (N_B + eps)
    objective = (np.log(ratio_energy + eps) - math.log(5.3)) ** 2

    best_idx = np.unravel_index(np.argmin(objective), objective.shape)
    best_ns = float(n_s_grid[best_idx[0]])
    best_s_sigma = float(s_sigma_grid[best_idx[1]])
    best_ratio_energy = float(ratio_energy[best_idx])
    best_ratio_count = float(ratio_count[best_idx])

    # Bootstrap
    rng_boot = np.random.default_rng(stable_seed(base_seed, 9999))
    bootstrap_best_ns = np.zeros(500)
    bootstrap_best_s_sigma = np.zeros(500)
    bootstrap_ratio_energy = np.zeros(500)

    for b in range(500):
        resample = rng_boot.integers(0, num_realizations, size=num_realizations)
        E_B_boot = np.sum(E_B_real[:, :, resample], axis=2)
        E_DM_boot = np.sum(E_DM_real[:, :, resample], axis=2)
        ratio_boot = E_DM_boot / (E_B_boot + eps)
        obj_boot = (np.log(ratio_boot + eps) - math.log(5.3)) ** 2
        idx = np.unravel_index(np.argmin(obj_boot), obj_boot.shape)
        bootstrap_best_ns[b] = n_s_grid[idx[0]]
        bootstrap_best_s_sigma[b] = s_sigma_grid[idx[1]]
        bootstrap_ratio_energy[b] = ratio_boot[idx]

    def ci_bounds(values: np.ndarray, low: float, high: float) -> list[float]:
        return [float(np.quantile(values, low)), float(np.quantile(values, high))]

    stage2_json = {
        "stage": 2,
        "spectrum_family": "P_phi(k)=k^{-n_s}*1_[kmin,kmax], normalized so RMS(|∇phi|)=1 then scaled by s_sigma",
        "best_fit": {
            "n_s": best_ns,
            "s_sigma": best_s_sigma,
            "k_min": float(k_min),
            "k_max": float(k_max),
        },
        "yield_ratio_energy": best_ratio_energy,
        "yield_ratio_count": best_ratio_count,
        "bootstrap": {
            "resamples": 500,
            "n_s_ci_95": ci_bounds(bootstrap_best_ns, 0.025, 0.975),
            "s_sigma_ci_95": ci_bounds(bootstrap_best_s_sigma, 0.025, 0.975),
            "yield_ratio_energy_ci_95": ci_bounds(bootstrap_ratio_energy, 0.025, 0.975),
            "n_s_ci_68": ci_bounds(bootstrap_best_ns, 0.16, 0.84),
            "s_sigma_ci_68": ci_bounds(bootstrap_best_s_sigma, 0.16, 0.84),
            "yield_ratio_energy_ci_68": ci_bounds(bootstrap_ratio_energy, 0.16, 0.84),
        },
        "protocol": {
            "domain_L": l,
            "grid_N": n,
            "num_realizations": num_realizations,
            "peak_threshold_sigma": 2.0,
            "min_peaks": 2000,
            "max_peaks": 5000,
            "radius_clip": [0.5, 8.0],
            "shell_width": 2.0,
            "energy_weight": "R^3",
            "n_s_grid": [0.5, 4.0, 0.05],
            "s_sigma_grid": [0.25, 4.0, 25, "logspace"],
        },
        "seed": int(base_seed),
    }

    write_json(output_dir / "stage2_fit.json", stage2_json)

    np.savez(
        output_dir / "raw_aggregates.npz",
        n_s_grid=n_s_grid,
        s_sigma_grid=s_sigma_grid,
        objective=objective,
        E_ratio_energy=ratio_energy,
        E_ratio_count=ratio_count,
        E_B=E_B,
        E_DM=E_DM,
        N_B=N_B,
        N_DM=N_DM,
        bootstrap_best_ns=bootstrap_best_ns,
        bootstrap_best_s_sigma=bootstrap_best_s_sigma,
        bootstrap_yield_ratio_energy=bootstrap_ratio_energy,
    )

    write_run_manifest(output_dir, Path.cwd(), sys.argv)

    fail_note = None
    if not (4.5 <= best_ratio_energy <= 6.5):
        fail_note = "STAGE2_FAIL_FIXED_MODEL_INCONSISTENT"

    if fail_note is not None:
        (output_dir / "stage2_failure_note.txt").write_text(fail_note)


if __name__ == "__main__":
    main()
