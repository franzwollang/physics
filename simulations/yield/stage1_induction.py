from __future__ import annotations

import argparse
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from common_io import (  # noqa: E402
    ensure_dir,
    setup_logging,
    stable_seed,
    write_json,
    write_run_manifest,
)
from common_topology import (  # noqa: E402
    sample_complex_from_polar_on_circle,
    winding_number_from_samples,
)


@dataclass
class SGLEWorkspace:
    lap_w: np.ndarray
    lap_phi: np.ndarray
    grad_w_x: np.ndarray
    grad_w_y: np.ndarray
    grad_phi_x: np.ndarray
    grad_phi_y: np.ndarray
    grad_phi_sq: np.ndarray
    coupling: np.ndarray
    tmp: np.ndarray
    tmp2: np.ndarray


class RunTimeout(Exception):
    pass


def _wilson_ci(successes: int, trials: int, z: float = 1.96) -> tuple[float, float]:
    if trials <= 0:
        return 0.0, 1.0
    p = successes / trials
    denom = 1.0 + (z * z) / trials
    center = (p + (z * z) / (2.0 * trials)) / denom
    margin = (
        z
        * math.sqrt((p * (1.0 - p) / trials) + (z * z) / (4.0 * trials * trials))
        / denom
    )
    return max(0.0, center - margin), min(1.0, center + margin)


def _boundary_indices_and_s(n: int, l: float) -> tuple[np.ndarray, np.ndarray, float]:
    indices: list[tuple[int, int]] = []
    # bottom edge
    for j in range(n):
        indices.append((0, j))
    # right edge
    for i in range(1, n):
        indices.append((i, n - 1))
    # top edge
    for j in range(n - 2, -1, -1):
        indices.append((n - 1, j))
    # left edge
    for i in range(n - 2, 0, -1):
        indices.append((i, 0))

    boundary = np.array(indices, dtype=int)
    l_b = 4.0 * l
    s_values = np.linspace(0.0, l_b, boundary.shape[0], endpoint=False)
    return boundary, s_values, l_b


def _eta_from_coeffs(a_m: np.ndarray, b_m: np.ndarray, s_values: np.ndarray, l_b: float) -> np.ndarray:
    m = np.arange(1, a_m.size + 1, dtype=float)[:, None]
    angles = 2.0 * np.pi * m * s_values[None, :] / l_b
    eta = (a_m[:, None] * np.cos(angles) + b_m[:, None] * np.sin(angles)).sum(axis=0)
    return eta.astype(np.float64, copy=False)


def _apply_boundary_phi(phi: np.ndarray, boundary: np.ndarray, values: np.ndarray) -> None:
    phi[boundary[:, 0], boundary[:, 1]] = values


def _apply_boundary_w(w: np.ndarray) -> None:
    w[0, :] = 1.0
    w[-1, :] = 1.0
    w[:, 0] = 1.0
    w[:, -1] = 1.0


def _compute_sigma_phi(
    phi: np.ndarray,
    shell_mask: np.ndarray,
    dx: float,
    workspace: SGLEWorkspace,
) -> float:
    inv2dx = 1.0 / (2.0 * dx)
    grad_phi_x = workspace.grad_phi_x
    grad_phi_y = workspace.grad_phi_y

    # Only fill where needed; boundaries are irrelevant for shell stats.
    grad_phi_x[:, 1:-1] = (phi[:, 2:] - phi[:, :-2]) * inv2dx
    grad_phi_y[1:-1, :] = (phi[2:, :] - phi[:-2, :]) * inv2dx

    values = grad_phi_x[shell_mask] ** 2 + grad_phi_y[shell_mask] ** 2
    if values.size == 0:
        return 0.0
    return float(np.sqrt(np.mean(values)))


def _phase_only_step(
    phi: np.ndarray,
    boundary: np.ndarray,
    boundary_phi: np.ndarray,
    dt: float,
    inv_dx2: float,
    kappa: float,
    workspace: SGLEWorkspace,
) -> None:
    lap_phi_int = workspace.lap_phi[1:-1, 1:-1]
    phi_int = phi[1:-1, 1:-1]

    lap_phi_int[...] = (
        phi[2:, 1:-1]
        + phi[:-2, 1:-1]
        + phi[1:-1, 2:]
        + phi[1:-1, :-2]
        - 4.0 * phi_int
    ) * inv_dx2

    phi_int += dt * kappa * lap_phi_int
    _apply_boundary_phi(phi, boundary, boundary_phi)


def _coupled_step(
    w: np.ndarray,
    phi: np.ndarray,
    boundary: np.ndarray,
    boundary_phi: np.ndarray,
    dt: float,
    inv_dx2: float,
    dx: float,
    alpha: float,
    beta: float,
    gamma: float,
    kappa: float,
    eps_w: float,
    workspace: SGLEWorkspace,
) -> None:
    inv2dx = 1.0 / (2.0 * dx)

    w_int = w[1:-1, 1:-1]
    phi_int = phi[1:-1, 1:-1]

    lap_w_int = workspace.lap_w[1:-1, 1:-1]
    lap_phi_int = workspace.lap_phi[1:-1, 1:-1]

    lap_w_int[...] = (
        w[2:, 1:-1]
        + w[:-2, 1:-1]
        + w[1:-1, 2:]
        + w[1:-1, :-2]
        - 4.0 * w_int
    ) * inv_dx2
    lap_phi_int[...] = (
        phi[2:, 1:-1]
        + phi[:-2, 1:-1]
        + phi[1:-1, 2:]
        + phi[1:-1, :-2]
        - 4.0 * phi_int
    ) * inv_dx2

    grad_w_x = workspace.grad_w_x
    grad_w_y = workspace.grad_w_y
    grad_phi_x = workspace.grad_phi_x
    grad_phi_y = workspace.grad_phi_y

    grad_w_x[:, 1:-1] = (w[:, 2:] - w[:, :-2]) * inv2dx
    grad_w_y[1:-1, :] = (w[2:, :] - w[:-2, :]) * inv2dx
    grad_phi_x[:, 1:-1] = (phi[:, 2:] - phi[:, :-2]) * inv2dx
    grad_phi_y[1:-1, :] = (phi[2:, :] - phi[:-2, :]) * inv2dx

    grad_phi_sq_int = workspace.grad_phi_sq[1:-1, 1:-1]
    coupling_int = workspace.coupling[1:-1, 1:-1]

    # grad_phi_sq
    gx = grad_phi_x[1:-1, 1:-1]
    gy = grad_phi_y[1:-1, 1:-1]
    grad_phi_sq_int[...] = gx * gx
    grad_phi_sq_int += gy * gy

    # coupling = 2 * (∇w · ∇phi) / (w + eps)
    coupling_int[...] = grad_w_x[1:-1, 1:-1] * gx
    coupling_int += grad_w_y[1:-1, 1:-1] * gy
    coupling_int *= 2.0
    coupling_int /= (w_int + eps_w)

    # w update
    tmp_int = workspace.tmp[1:-1, 1:-1]
    tmp_int[...] = lap_w_int
    tmp_int -= w_int * grad_phi_sq_int
    tmp_int *= alpha
    tmp_int += beta * w_int

    # tmp2 = w^3
    tmp2_int = workspace.tmp2[1:-1, 1:-1]
    tmp2_int[...] = w_int * w_int
    tmp2_int *= w_int

    tmp_int -= 2.0 * gamma * tmp2_int
    w_int += dt * tmp_int

    # phi update
    phi_int += dt * kappa * (lap_phi_int + coupling_int)

    _apply_boundary_w(w)
    _apply_boundary_phi(phi, boundary, boundary_phi)


def _format_eta(elapsed: float, progress: float, total: float) -> str:
    if progress <= 0.0:
        return "ETA: unknown"
    rate = elapsed / progress
    remaining = max(total - progress, 0.0)
    eta = rate * remaining
    return f"ETA: {eta:.1f}s"


def _make_workspace(n: int) -> SGLEWorkspace:
    zeros = lambda: np.zeros((n, n), dtype=np.float64)  # noqa: E731
    return SGLEWorkspace(
        lap_w=zeros(),
        lap_phi=zeros(),
        grad_w_x=zeros(),
        grad_w_y=zeros(),
        grad_phi_x=zeros(),
        grad_phi_y=zeros(),
        grad_phi_sq=zeros(),
        coupling=zeros(),
        tmp=zeros(),
        tmp2=zeros(),
    )


def _save_checkpoint(
    path: Path,
    next_flat: int,
    success_counts: np.ndarray,
    trial_counts: np.ndarray,
    sigma_phi_measured_calibration: np.ndarray,
    sigma_phi_measured_start_drive: np.ndarray,
    n_final: np.ndarray,
    min_absPsi_min_over_time: np.ndarray,
    state: dict | None,
    r_values: np.ndarray,
    sigma_phi_values: np.ndarray,
) -> None:
    payload = {
        "next_flat": int(next_flat),
        "R_values": r_values,
        "sigma_phi_values": sigma_phi_values,
        "success_counts": success_counts,
        "trial_counts": trial_counts,
        "sigma_phi_measured_calibration": sigma_phi_measured_calibration,
        "sigma_phi_measured_start_drive": sigma_phi_measured_start_drive,
        "n_final": n_final,
        "min_absPsi_min_over_time": min_absPsi_min_over_time,
    }
    if state is not None:
        payload.update(state)
    np.savez(path, **payload)


def _load_checkpoint(path: Path) -> dict:
    data = np.load(path, allow_pickle=True)
    return {k: data[k] for k in data.keys()}


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage 1: Induction test for primordial yield.")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--checkpoint-dir", type=Path, default=Path("checkpoints/stage1"))
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--max-wall-seconds", type=int, default=1800)
    parser.add_argument("--log-dir", type=Path, default=Path("logs"))
    args = parser.parse_args()

    output_dir = args.output_dir
    ensure_dir(output_dir)
    checkpoint_dir = args.checkpoint_dir
    ensure_dir(checkpoint_dir)

    logger, log_path = setup_logging(args.log_dir, "stage1")
    logger.info("Stage 1 starting. Logs at %s", log_path)

    # Fixed PDE parameters (dimensionless units).
    alpha = 2.0
    beta = 1.0
    gamma = 0.5
    kappa = 1.0
    eps_w = 1e-6

    # Fixed domain and discretization.
    dim = 2
    l = 80.0
    n = 512
    dx = l / n
    inv_dx2 = 1.0 / (dx * dx)

    dt = 5e-4
    t_relax = 20.0
    t_drive = 80.0
    dt_sample = 0.2
    t_persist = 10.0
    t_cal = 5.0

    steps_relax = int(round(t_relax / dt))
    steps_cal = int(round(t_cal / dt))
    steps_drive = int(round(t_drive / dt))
    sample_stride = int(round(dt_sample / dt))

    m_loop = 512
    epsilon_core = 0.2

    r_values = np.logspace(np.log10(0.5), np.log10(8.0), 30)
    sigma_targets = np.linspace(0.0, 3.0, 30)
    trials_per_point = 200

    # Use endpoint=False so dx matches L/N exactly.
    x_min = -l / 2.0
    y_min = -l / 2.0
    x = x_min + dx * np.arange(n)
    y = y_min + dx * np.arange(n)
    xx, yy = np.meshgrid(x, y, indexing="xy")
    r_grid = np.sqrt(xx**2 + yy**2)

    boundary, s_values, l_b = _boundary_indices_and_s(n, l)
    zero_boundary_phi = np.zeros(boundary.shape[0], dtype=np.float64)

    shell_masks = [(r_grid >= r) & (r_grid <= r + 2.0) for r in r_values]
    disk_masks = [r_grid <= (r / 2.0) for r in r_values]

    workspace = _make_workspace(n)

    success_counts = np.zeros((r_values.size, sigma_targets.size), dtype=int)
    trial_counts = np.full((r_values.size, sigma_targets.size), trials_per_point, dtype=int)
    sigma_phi_measured_cal = np.zeros((r_values.size, sigma_targets.size, trials_per_point), dtype=np.float64)
    sigma_phi_measured_start = np.zeros((r_values.size, sigma_targets.size, trials_per_point), dtype=np.float64)
    n_final = np.zeros((r_values.size, sigma_targets.size, trials_per_point), dtype=np.float64)
    min_abs_min = np.zeros((r_values.size, sigma_targets.size, trials_per_point), dtype=np.float64)

    checkpoint_path = checkpoint_dir / "stage1_checkpoint.npz"

    next_flat = 0
    active_state: dict | None = None

    if args.resume and checkpoint_path.exists() and not args.force_recompute:
        ckpt = _load_checkpoint(checkpoint_path)
        next_flat = int(ckpt.get("next_flat", 0))
        success_counts = ckpt.get("success_counts", success_counts)
        trial_counts = ckpt.get("trial_counts", trial_counts)
        sigma_phi_measured_cal = ckpt.get("sigma_phi_measured_calibration", sigma_phi_measured_cal)
        sigma_phi_measured_start = ckpt.get("sigma_phi_measured_start_drive", sigma_phi_measured_start)
        n_final = ckpt.get("n_final", n_final)
        min_abs_min = ckpt.get("min_absPsi_min_over_time", min_abs_min)

        if "active_flat" in ckpt and int(ckpt["active_flat"]) >= 0:
            active_state = {
                "active_flat": int(ckpt["active_flat"]),
                "active_phase": int(ckpt["active_phase"]),
                "active_step": int(ckpt["active_step"]),
                "w": ckpt["w"],
                "phi": ckpt["phi"],
                "eta_values": ckpt.get("eta_values", None),
                "drive_n_series": ckpt.get("drive_n_series", None),
                "drive_min_abs_series": ckpt.get("drive_min_abs_series", None),
                "drive_min_abs_min": float(ckpt.get("drive_min_abs_min", np.inf)),
            }
            next_flat = active_state["active_flat"]
            logger.info(
                "Resuming mid-trial at flat=%d phase=%d step=%d",
                active_state["active_flat"],
                active_state["active_phase"],
                active_state["active_step"],
            )
        else:
            logger.info("Resuming from checkpoint at next_flat=%d", next_flat)
    elif args.resume and not checkpoint_path.exists():
        raise SystemExit(f"Checkpoint not found at {checkpoint_path}")

    # Representative regimes for debug artifacts.
    repr_sigma_indices = {
        "low": 0,
        "mid": sigma_targets.size // 2,
        "high": sigma_targets.size - 1,
    }
    repr_r_index = r_values.size // 2

    total_trials = r_values.size * sigma_targets.size * trials_per_point

    start_time = time.monotonic()
    last_log = start_time
    last_checkpoint = start_time

    log_interval = 30.0
    checkpoint_interval = 60.0

    current = {
        "flat": next_flat,
        "phase": 0,
        "step": 0,
        "trial_fraction": 0.0,
        "completed_trials": 0,
    }

    def compute_trial_fraction(phase: int, step: int) -> float:
        trial_total_steps = steps_relax + steps_cal + steps_drive
        if phase == 0:
            done = step
        elif phase == 1:
            done = steps_relax + step
        else:
            done = steps_relax + steps_cal + step
        return float(done) / float(trial_total_steps)

    def checkpoint_now(flat_for_checkpoint: int, state: dict | None) -> None:
        state_payload = None
        if state is not None:
            state_payload = {
                "active_flat": int(state["active_flat"]),
                "active_phase": int(state["active_phase"]),
                "active_step": int(state["active_step"]),
                "w": state["w"],
                "phi": state["phi"],
                "eta_values": state.get("eta_values", None) if state.get("eta_values", None) is not None else np.array([]),
                "drive_n_series": state.get("drive_n_series", None) if state.get("drive_n_series", None) is not None else np.array([]),
                "drive_min_abs_series": state.get("drive_min_abs_series", None) if state.get("drive_min_abs_series", None) is not None else np.array([]),
                "drive_min_abs_min": float(state.get("drive_min_abs_min", np.inf)),
            }
        _save_checkpoint(
            checkpoint_path,
            next_flat=flat_for_checkpoint,
            success_counts=success_counts,
            trial_counts=trial_counts,
            sigma_phi_measured_calibration=sigma_phi_measured_cal,
            sigma_phi_measured_start_drive=sigma_phi_measured_start,
            n_final=n_final,
            min_absPsi_min_over_time=min_abs_min,
            state=state_payload,
            r_values=r_values,
            sigma_phi_values=sigma_targets,
        )

    def maybe_log_and_checkpoint(state: dict | None) -> None:
        nonlocal last_log, last_checkpoint
        now = time.monotonic()
        elapsed = now - start_time

        flat = current["flat"]
        i_r, i_s, trial_idx = _flat_to_indices(flat, sigma_targets.size, trials_per_point)

        progress = (flat - next_flat) + current["trial_fraction"]
        total = float(total_trials - next_flat)

        if now - last_log >= log_interval:
            eta = _format_eta(elapsed, max(progress, 1e-9), total)
            phase_name = {0: "relax", 1: "cal", 2: "drive"}.get(current["phase"], "?")
            logger.info(
                "Progress ~%.4f/%d (R_idx=%d, sigma_idx=%d, trial_idx=%d, phase=%s step=%d). Elapsed %.1fs. %s",
                progress,
                total_trials - next_flat,
                i_r,
                i_s,
                trial_idx,
                phase_name,
                current["step"],
                elapsed,
                eta,
            )
            last_log = now

        if now - last_checkpoint >= checkpoint_interval:
            checkpoint_now(flat, state)
            logger.info("Checkpoint saved at %s (flat=%d phase=%d step=%d)", checkpoint_path, flat, current["phase"], current["step"])
            last_checkpoint = now

        if elapsed >= args.max_wall_seconds:
            checkpoint_now(flat, state)
            logger.error("RUN_TIMEOUT_CHECKPOINT_WRITTEN")
            print("RUN_TIMEOUT_CHECKPOINT_WRITTEN")
            raise RunTimeout()

    def _flat_to_indices(flat: int, n_sigma: int, n_trials: int) -> tuple[int, int, int]:
        per_r = n_sigma * n_trials
        i_r = flat // per_r
        rem = flat % per_r
        i_s = rem // n_trials
        trial = rem % n_trials
        return i_r, i_s, trial

    base_seed = int(args.seed)

    # Main loop over trials.
    try:
        flat = next_flat
        while flat < total_trials:
            current["flat"] = flat
            i_r, i_s, trial_idx = _flat_to_indices(flat, sigma_targets.size, trials_per_point)
            radius = float(r_values[i_r])
            sigma_target = float(sigma_targets[i_s])
            shell_mask = shell_masks[i_r]
            disk_mask = disk_masks[i_r]
            disk_radius = radius / 2.0
            r_loop = radius + 1.0

            disk_min_seen = float("inf")

            # Restore or initialize trial state.
            if active_state is not None and active_state["active_flat"] == flat:
                phase = int(active_state["active_phase"])
                phase_step = int(active_state["active_step"])
                w = active_state["w"].astype(np.float64, copy=False)
                phi = active_state["phi"].astype(np.float64, copy=False)
                eta_values = active_state.get("eta_values", None)
                if isinstance(eta_values, np.ndarray) and eta_values.size == 0:
                    eta_values = None
                n_series = active_state.get("drive_n_series", None)
                min_abs_series = active_state.get("drive_min_abs_series", None)
                if isinstance(n_series, np.ndarray) and n_series.size == 0:
                    n_series = None
                if isinstance(min_abs_series, np.ndarray) and min_abs_series.size == 0:
                    min_abs_series = None
                if n_series is None:
                    n_series_list: list[int] = []
                else:
                    n_series_list = [int(v) for v in n_series.tolist()]
                if min_abs_series is None:
                    min_abs_series_list: list[float] = []
                else:
                    min_abs_series_list = [float(v) for v in min_abs_series.tolist()]
                disk_min_seen = float(active_state.get("drive_min_abs_min", np.inf))
            else:
                phase = 0
                phase_step = 0
                w = 1.0 - 0.6 * np.exp(-(r_grid**2) / (2.0 * radius * radius))
                phi = np.zeros_like(w)
                eta_values = None
                n_series_list = []
                min_abs_series_list = []

            # Phase 0: relax
            if phase == 0:
                _apply_boundary_w(w)
                _apply_boundary_phi(phi, boundary, zero_boundary_phi)
                for step in range(phase_step, steps_relax):
                    current["phase"] = 0
                    current["step"] = step
                    current["trial_fraction"] = compute_trial_fraction(0, step)

                    _coupled_step(
                        w,
                        phi,
                        boundary,
                        zero_boundary_phi,
                        dt,
                        inv_dx2,
                        dx,
                        alpha,
                        beta,
                        gamma,
                        kappa,
                        eps_w,
                        workspace,
                    )

                    if step % 200 == 0:
                        state = {
                            "active_flat": flat,
                            "active_phase": 0,
                            "active_step": step + 1,
                            "w": w,
                            "phi": phi,
                        }
                        maybe_log_and_checkpoint(state)

                phase = 1
                phase_step = 0

            # Phase 1: calibration (always restartable from scratch).
            if phase == 1:
                current["phase"] = 1
                current["step"] = 0
                current["trial_fraction"] = compute_trial_fraction(1, 0)
                state = {
                    "active_flat": flat,
                    "active_phase": 1,
                    "active_step": 0,
                    "w": w,
                    "phi": phi,
                }
                maybe_log_and_checkpoint(state)

                seed = stable_seed(base_seed, i_r, i_s, trial_idx)
                rng = np.random.default_rng(seed)

                m_max = 64
                q = 1.0
                m_arr = np.arange(1, m_max + 1, dtype=float)
                sigma_m = m_arr ** (-q)
                a_m = rng.normal(0.0, sigma_m)
                b_m = rng.normal(0.0, sigma_m)
                eta_unscaled = _eta_from_coeffs(a_m, b_m, s_values, l_b)

                # Phase-only diffusion from phi=0 with Dirichlet boundary = eta.
                phi_cal = np.zeros_like(phi)
                _apply_boundary_phi(phi_cal, boundary, eta_unscaled)
                for step in range(steps_cal):
                    current["phase"] = 1
                    current["step"] = step
                    current["trial_fraction"] = compute_trial_fraction(1, step)

                    _phase_only_step(phi_cal, boundary, eta_unscaled, dt, inv_dx2, kappa, workspace)

                    if step % 200 == 0:
                        state = {
                            "active_flat": flat,
                            "active_phase": 1,
                            "active_step": 0,
                            "w": w,
                            "phi": phi,
                        }
                        maybe_log_and_checkpoint(state)

                sigma_meas = _compute_sigma_phi(phi_cal, shell_mask, dx, workspace)
                sigma_phi_measured_cal[i_r, i_s, trial_idx] = sigma_meas

                scale = sigma_target / (sigma_meas + 1e-12)
                a_m *= scale
                b_m *= scale
                eta_values = _eta_from_coeffs(a_m, b_m, s_values, l_b)

                # Sigma at start of drive
                phi_drive = phi
                _apply_boundary_phi(phi_drive, boundary, eta_values)
                sigma_phi_measured_start[i_r, i_s, trial_idx] = _compute_sigma_phi(
                    phi_drive, shell_mask, dx, workspace
                )

                phase = 2
                phase_step = 0

            # Phase 2: drive
            if phase == 2:
                assert eta_values is not None
                _apply_boundary_w(w)
                _apply_boundary_phi(phi, boundary, eta_values)

                # If we are resuming mid-drive, lists already contain samples.
                start_step = phase_step
                if disk_min_seen == float("inf"):
                    disk_min_seen = float(np.min(w[disk_mask]))
                for step in range(start_step, steps_drive + 1):
                    current["phase"] = 2
                    current["step"] = step
                    current["trial_fraction"] = compute_trial_fraction(2, step)

                    # Sample first (step=0) and then after update at subsequent steps.
                    if step % sample_stride == 0:
                        samples = sample_complex_from_polar_on_circle(
                            w,
                            phi,
                            0.0,
                            0.0,
                            r_loop,
                            m_loop,
                            x_min,
                            y_min,
                            dx,
                        )
                        n_series_list.append(winding_number_from_samples(samples))
                        min_abs_series_list.append(float(np.min(w[disk_mask])))

                    if step == steps_drive:
                        break

                    _coupled_step(
                        w,
                        phi,
                        boundary,
                        eta_values,
                        dt,
                        inv_dx2,
                        dx,
                        alpha,
                        beta,
                        gamma,
                        kappa,
                        eps_w,
                        workspace,
                    )

                    # Track absolute minimum over time inside the disk.
                    disk_min_seen = min(disk_min_seen, float(np.min(w[disk_mask])))

                    if step % 200 == 0:
                        state = {
                            "active_flat": flat,
                            "active_phase": 2,
                            "active_step": step + 1,
                            "w": w,
                            "phi": phi,
                            "eta_values": eta_values,
                            "drive_n_series": np.array(n_series_list, dtype=np.int64),
                            "drive_min_abs_series": np.array(min_abs_series_list, dtype=np.float64),
                            "drive_min_abs_min": float(disk_min_seen),
                        }
                        maybe_log_and_checkpoint(state)

                # Record trial outputs.
                if n_series_list:
                    n_final[i_r, i_s, trial_idx] = float(n_series_list[-1])
                min_abs_min[i_r, i_s, trial_idx] = float(disk_min_seen)

                # Persistence check.
                n_array = np.array(n_series_list, dtype=int)
                core_array = np.array(min_abs_series_list, dtype=float) < epsilon_core
                times = np.arange(0.0, dt_sample * len(n_array), dt_sample)

                success = False
                start_idx = None
                for idx, n_val in enumerate(n_array):
                    if n_val != 0:
                        if start_idx is None:
                            start_idx = idx
                    else:
                        if start_idx is not None:
                            if times[idx - 1] - times[start_idx] >= t_persist:
                                if core_array[start_idx:idx].any():
                                    success = True
                                    break
                            start_idx = None
                if not success and start_idx is not None:
                    if times[-1] - times[start_idx] >= t_persist:
                        if core_array[start_idx:].any():
                            success = True

                if success:
                    success_counts[i_r, i_s] += 1

                # Save representative debug artifacts.
                for label, sigma_idx in repr_sigma_indices.items():
                    if i_r == repr_r_index and i_s == sigma_idx and trial_idx == 0:
                        debug_dir = output_dir / f"debug_{label}"
                        ensure_dir(debug_dir)
                        np.save(debug_dir / "n_time_series.npy", n_array)
                        np.save(debug_dir / "min_absPsi_time_series.npy", np.array(min_abs_series_list, dtype=float))
                        np.save(debug_dir / "snapshot_w.npy", w)
                        np.save(debug_dir / "snapshot_phi.npy", phi)

                # Trial complete.
                flat += 1
                active_state = None

                # Update checkpoint to reflect completed trial.
                checkpoint_now(flat, None)
                continue

            raise RuntimeError(f"Unexpected phase state: {phase}")

    except RunTimeout:
        sys.exit(1)

    # Confidence intervals.
    ci_low = np.zeros_like(success_counts, dtype=float)
    ci_high = np.zeros_like(success_counts, dtype=float)
    p_hat = success_counts / trial_counts
    for i in range(success_counts.shape[0]):
        for j in range(success_counts.shape[1]):
            low, high = _wilson_ci(int(success_counts[i, j]), int(trial_counts[i, j]))
            ci_low[i, j] = low
            ci_high[i, j] = high

    # Stage 1 JSON (match Section 6.1 schema).
    stage1_json = {
        "stage": 1,
        "map_resolution": [int(r_values.size), int(sigma_targets.size)],
        "param_R": [float(r_values.min()), float(r_values.max())],
        "param_sigma": [float(sigma_targets.min()), float(sigma_targets.max())],
        "transition_matrix": p_hat.tolist(),
        "transition_matrix_ci_low": ci_low.tolist(),
        "transition_matrix_ci_high": ci_high.tolist(),
        "solver": {
            "dimension": dim,
            "L": l,
            "N": n,
            "dt": dt,
            "t_relax": t_relax,
            "t_drive": t_drive,
            "forcing": {
                "k_vec": [0.0, 0.0],
                "boundary_noise_powerlaw_q": 1.0,
                "boundary_noise_mmax": 64,
                "boundary_noise_calibration": (
                    "two-pass: phase-only with w=1 for T_cal=5, then rescale coefficients to hit sigma_phi_target"
                ),
            },
        },
        "trials_per_point": trials_per_point,
    }

    write_json(output_dir / "stage1_selection.json", stage1_json)

    np.savez(
        output_dir / "raw_aggregates.npz",
        R_values=r_values,
        sigma_phi_values=sigma_targets,
        success_counts=success_counts,
        trial_counts=trial_counts,
        sigma_phi_measured_calibration=sigma_phi_measured_cal,
        sigma_phi_measured_start_drive=sigma_phi_measured_start,
        n_final=n_final,
        min_absPsi_min_over_time=min_abs_min,
    )

    write_run_manifest(output_dir, Path.cwd(), sys.argv)

    # Acceptance tests.
    fail_note = None
    if np.any(p_hat[:, 0] >= 0.05):
        fail_note = "STAGE1_FAIL_ALWAYS_NUCLEATES"
    elif not np.any((sigma_targets >= 2.0)[None, :] & (p_hat > 0.30)):
        fail_note = "STAGE1_FAIL_NO_NUCLEATION"

    if fail_note is not None:
        (output_dir / "stage1_failure_note.txt").write_text(fail_note)


if __name__ == "__main__":
    main()
