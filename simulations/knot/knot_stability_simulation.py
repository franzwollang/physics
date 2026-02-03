"""Knot stability simulation for phase-dark soliton effective strings."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


@dataclass
class CanonicalParams:
    alpha_grad: float = 1.0
    beta_pot: float = 0.5
    gamma: float = 1.0
    kappa: float = 1.0

    @property
    def w_star(self) -> float:
        return math.sqrt(self.beta_pot / (2 * self.gamma))

    @property
    def ell(self) -> float:
        return math.sqrt(self.alpha_grad / (2 * self.beta_pot))

    @property
    def tension(self) -> float:
        return math.pi * self.alpha_grad * self.w_star**2 / 2

    @property
    def bending(self) -> float:
        return self.tension * self.ell**2

    @property
    def barrier_estimate(self) -> float:
        return (self.ell**3) * self.beta_pot * self.w_star**2


@dataclass
class SimulationParams:
    knot_type: str = "Trefoil"
    n_beads: int = 100
    k_T: float = 100.0
    k_B: float = 1.0
    sigma: float = 1.0
    epsilon_LJ: float = 10.0
    l0: float = 1.0


def trefoil_coordinates(n_beads: int) -> np.ndarray:
    t = np.linspace(0.0, 2 * np.pi, n_beads, endpoint=False)
    x = np.sin(t) + 2 * np.sin(2 * t)
    y = np.cos(t) - 2 * np.cos(2 * t)
    z = -np.sin(3 * t)
    return np.column_stack([x, y, z])


def rescale_to_length(coords: np.ndarray, target_length: float) -> np.ndarray:
    segment_lengths = np.linalg.norm(np.roll(coords, -1, axis=0) - coords, axis=1)
    current_length = segment_lengths.sum()
    scale = target_length / current_length
    return coords * scale


def energy_function(flat_coords: np.ndarray, params: SimulationParams) -> float:
    coords = flat_coords.reshape((-1, 3))
    n_beads = coords.shape[0]

    diffs = np.roll(coords, -1, axis=0) - coords
    seg_lengths = np.linalg.norm(diffs, axis=1)
    seg_lengths = np.maximum(seg_lengths, 1e-8)
    tension_energy = 0.5 * params.k_T * np.sum((seg_lengths - params.l0) ** 2)

    tangents = diffs / seg_lengths[:, None]
    prev_tangents = np.roll(tangents, 1, axis=0)
    cos_angles = np.einsum("ij,ij->i", tangents, prev_tangents)
    cos_angles = np.clip(cos_angles, -1.0, 1.0)
    bending_energy = params.k_B * np.sum(1.0 - cos_angles)

    cutoff = 2.5 * params.sigma
    diffs_matrix = coords[:, None, :] - coords[None, :, :]
    dist2 = np.sum(diffs_matrix**2, axis=2)
    dist = np.sqrt(dist2)
    indices = np.arange(n_beads)
    index_diff = np.abs(indices[:, None] - indices[None, :])
    separation = np.minimum(index_diff, n_beads - index_diff)
    upper = np.triu(np.ones((n_beads, n_beads), dtype=bool), k=1)
    lj_mask = (separation > 2) & (dist < cutoff) & upper
    lj_energy = np.sum(
        4
        * params.epsilon_LJ
        * (params.sigma**12)
        / np.maximum(dist2, 1e-12) ** 6
        * lj_mask
    )

    return tension_energy + bending_energy + lj_energy


def radius_of_gyration(coords: np.ndarray) -> float:
    center = coords.mean(axis=0)
    return math.sqrt(np.mean(np.sum((coords - center) ** 2, axis=1)))


def min_nonadjacent_distance(coords: np.ndarray) -> float:
    n_beads = coords.shape[0]
    diffs_matrix = coords[:, None, :] - coords[None, :, :]
    dist2 = np.sum(diffs_matrix**2, axis=2)
    dist = np.sqrt(dist2)
    indices = np.arange(n_beads)
    index_diff = np.abs(indices[:, None] - indices[None, :])
    separation = np.minimum(index_diff, n_beads - index_diff)
    mask = (separation > 1) & np.triu(np.ones((n_beads, n_beads), dtype=bool), k=1)
    if not np.any(mask):
        return float("inf")
    return float(np.min(dist[mask]))


def approx_writhe(coords: np.ndarray) -> float | None:
    try:
        from pyknotid.spacecurves import Knot
    except ImportError:
        return None

    knot = Knot(coords)
    return float(knot.writhe())


def run_simulation(
    params: SimulationParams,
    canonical: CanonicalParams,
    target_length: float = 50.0,
    maxiter: int = 300,
    maxfun: int | None = None,
    method: str = "L-BFGS-B",
    ftol: float = 1e-9,
    gtol: float = 1e-6,
    maxls: int = 40,
) -> dict:
    coords = trefoil_coordinates(params.n_beads)
    coords = rescale_to_length(coords, target_length)

    initial_energy = energy_function(coords.flatten(), params)

    options = {"maxiter": maxiter}
    if maxfun is not None:
        options["maxfun"] = maxfun
    if method.upper() == "L-BFGS-B":
        options["ftol"] = ftol
        options["gtol"] = gtol
        options["maxls"] = maxls
    elif method.upper() == "POWELL":
        options["ftol"] = ftol

    result = minimize(
        energy_function,
        coords.flatten(),
        args=(params,),
        method=method,
        options=options,
    )

    final_coords = result.x.reshape((-1, 3))
    final_energy = float(result.fun)
    rg = radius_of_gyration(final_coords)
    total_length = np.linalg.norm(np.roll(final_coords, -1, axis=0) - final_coords, axis=1).sum()
    r_star = math.sqrt(params.k_B / 2.0)
    barrier = canonical.barrier_estimate
    barrier_ratio = barrier / final_energy if final_energy > 0 else float("inf")
    frozen = barrier > 1.0
    min_distance = min_nonadjacent_distance(final_coords)
    topology_preserved = min_distance > 0.8 * params.sigma
    writhe_value = approx_writhe(final_coords)
    sigma_over_m = math.pi * rg**2 / (canonical.tension * total_length)

    return {
        "simulation": "KnotStability",
        "version": "1.0",
        "params": asdict(params),
        "canonical_params": {
            "alpha_grad": canonical.alpha_grad,
            "beta_pot": canonical.beta_pot,
            "gamma": canonical.gamma,
            "ell": canonical.ell,
            "T_string": canonical.tension,
            "B_string": canonical.bending,
        },
        "outputs": {
            "converged": bool(result.success),
            "optimizer_status": int(result.status),
            "optimizer_message": result.message,
            "n_iterations": int(result.nit),
            "n_function_evals": int(result.nfev),
            "initial_energy": float(initial_energy),
            "final_energy": final_energy,
            "radius_of_gyration": float(rg),
            "metastable_radius_R_star": float(r_star),
            "topology_preserved": topology_preserved,
            "final_writhe": writhe_value,
            "sigma_over_m": float(sigma_over_m),
        },
        "stability_analysis": {
            "barrier_estimate": float(barrier),
            "barrier_to_energy_ratio": float(barrier_ratio),
            "frozen": frozen,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("simulations/results"))
    parser.add_argument("--kB", type=float, nargs="+", default=[0.1, 1.0, 10.0])
    parser.add_argument("--n-beads", type=int, default=100)
    parser.add_argument("--kT", type=float, default=100.0)
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--epsilon", type=float, default=10.0)
    parser.add_argument("--l0", type=float, default=1.0)
    parser.add_argument("--maxiter", type=int, default=300)
    parser.add_argument("--maxfun", type=int, default=None)
    parser.add_argument("--method", type=str, default="L-BFGS-B")
    parser.add_argument("--ftol", type=float, default=1e-9)
    parser.add_argument("--gtol", type=float, default=1e-6)
    parser.add_argument("--maxls", type=int, default=40)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    canonical = CanonicalParams()
    target_length = 50.0 * canonical.ell

    for k_B in args.kB:
        params = SimulationParams(
            n_beads=args.n_beads,
            k_T=args.kT,
            k_B=k_B,
            sigma=args.sigma,
            epsilon_LJ=args.epsilon,
            l0=args.l0,
        )
        payload = run_simulation(
            params,
            canonical,
            target_length=target_length,
            maxiter=args.maxiter,
            maxfun=args.maxfun,
            method=args.method,
            ftol=args.ftol,
            gtol=args.gtol,
            maxls=args.maxls,
        )
        output_path = args.output_dir / f"knot_stability_kB_{k_B:.2f}.json"
        output_path.write_text(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
