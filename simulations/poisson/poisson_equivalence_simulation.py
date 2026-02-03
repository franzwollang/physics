#!/usr/bin/env python3
"""Simulate Poisson-equivalent noise variance from stiffness perturbations.

Implements the plan in simulations/simulation_plan_poisson_equivalence.md.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np
from scipy.optimize import curve_fit
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator, cg, spilu


@dataclass
class SimulationParams:
    N: int
    epsilon: float
    R_source: int
    K_samples: int
    boundary: str = "Dirichlet"


@dataclass
class SimulationOutputs:
    fit_A: float
    fit_B: float
    fit_R_squared: float
    C_eff: float
    linearity_check: bool


@dataclass
class Diagnostics:
    cg_iterations_avg: float
    cg_nonconverged: int
    hutchinson_variance: float
    exact_diag_rmse: Optional[float] = None


def build_kappa_field(N: int, epsilon: float, R_source: int) -> np.ndarray:
    grid = np.indices((N, N, N), dtype=float)
    center = (N - 1) / 2.0
    r = np.sqrt(
        (grid[0] - center) ** 2 + (grid[1] - center) ** 2 + (grid[2] - center) ** 2
    )
    kappa = np.ones((N, N, N), dtype=float)
    kappa[r < R_source] += epsilon
    return kappa


def harmonic_mean(a: float, b: float) -> float:
    return 2.0 * a * b / (a + b)


def build_laplacian(kappa: np.ndarray) -> csr_matrix:
    N = kappa.shape[0]
    n = N**3
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []

    def index(i: int, j: int, k: int) -> int:
        return (i * N + j) * N + k

    for i in range(N):
        for j in range(N):
            for k in range(N):
                idx = index(i, j, k)
                if i in (0, N - 1) or j in (0, N - 1) or k in (0, N - 1):
                    rows.append(idx)
                    cols.append(idx)
                    data.append(1.0)
                    continue
                diag = 0.0
                for di, dj, dk in (
                    (1, 0, 0),
                    (-1, 0, 0),
                    (0, 1, 0),
                    (0, -1, 0),
                    (0, 0, 1),
                    (0, 0, -1),
                ):
                    ni, nj, nk = i + di, j + dj, k + dk
                    edge = harmonic_mean(kappa[i, j, k], kappa[ni, nj, nk])
                    diag += edge
                    nidx = index(ni, nj, nk)
                    rows.append(idx)
                    cols.append(nidx)
                    data.append(-edge)
                rows.append(idx)
                cols.append(idx)
                data.append(diag)
    return csr_matrix((data, (rows, cols)), shape=(n, n))


def build_preconditioner(L: csr_matrix) -> Optional[LinearOperator]:
    try:
        ilu = spilu(L.tocsc(), drop_tol=1e-4, fill_factor=10)
    except Exception:
        return None

    def solve(x: np.ndarray) -> np.ndarray:
        return ilu.solve(x)

    return LinearOperator(L.shape, matvec=solve)


def solve_cg(
    L: csr_matrix,
    v: np.ndarray,
    *,
    preconditioner: Optional[LinearOperator],
    atol: float = 1e-8,
    rtol: float = 1e-5,
    maxiter: int = 2000,
) -> tuple[np.ndarray, int, int]:
    iters = 0

    def callback(_: np.ndarray) -> None:
        nonlocal iters
        iters += 1

    u, info = cg(L, v, M=preconditioner, atol=atol, rtol=rtol, maxiter=maxiter, callback=callback)
    if info != 0 and preconditioner is not None:
        retry_iters = 0

        def callback_retry(_: np.ndarray) -> None:
            nonlocal retry_iters
            retry_iters += 1

        u, info = cg(
            L, v, M=None, atol=atol, rtol=rtol, maxiter=maxiter * 2, callback=callback_retry
        )
        iters += retry_iters
    return u, info, iters


def hutchinson_diag(
    L: csr_matrix,
    K: int,
    rng: np.random.Generator,
    boundary_mask: np.ndarray,
    preconditioner: Optional[LinearOperator] = None,
) -> tuple[np.ndarray, float, int, float]:
    n = L.shape[0]
    accumulator = np.zeros(n)
    iter_counts: list[int] = []
    sample_means: list[float] = []
    nonconverged = 0

    for _ in range(K):
        v = rng.choice([-1.0, 1.0], size=n)
        v[boundary_mask] = 0.0
        u, info, iters = solve_cg(L, v, preconditioner=preconditioner)
        if info != 0:
            nonconverged += 1
        accumulator += v * u
        iter_counts.append(iters)
        sample_means.append(float(np.mean(v * u)))

    diag_estimate = accumulator / K
    variance = float(np.var(sample_means, ddof=1)) if K > 1 else 0.0
    return diag_estimate, float(np.mean(iter_counts)), nonconverged, variance


def hutchinson_delta_diag(
    L: csr_matrix,
    L0: csr_matrix,
    K: int,
    rng: np.random.Generator,
    boundary_mask: np.ndarray,
    preconditioner: Optional[LinearOperator] = None,
    preconditioner0: Optional[LinearOperator] = None,
) -> tuple[np.ndarray, float, int, float]:
    n = L.shape[0]
    accumulator = np.zeros(n)
    iter_counts: list[int] = []
    sample_means: list[float] = []
    nonconverged = 0

    for _ in range(K):
        v = rng.choice([-1.0, 1.0], size=n)
        v[boundary_mask] = 0.0
        u, info, iters = solve_cg(L, v, preconditioner=preconditioner)
        if info != 0:
            nonconverged += 1
        u0, info0, iters0 = solve_cg(L0, v, preconditioner=preconditioner0)
        if info0 != 0:
            nonconverged += 1
        iters += iters0
        delta = v * (u - u0)
        accumulator += delta
        iter_counts.append(iters)
        sample_means.append(float(np.mean(delta)))

    diag_estimate = accumulator / K
    variance = float(np.var(sample_means, ddof=1)) if K > 1 else 0.0
    return diag_estimate, float(np.mean(iter_counts)), nonconverged, variance


def radial_profile(values: np.ndarray, center: float) -> tuple[np.ndarray, np.ndarray]:
    N = values.shape[0]
    grid = np.indices((N, N, N), dtype=float)
    r = np.sqrt(
        (grid[0] - center) ** 2 + (grid[1] - center) ** 2 + (grid[2] - center) ** 2
    )
    r_flat = r.ravel()
    values_flat = values.ravel()
    bins = np.floor(r_flat).astype(int)
    max_bin = bins.max()
    sums = np.bincount(bins, weights=values_flat, minlength=max_bin + 1)
    counts = np.bincount(bins, minlength=max_bin + 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        profile = np.where(counts > 0, sums / counts, 0.0)
    radii = np.arange(len(profile)) + 0.5
    return radii, profile


def fit_profile(radii: np.ndarray, profile: np.ndarray, R_source: int, N: int) -> tuple[float, float, float]:
    mask = (radii > R_source) & (radii < N / 4.0) & np.isfinite(profile)
    fit_r = radii[mask]
    fit_y = profile[mask]

    if fit_r.size < 2:
        return float("nan"), float("nan"), 0.0

    def model(r: np.ndarray, A: float, B: float) -> np.ndarray:
        return A / r + B

    try:
        popt, _ = curve_fit(model, fit_r, fit_y, p0=(1.0, 0.0), maxfev=10000)
    except Exception:
        return float("nan"), float("nan"), 0.0
    pred = model(fit_r, *popt)
    ss_res = float(np.sum((fit_y - pred) ** 2))
    ss_tot = float(np.sum((fit_y - np.mean(fit_y)) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return float(popt[0]), float(popt[1]), r_squared


def compute_exact_diag(L: csr_matrix) -> np.ndarray:
    dense = L.toarray()
    inv = np.linalg.inv(dense)
    return np.diag(inv)


def run_simulation(
    params: SimulationParams,
    *,
    seed: int,
    linearity_epsilons: Optional[Iterable[float]] = None,
    validate_exact: bool = False,
) -> tuple[SimulationOutputs, Diagnostics]:
    rng = np.random.default_rng(seed)
    kappa = build_kappa_field(params.N, params.epsilon, params.R_source)
    L = build_laplacian(kappa)
    kappa0 = np.ones_like(kappa)
    L0 = build_laplacian(kappa0)

    boundary_mask = np.zeros(params.N**3, dtype=bool)
    N = params.N
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if i in (0, N - 1) or j in (0, N - 1) or k in (0, N - 1):
                    boundary_mask[(i * N + j) * N + k] = True

    preconditioner = build_preconditioner(L)
    preconditioner0 = build_preconditioner(L0)

    delta_diag, cg_iters, nonconverged, variance = hutchinson_delta_diag(
        L,
        L0,
        params.K_samples,
        rng,
        boundary_mask,
        preconditioner=preconditioner,
        preconditioner0=preconditioner0,
    )

    delta_field = delta_diag.reshape((params.N, params.N, params.N))
    center = (params.N - 1) / 2.0
    radii, profile = radial_profile(delta_field, center)
    fit_A, fit_B, fit_r2 = fit_profile(radii, profile, params.R_source, params.N)

    grid = np.indices((params.N, params.N, params.N), dtype=float)
    r = np.sqrt(
        (grid[0] - center) ** 2 + (grid[1] - center) ** 2 + (grid[2] - center) ** 2
    )
    source_volume = float(np.sum(r < params.R_source))
    C_eff = 2.0 * fit_A / (params.epsilon * source_volume) if params.epsilon != 0 else 0.0

    linearity_check = False
    if linearity_epsilons:
        fit_As = []
        for eps in linearity_epsilons:
            local_params = SimulationParams(
                N=params.N,
                epsilon=eps,
                R_source=params.R_source,
                K_samples=params.K_samples,
                boundary=params.boundary,
            )
            local_outputs, _ = run_simulation(
                local_params,
                seed=seed,
                linearity_epsilons=None,
                validate_exact=False,
            )
            fit_As.append(local_outputs.fit_A)
        slope, intercept = np.polyfit(list(linearity_epsilons), fit_As, 1)
        pred = slope * np.array(list(linearity_epsilons)) + intercept
        ss_res = float(np.sum((np.array(fit_As) - pred) ** 2))
        ss_tot = float(np.sum((np.array(fit_As) - np.mean(fit_As)) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        linearity_check = r2 > 0.95

    exact_rmse = None
    if validate_exact and params.N <= 20:
        exact_diag = compute_exact_diag(L)
        est_diag, _, _, _ = hutchinson_diag(L, params.K_samples, rng, boundary_mask, preconditioner)
        exact_rmse = float(np.sqrt(np.mean((exact_diag - est_diag) ** 2)))

    outputs = SimulationOutputs(
        fit_A=fit_A,
        fit_B=fit_B,
        fit_R_squared=fit_r2,
        C_eff=C_eff,
        linearity_check=linearity_check,
    )
    diagnostics = Diagnostics(
        cg_iterations_avg=cg_iters,
        cg_nonconverged=nonconverged,
        hutchinson_variance=variance,
        exact_diag_rmse=exact_rmse,
    )
    return outputs, diagnostics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Poisson equivalence simulation")
    parser.add_argument("--N", type=int, default=40, help="Lattice size")
    parser.add_argument("--epsilon", type=float, default=0.5, help="Coupling strength")
    parser.add_argument("--R-source", type=int, default=3, dest="R_source", help="Source radius")
    parser.add_argument("--K-samples", type=int, default=100, dest="K_samples", help="Hutchinson samples")
    parser.add_argument("--output", type=str, default="poisson_equivalence_results.json")
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--linearity-eps", type=float, nargs="*", default=None)
    parser.add_argument("--validate-exact", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    params = SimulationParams(
        N=args.N,
        epsilon=args.epsilon,
        R_source=args.R_source,
        K_samples=args.K_samples,
    )
    outputs, diagnostics = run_simulation(
        params,
        seed=args.seed,
        linearity_epsilons=args.linearity_eps,
        validate_exact=args.validate_exact,
    )
    result = {
        "simulation": "PoissonEquivalence",
        "version": "1.0",
        "params": {
            "N": params.N,
            "epsilon": params.epsilon,
            "R_source": params.R_source,
            "K_samples": params.K_samples,
            "boundary": params.boundary,
        },
        "canonical_params": {
            "kappa_0": 1.0,
            "alpha_grad": 1.0,
            "beta_pot": 0.5,
        },
        "outputs": {
            "fit_A": outputs.fit_A,
            "fit_B": outputs.fit_B,
            "fit_R_squared": outputs.fit_R_squared,
            "C_eff": outputs.C_eff,
            "linearity_check": outputs.linearity_check,
        },
        "diagnostics": {
            "cg_iterations_avg": diagnostics.cg_iterations_avg,
            "cg_nonconverged": diagnostics.cg_nonconverged,
            "hutchinson_variance": diagnostics.hutchinson_variance,
        },
    }
    if diagnostics.exact_diag_rmse is not None:
        result["diagnostics"]["exact_diag_rmse"] = diagnostics.exact_diag_rmse

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
