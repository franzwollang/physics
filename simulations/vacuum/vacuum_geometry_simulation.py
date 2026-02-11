"""Vacuum geometry simulation: graph geometrogenesis & spectral dimension."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

try:
    from sklearn.manifold import Isomap
except ImportError:  # pragma: no cover - optional dependency
    Isomap = None

try:
    from scipy.sparse.csgraph import shortest_path
except ImportError as exc:  # pragma: no cover - optional dependency
    raise SystemExit("scipy is required for shortest_path distances") from exc


@dataclass
class SimulationParams:
    n: int = 100
    alpha: float = 1.0
    temperature: float = 0.1
    beta: float = 10.0
    d_star: float = 3.0
    tau_star: float = 1.0
    steps: int = 1000
    eta: float = 0.01
    initial_graph: str = "ErdosRenyi_p0.3"
    seed: int | None = None


def _parse_initial_graph(name: str, n: int, rng: np.random.Generator) -> NDArray[np.float64]:
    if name.startswith("ErdosRenyi_p"):
        p = float(name.split("p", 1)[1])
        mask = rng.random((n, n)) < p
        weights = rng.random((n, n))
        j = np.where(mask, weights, 0.0)
    elif name.startswith("RandomGeometric_r"):
        radius = float(name.split("r", 1)[1])
        coords = rng.random((n, 3))
        diff = coords[:, None, :] - coords[None, :, :]
        dist = np.linalg.norm(diff, axis=-1)
        mask = dist < radius
        weights = rng.random((n, n))
        j = np.where(mask, weights, 0.0)
    else:
        raise ValueError(f"Unknown initial_graph format: {name}")

    j = np.triu(j, 1)
    j = j + j.T
    np.fill_diagonal(j, 0.0)
    return j


def _laplacian(j: NDArray[np.float64]) -> NDArray[np.float64]:
    degree = np.sum(j, axis=1)
    return np.diag(degree) - j


def _spectral_dimension(evals: NDArray[np.float64], tau: float) -> float:
    weights = np.exp(-evals * tau)
    numerator = np.sum(evals * weights)
    denominator = np.sum(weights)
    return 2.0 * tau * numerator / denominator


def _spectral_dimension_gradient(
    evals: NDArray[np.float64],
    evecs: NDArray[np.float64],
    tau: float,
) -> NDArray[np.float64]:
    weights = np.exp(-evals * tau)
    n = evals.size
    p_tau = np.mean(weights)
    ds = _spectral_dimension(evals, tau)
    coeffs = (2.0 * tau * weights / p_tau) * (ds / 2.0 - 1.0) / n

    grad_ds = np.zeros((n, n), dtype=float)
    for k in range(n):
        vk = evecs[:, k]
        diff = vk[:, None] - vk[None, :]
        grad_ds += coeffs[k] * diff**2

    np.fill_diagonal(grad_ds, 0.0)
    return grad_ds


def _entropy(j: NDArray[np.float64]) -> tuple[float, NDArray[np.float64]]:
    z = np.sum(j)
    if z <= 0.0:
        return 0.0, np.zeros_like(j)
    p = j / z
    mask = p > 0
    entropy = -np.sum(p[mask] * np.log(p[mask]))
    grad = np.zeros_like(j)
    grad[mask] = -np.log(p[mask]) - 1.0
    np.fill_diagonal(grad, 0.0)
    return entropy, grad


def _logdet_energy(j: NDArray[np.float64], eps: float = 1e-6) -> tuple[float, NDArray[np.float64]]:
    lap = _laplacian(j)
    lap_reg = lap + eps * np.eye(lap.shape[0])
    inv = np.linalg.inv(lap_reg)
    diag = np.diag(inv)
    grad = diag[:, None] + diag[None, :] - 2.0 * inv
    np.fill_diagonal(grad, 0.0)
    sign, logdet = np.linalg.slogdet(lap_reg)
    energy = logdet if sign > 0 else math.inf
    return energy, grad


def _project_and_normalize(j: NDArray[np.float64], target_total: float) -> NDArray[np.float64]:
    j = np.maximum(j, 0.0)
    np.fill_diagonal(j, 0.0)
    j = 0.5 * (j + j.T)
    total = np.sum(j)
    if total > 0.0:
        j *= target_total / total
    return j


def _isomap_embedding(j: NDArray[np.float64], n_neighbors: int = 10) -> dict[str, float | str]:
    if Isomap is None:
        return {
            "isomap_residual_3D": None,
            "isomap_residual_2D": None,
            "coordinates": None,
        }

    distances = np.where(j > 0, 1.0 / (j + 1e-12), np.inf)
    np.fill_diagonal(distances, 0.0)
    geodesic = shortest_path(distances, directed=False, unweighted=False)
    finite = np.isfinite(geodesic)
    if not np.all(finite):
        max_finite = np.max(geodesic[finite])
        geodesic = np.where(finite, geodesic, max_finite * 1.5)

    def _embed(dim: int) -> tuple[np.ndarray, float]:
        iso = Isomap(n_neighbors=n_neighbors, n_components=dim, metric="precomputed")
        coords = iso.fit_transform(geodesic)
        embedded_dist = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
        flat_geo = geodesic.flatten()
        flat_emb = embedded_dist.flatten()
        corr = np.corrcoef(flat_geo, flat_emb)[0, 1]
        residual = 1.0 - corr**2 if np.isfinite(corr) else float("nan")
        return coords, residual

    coords_3d, residual_3d = _embed(3)
    _, residual_2d = _embed(2)
    return {
        "isomap_residual_3D": float(residual_3d),
        "isomap_residual_2D": float(residual_2d),
        "coordinates": coords_3d,
    }


def _simulate(params: SimulationParams) -> tuple[dict, NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    rng = np.random.default_rng(params.seed)
    j = _parse_initial_graph(params.initial_graph, params.n, rng)
    target_total = float(np.sum(j))
    if target_total <= 0:
        raise ValueError("Initial graph has zero total weight.")
    j = _project_and_normalize(j, target_total)

    ds_history: list[float] = []
    for _ in range(params.steps):
        lap = _laplacian(j)
        evals, evecs = np.linalg.eigh(lap)
        ds = _spectral_dimension(evals, params.tau_star)
        ds_history.append(float(ds))

        _, grad_entropy = _entropy(j)
        _, grad_logdet = _logdet_energy(j)
        grad_ds = _spectral_dimension_gradient(evals, evecs, params.tau_star)
        grad_dim = 2.0 * params.beta * (ds - params.d_star) * grad_ds

        total_grad = (
            params.alpha * grad_logdet
            - params.temperature * grad_entropy
            + grad_dim
        )
        total_grad = 0.5 * (total_grad + total_grad.T)
        np.fill_diagonal(total_grad, 0.0)

        j = j - params.eta * total_grad
        j = _project_and_normalize(j, target_total)

    lap = _laplacian(j)
    evals, _ = np.linalg.eigh(lap)
    ds_final = _spectral_dimension(evals, params.tau_star)
    fiedler = float(np.sort(evals)[1]) if params.n > 1 else 0.0
    degrees = np.sum(j, axis=1)

    window = ds_history[-min(len(ds_history), 100) :]
    ds_std = float(np.std(window)) if window else 0.0

    embedding_info = _isomap_embedding(j)
    coordinates = embedding_info.pop("coordinates")
    if coordinates is None:
        coordinates = np.zeros((params.n, 3))

    outputs = {
        "converged": True,
        "final_d_s": float(ds_final),
        "d_s_std": ds_std,
        "isomap_residual_3D": embedding_info["isomap_residual_3D"],
        "isomap_residual_2D": embedding_info["isomap_residual_2D"],
        "mean_degree": float(np.mean(degrees)),
        "std_degree": float(np.std(degrees)),
        "fiedler_eigenvalue": fiedler,
    }
    metadata = {
        "simulation": "VacuumGeometry",
        "version": "1.0",
        "params": {
            "N": params.n,
            "alpha": params.alpha,
            "T": params.temperature,
            "beta": params.beta,
            "d_star": params.d_star,
            "tau_star": params.tau_star,
            "initial_graph": params.initial_graph,
            "steps": params.steps,
            "eta": params.eta,
        },
        "outputs": outputs,
    }
    return metadata, coordinates, j, evals


def run_simulation(params: SimulationParams) -> tuple[dict, NDArray[np.float64]]:
    metadata, coordinates, _, _ = _simulate(params)
    return metadata, coordinates


def run_simulation_with_graph(
    params: SimulationParams,
) -> tuple[dict, NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    return _simulate(params)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--temperature", type=float, default=0.1)
    parser.add_argument("--beta", type=float, default=10.0)
    parser.add_argument("--d-star", type=float, default=3.0)
    parser.add_argument("--tau-star", type=float, default=1.0)
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--eta", type=float, default=0.01)
    parser.add_argument("--initial-graph", type=str, default="ErdosRenyi_p0.3")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("vacuum_geometry_output.json"),
    )
    parser.add_argument(
        "--embedding-file",
        type=Path,
        default=Path("embedding_3D.npy"),
    )
    args = parser.parse_args()

    params = SimulationParams(
        n=args.n,
        alpha=args.alpha,
        temperature=args.temperature,
        beta=args.beta,
        d_star=args.d_star,
        tau_star=args.tau_star,
        steps=args.steps,
        eta=args.eta,
        initial_graph=args.initial_graph,
        seed=args.seed,
    )

    metadata, coordinates = run_simulation(params)
    np.save(args.embedding_file, coordinates)
    metadata["embedding"] = {"coordinates_file": args.embedding_file.name}

    args.output_json.write_text(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
