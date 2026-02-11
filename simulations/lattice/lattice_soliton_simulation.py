"""Lattice formation simulation for phase-dark solitons.

Implements an overdamped Langevin quench for a complex scalar field on a
periodic cubic lattice and extracts vortex loop statistics via plaquette
winding on the dual lattice.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np


@dataclass
class SimulationParams:
    size: int
    hardness: float
    steps: int
    dt: float
    t_initial: float
    t_final: float
    seed: int
    output_path: Path
    progress_interval: float


def parse_args() -> SimulationParams:
    parser = argparse.ArgumentParser(
        description="Lattice formation and vortex loop statistics simulation.",
    )
    parser.add_argument("--size", type=int, default=64, help="Lattice size N.")
    parser.add_argument("--lambda", dest="hardness", type=float, default=1.0)
    parser.add_argument("--steps", type=int, default=5000)
    parser.add_argument("--dt", type=float, default=0.02)
    parser.add_argument("--t-initial", type=float, default=1.0)
    parser.add_argument("--t-final", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("lattice_simulation_output.json"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="Directory to write sweep outputs.",
    )
    parser.add_argument(
        "--progress-interval",
        type=float,
        default=30.0,
        help="Seconds between progress updates.",
    )
    parser.add_argument(
        "--sweep",
        action="store_true",
        help="Run a parameter sweep over lambda values and seeds.",
    )
    parser.add_argument(
        "--lambda-values",
        type=str,
        default="0.5,1.0,2.0,3.0,5.0",
        help="Comma-separated lambda values for sweep.",
    )
    parser.add_argument(
        "--seed-start",
        type=int,
        default=1000,
        help="Starting seed for sweep.",
    )
    parser.add_argument(
        "--seed-count",
        type=int,
        default=10,
        help="Number of seeds to cycle through in sweep.",
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=50,
        help="Maximum number of sweep runs to execute.",
    )
    parser.add_argument(
        "--max-runtime",
        type=float,
        default=3600.0,
        help="Maximum runtime for sweep in seconds.",
    )
    args = parser.parse_args()
    return SimulationParams(
        size=args.size,
        hardness=args.hardness,
        steps=args.steps,
        dt=args.dt,
        t_initial=args.t_initial,
        t_final=args.t_final,
        seed=args.seed,
        output_path=args.output,
        progress_interval=args.progress_interval,
    )


def wrap_phase(delta: np.ndarray) -> np.ndarray:
    return (delta + np.pi) % (2 * np.pi) - np.pi


def initialize_field(size: int, rng: np.random.Generator) -> np.ndarray:
    theta = rng.uniform(0.0, 2.0 * np.pi, size=(size, size, size))
    rho = rng.rayleigh(scale=0.1, size=(size, size, size))
    return rho * np.exp(1j * theta)


def compute_laplacian(field: np.ndarray) -> np.ndarray:
    return (
        np.roll(field, 1, axis=0)
        + np.roll(field, -1, axis=0)
        + np.roll(field, 1, axis=1)
        + np.roll(field, -1, axis=1)
        + np.roll(field, 1, axis=2)
        + np.roll(field, -1, axis=2)
        - 6.0 * field
    )


def simulate(params: SimulationParams) -> np.ndarray:
    rng = np.random.default_rng(params.seed)
    field = initialize_field(params.size, rng)
    start = time.perf_counter()
    last_report = start
    for step in range(params.steps):
        if params.steps > 1:
            temperature = params.t_initial + (
                params.t_final - params.t_initial
            ) * step / (params.steps - 1)
        else:
            temperature = params.t_final
        laplacian = compute_laplacian(field)
        potential_force = -2.0 * params.hardness * (
            np.abs(field) ** 2 - 1.0
        ) * field
        noise = (
            rng.normal(size=field.shape) + 1j * rng.normal(size=field.shape)
        )
        field = (
            field
            + params.dt * (laplacian + potential_force)
            + math.sqrt(2.0 * temperature * params.dt) * noise
        )
        now = time.perf_counter()
        if params.progress_interval > 0 and (now - last_report) >= params.progress_interval:
            elapsed = now - start
            pct = 100.0 * (step + 1) / params.steps
            print(
                f"[progress] step {step + 1}/{params.steps} ({pct:.1f}%) "
                f"T={temperature:.3f} elapsed={elapsed:.1f}s",
                flush=True,
            )
            last_report = now
    return field


def compute_energy(field: np.ndarray, hardness: float) -> float:
    energy = 0.0
    for axis in range(3):
        shifted = np.roll(field, -1, axis=axis)
        energy += np.sum(np.abs(shifted - field) ** 2)
    energy += np.sum(hardness * (np.abs(field) ** 2 - 1.0) ** 2)
    return float(energy)


def plaquette_winding_xy(theta: np.ndarray) -> np.ndarray:
    theta_x = np.roll(theta, -1, axis=0)
    theta_y = np.roll(theta, -1, axis=1)
    theta_xy = np.roll(theta_x, -1, axis=1)
    d1 = wrap_phase(theta_x - theta)
    d2 = wrap_phase(theta_xy - theta_x)
    d3 = wrap_phase(theta_y - theta_xy)
    d4 = wrap_phase(theta - theta_y)
    return d1 + d2 + d3 + d4


def plaquette_winding_yz(theta: np.ndarray) -> np.ndarray:
    theta_y = np.roll(theta, -1, axis=1)
    theta_z = np.roll(theta, -1, axis=2)
    theta_yz = np.roll(theta_y, -1, axis=2)
    d1 = wrap_phase(theta_y - theta)
    d2 = wrap_phase(theta_yz - theta_y)
    d3 = wrap_phase(theta_z - theta_yz)
    d4 = wrap_phase(theta - theta_z)
    return d1 + d2 + d3 + d4


def plaquette_winding_zx(theta: np.ndarray) -> np.ndarray:
    theta_z = np.roll(theta, -1, axis=2)
    theta_x = np.roll(theta, -1, axis=0)
    theta_zx = np.roll(theta_z, -1, axis=0)
    d1 = wrap_phase(theta_z - theta)
    d2 = wrap_phase(theta_zx - theta_z)
    d3 = wrap_phase(theta_x - theta_zx)
    d4 = wrap_phase(theta - theta_x)
    return d1 + d2 + d3 + d4


def winding_to_edges(
    winding: np.ndarray,
    orientation: str,
) -> Iterable[Tuple[Tuple[int, int, int], Tuple[int, int, int]]]:
    size = winding.shape[0]
    winding_number = np.rint(winding / (2 * np.pi)).astype(int)
    winding_number[np.abs(winding) < np.pi] = 0
    indices = np.argwhere(winding_number != 0)
    for i, j, k in indices:
        magnitude = abs(int(winding_number[i, j, k]))
        if magnitude == 0:
            continue
        if orientation == "xy":
            node_a = (i, j, (k - 1) % size)
            node_b = (i, j, k)
        elif orientation == "yz":
            node_a = ((i - 1) % size, j, k)
            node_b = (i, j, k)
        elif orientation == "zx":
            node_a = (i, (j - 1) % size, k)
            node_b = (i, j, k)
        else:
            raise ValueError(f"Unknown orientation: {orientation}")
        for _ in range(magnitude):
            yield node_a, node_b


def build_graph(edges: Iterable[Tuple[Tuple[int, int, int], Tuple[int, int, int]]]) -> Dict[
    Tuple[int, int, int], List[Tuple[int, int, int]]
]:
    adjacency: Dict[Tuple[int, int, int], List[Tuple[int, int, int]]] = {}
    for node_a, node_b in edges:
        adjacency.setdefault(node_a, []).append(node_b)
        adjacency.setdefault(node_b, []).append(node_a)
    return adjacency


def analyze_components(
    adjacency: Dict[Tuple[int, int, int], List[Tuple[int, int, int]]]
) -> Tuple[List[List[Tuple[int, int, int]]], List[List[Tuple[int, int, int]]], int]:
    visited = set()
    loops = []
    junction_clusters = []
    total_edges = 0
    
    for node in adjacency:
        if node in visited:
            continue
        stack = [node]
        component_nodes: List[Tuple[int, int, int]] = []
        visited.add(node)
        while stack:
            current = stack.pop()
            component_nodes.append(current)
            for neighbor in adjacency[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        
        degrees = [len(adjacency[n]) for n in component_nodes]
        component_edges = sum(degrees) // 2
        total_edges += component_edges
        
        if all(degree == 2 for degree in degrees):
            loops.append(component_nodes)
        else:
            junction_clusters.append(component_nodes)
            
    return loops, junction_clusters, total_edges


def compute_loop_center(loop: List[Tuple[int, int, int]], size: int) -> np.ndarray:
    # Handle periodic boundary conditions for center calculation
    # We unwrap coordinates relative to the first point
    points = np.array(loop)
    ref = points[0]
    # Calculate offset relative to ref, handling periodicity
    diffs = points - ref
    diffs = (diffs + size // 2) % size - size // 2
    unwrapped = ref + diffs
    return np.mean(unwrapped, axis=0)


def cluster_loops(
    loops: List[List[Tuple[int, int, int]]], 
    threshold: float,
    size: int
) -> Tuple[List[List[Tuple[int, int, int]]], List[List[List[Tuple[int, int, int]]]]]:
    """Cluster loops based on minimum distance between their constituent points."""
    if not loops:
        return [], []
        
    n = len(loops)
    # 1. Compute centers for coarse filtering (optional, but good for large N)
    centers = np.array([compute_loop_center(l, size) for l in loops])
    
    # Simple adjacency matrix for clusters
    adj = {i: [] for i in range(n)}
    
    # Brute force distance check (optimize if N > 1000)
    # We consider two loops connected if ANY point in loop A is close to ANY point in loop B
    # To speed this up, we use the center distance first, then detailed check if close
    
    # Pre-convert all loops to numpy arrays for faster broadcasting
    loop_points = [np.array(l) for l in loops]
    
    for i in range(n):
        for j in range(i + 1, n):
            # Fast check: Center distance minus radii
            # For simplicity in this script, we'll just check min dist between points
            # This is O(N^2 * M^2) worst case, where M is loop length.
            # Given M ~ 100, N ~ 10-50, this is acceptable.
            
            # Vectorized distance calculation with PBC
            # Expand dims to form (M_i, 1, 3) and (1, M_j, 3)
            pts_i = loop_points[i][:, np.newaxis, :]
            pts_j = loop_points[j][np.newaxis, :, :]
            
            # Distances with PBC
            diff = np.abs(pts_i - pts_j)
            diff = np.minimum(diff, size - diff)
            dist_sq = np.sum(diff**2, axis=2)
            min_dist = np.sqrt(np.min(dist_sq))
            
            if min_dist < threshold:
                adj[i].append(j)
                adj[j].append(i)
    
    # Find connected components of loops
    visited = set()
    isolated_loops = []
    clustered_loops_groups = []
    
    for i in range(n):
        if i in visited:
            continue
        
        stack = [i]
        visited.add(i)
        group_indices = []
        
        while stack:
            curr = stack.pop()
            group_indices.append(curr)
            for neighbor in adj[curr]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        
        if len(group_indices) == 1:
            isolated_loops.append(loops[group_indices[0]])
        else:
            group = [loops[idx] for idx in group_indices]
            clustered_loops_groups.append(group)
            
    return isolated_loops, clustered_loops_groups


def extract_vortex_statistics(field: np.ndarray, hardness: float) -> Dict[str, object]:
    theta = np.angle(field)
    winding_xy = plaquette_winding_xy(theta)
    winding_yz = plaquette_winding_yz(theta)
    winding_zx = plaquette_winding_zx(theta)

    edges = list(winding_to_edges(winding_xy, "xy"))
    edges.extend(winding_to_edges(winding_yz, "yz"))
    edges.extend(winding_to_edges(winding_zx, "zx"))

    adjacency = build_graph(edges)
    loops, junction_clusters, total_edges = analyze_components(adjacency)
    
    # Calculate derived lengths
    loop_lengths = [len(l) for l in loops] # Actually edges = nodes for loops
    junction_lengths = [len(c) for c in junction_clusters] # Approx edges
    
    # CLUSTERING LOGIC
    # Interaction range is approx screening length ell
    # beta_pot = hardness (from V = hardness * (|w|^2 - 1)^2)
    # ell = 1.0 / sqrt(2 * beta_pot) = 1.0 / sqrt(2 * hardness)
    # We define threshold as 2.0 * ell (touching cores)
    ell = 1.0 / math.sqrt(2.0 * hardness) if hardness > 0 else 1.0
    threshold = 2.0 * ell
    
    size = field.shape[0]
    isolated, clustered_groups = cluster_loops(loops, threshold, size)
    
    # Statistics
    len_isolated = sum(len(l) for l in isolated)
    len_clustered_loops = sum(sum(len(l) for l in group) for group in clustered_groups)
    len_junctions = sum(len(c) for c in junction_clusters)
    
    # Total "Baryonic" length = Junctions + Clustered Loops
    len_baryonic = len_junctions + len_clustered_loops
    # Total "DM" length = Isolated Loops
    len_dm = len_isolated
    
    ratio_dm_baryon = len_dm / len_baryonic if len_baryonic > 0 else None
    
    return {
        "num_loops": len(loops),
        "num_junction_clusters": len(junction_clusters),
        "loop_lengths": loop_lengths,
        "total_vortex_length": total_edges,
        "num_isolated_loops": len(isolated),
        "num_clustered_groups": len(clustered_groups),
        "length_isolated": len_isolated,
        "length_clustered": len_clustered_loops,
        "length_junctions": len_junctions,
        "ratio_dm_baryon": ratio_dm_baryon,
        "mean_loop_length": float(np.mean(loop_lengths)) if loop_lengths else 0.0,
        "std_loop_length": float(np.std(loop_lengths)) if loop_lengths else 0.0,
    }


def build_output(
    params: SimulationParams, field: np.ndarray, stats: Dict[str, object]
) -> Dict[str, object]:
    beta_pot = params.hardness
    ell_over_a = 1.0 / math.sqrt(2.0 * beta_pot) if beta_pot > 0 else 1.0
    output = {
        "simulation": "LatticeFormation",
        "version": "1.0",
        "params": {
            "N": params.size,
            "lambda": params.hardness,
            "T_initial": params.t_initial,
            "T_final": params.t_final,
            "steps": params.steps,
            "seed": params.seed,
        },
        "canonical_params": {
            "alpha_grad": 1.0,
            "beta_pot": beta_pot,
            "gamma": params.hardness,
            "ell_over_a": ell_over_a,
        },
        "outputs": stats,
        "diagnostics": {
            "final_energy": compute_energy(field, params.hardness),
            "reconnection_events": 0,
        },
    }
    return output


def _parse_lambda_values(raw: str) -> List[float]:
    values: List[float] = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        values.append(float(token))
    return values


def run_single(params: SimulationParams, output_path: Path) -> None:
    field = simulate(params)
    stats = extract_vortex_statistics(field, params.hardness)
    output = build_output(params, field, stats)
    output_path.write_text(json.dumps(output, indent=2))
    print(f"Wrote results to {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Lattice formation and vortex loop statistics simulation.",
    )
    parser.add_argument("--size", type=int, default=64, help="Lattice size N.")
    parser.add_argument("--lambda", dest="hardness", type=float, default=1.0)
    parser.add_argument("--steps", type=int, default=5000)
    parser.add_argument("--dt", type=float, default=0.02)
    parser.add_argument("--t-initial", type=float, default=1.0)
    parser.add_argument("--t-final", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("lattice_simulation_output.json"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="Directory to write sweep outputs.",
    )
    parser.add_argument(
        "--progress-interval",
        type=float,
        default=30.0,
        help="Seconds between progress updates.",
    )
    parser.add_argument(
        "--sweep",
        action="store_true",
        help="Run a parameter sweep over lambda values and seeds.",
    )
    parser.add_argument(
        "--lambda-values",
        type=str,
        default="0.5,1.0,2.0,3.0,5.0",
        help="Comma-separated lambda values for sweep.",
    )
    parser.add_argument(
        "--seed-start",
        type=int,
        default=1000,
        help="Starting seed for sweep.",
    )
    parser.add_argument(
        "--seed-count",
        type=int,
        default=10,
        help="Number of seeds to cycle through in sweep.",
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=50,
        help="Maximum number of sweep runs to execute.",
    )
    parser.add_argument(
        "--max-runtime",
        type=float,
        default=3600.0,
        help="Maximum runtime for sweep in seconds.",
    )
    args = parser.parse_args()

    if not args.sweep:
        params = SimulationParams(
            size=args.size,
            hardness=args.hardness,
            steps=args.steps,
            dt=args.dt,
            t_initial=args.t_initial,
            t_final=args.t_final,
            seed=args.seed,
            output_path=args.output,
            progress_interval=args.progress_interval,
        )
        run_single(params, args.output)
        return

    lambda_values = _parse_lambda_values(args.lambda_values)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    start_time = time.perf_counter()
    run_idx = 0
    for lam in lambda_values:
        for offset in range(args.seed_count):
            if run_idx >= args.runs:
                return
            elapsed = time.perf_counter() - start_time
            if elapsed >= args.max_runtime:
                print(
                    f"[sweep] Reached max runtime {args.max_runtime:.0f}s, stopping.",
                    flush=True,
                )
                return
            seed = args.seed_start + offset
            params = SimulationParams(
                size=args.size,
                hardness=lam,
                steps=args.steps,
                dt=args.dt,
                t_initial=args.t_initial,
                t_final=args.t_final,
                seed=seed,
                output_path=args.output,
                progress_interval=args.progress_interval,
            )
            output_path = args.output_dir / f"lattice_run_{run_idx:03d}_lam_{lam:.3f}_seed_{seed}.json"
            print(
                f"[sweep] Run {run_idx + 1}/{args.runs} "
                f"lambda={lam:.3f} seed={seed}",
                flush=True,
            )
            run_single(params, output_path)
            run_idx += 1


if __name__ == "__main__":
    main()
