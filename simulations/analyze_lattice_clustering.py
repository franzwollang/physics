import json
import statistics
import sys
from pathlib import Path
from typing import Dict, List


def analyze_results(results_dir: Path) -> None:
    paths = sorted(results_dir.glob("lattice_run_*.json"))
    if not paths:
        print(f"No results found in {results_dir}")
        return

    # Group by lambda
    by_lambda: Dict[float, List[Dict]] = {}

    print(f"Found {len(paths)} run files.")

    for p in paths:
        try:
            data = json.loads(p.read_text())
            lam = data["params"]["lambda"]
            outputs = data["outputs"]

            if lam not in by_lambda:
                by_lambda[lam] = []

            by_lambda[lam].append(outputs)
        except Exception as e:
            print(f"Error reading {p}: {e}")

    print(f"\nSummary Statistics (Clustering Analysis):")
    print(
        f"{'Lambda':<8} | {'Beta':<6} | {'Ell':<6} | {'Runs':<4} | {'Ratio(I/C)':<12} | {'Iso Length':<12} | {'Clus Length':<12} | {'# Iso':<8} | {'# Clus':<8}"
    )
    print("-" * 110)

    for lam in sorted(by_lambda.keys()):
        runs = by_lambda[lam]
        n = len(runs)

        # Extract canonical params from first run (assume constant for fixed lambda)
        canonical = runs[0].get("canonical_params", {}) if n > 0 else {}
        # In the new format, canonical_params is at root level if I saved it right?
        # Wait, the JSON structure is:
        # { "params": {...}, "canonical_params": {...}, "outputs": {...} }
        # My script reads data = json.loads(), then outputs = data['outputs'].
        # I need to read canonical_params from data.

        # RE-READING LOGIC NEEDED
        pass

    # Re-loop to get full data structure
    by_lambda_full = {}
    for p in paths:
        try:
            data = json.loads(p.read_text())
            lam = data["params"]["lambda"]
            if lam not in by_lambda_full:
                by_lambda_full[lam] = []
            by_lambda_full[lam].append(data)
        except:
            pass

    for lam in sorted(by_lambda_full.keys()):
        runs = by_lambda_full[lam]
        n = len(runs)
        first = runs[0]
        canon = first.get("canonical_params", {})
        beta = canon.get("beta_pot", 0.0)
        ell = canon.get("ell_over_a", 0.0)

        ratios = []
        iso_lens = []
        clus_lens = []
        num_isos = []
        num_clus = []

        for d in runs:
            out = d["outputs"]
            if out.get("ratio_dm_baryon") is not None:
                ratios.append(out["ratio_dm_baryon"])
            iso_lens.append(out.get("length_isolated", 0))
            clus_lens.append(
                out.get("length_clustered", 0) + out.get("length_junctions", 0)
            )
            num_isos.append(out.get("num_isolated_loops", 0))
            num_clus.append(
                out.get("num_clustered_groups", 0) + out.get("num_junction_clusters", 0)
            )

        avg_ratio = statistics.mean(ratios) if ratios else 0.0
        avg_iso = statistics.mean(iso_lens)
        avg_clus = statistics.mean(clus_lens)
        avg_n_iso = statistics.mean(num_isos)
        avg_n_clus = statistics.mean(num_clus)

        print(
            f"{lam:<8.3f} | {beta:<6.3f} | {ell:<6.2f} | {n:<4} | {avg_ratio:<12.3f} | {avg_iso:<12.1f} | {avg_clus:<12.1f} | {avg_n_iso:<8.1f} | {avg_n_clus:<8.1f}"
        )


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python analyze_lattice_clustering.py <results_dir>")
        sys.exit(1)

    analyze_results(Path(sys.argv[1]))
