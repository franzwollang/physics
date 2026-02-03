# Knot Stability Simulation

This script implements the effective-string minimization described in
`soliton_model_knot_theory_analysis.md`. It initializes a trefoil knot,
minimizes the discrete energy functional, and writes structured JSON outputs
for a scan of bending stiffness values.

## Usage

```bash
python simulations/knot_stability_simulation.py --kB 0.1 1.0 10.0
# Optional: cap optimization evaluations for quicker runs
python simulations/knot_stability_simulation.py --kB 1.0 --maxiter 50 --maxfun 2000
# Optional: switch optimizer (Powell can be more stable but slower)
python simulations/knot_stability_simulation.py --kB 1.0 --method Powell --maxiter 200
# Optional: adjust L-BFGS-B line search/gradient tolerances
python simulations/knot_stability_simulation.py --kB 1.0 --gtol 1e-5 --maxls 80
# Example: faster convergent run for sanity-checking
python simulations/knot_stability_simulation.py --kB 1.0 --n-beads 30 --kT 10 --epsilon 1.0 --maxiter 200
```

The outputs are written to `simulations/results/` by default. Each run produces
one JSON file per `k_B` value:

```
simulations/results/knot_stability_kB_0.10.json
simulations/results/knot_stability_kB_1.00.json
simulations/results/knot_stability_kB_10.00.json
```

## Notes

- The script uses a Lennard-Jones repulsion to discourage self-intersections.
- Topology preservation is estimated by checking the minimum non-adjacent
  bead distance against the core radius.
- If `pyknotid` is available, a writhe estimate is included in the output.
- Optimization metadata (`n_iterations`, `n_function_evals`, and the optimizer
  message) are included in the JSON outputs to help diagnose convergence.
- Use `--method Powell` if L-BFGS-B stalls or returns an abnormal termination.
- `--gtol` and `--maxls` can help the line search converge when the Lennard-Jones
  term makes the landscape stiff.
- If `--maxfun` is set and the evaluation budget is exhausted, the optimizer will
  report a non-success status with a message indicating the evaluation limit was hit.
- Lowering `n_beads` or the repulsion strength can be useful for quick convergence
  checks before running higher-resolution scans.
