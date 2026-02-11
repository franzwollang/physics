# Lattice Formation Simulation (Phase-Dark Solitons)

This document implements the simulation protocol described in
`soliton_model_lattice_simulation_plan.md`. It provides a concrete workflow
for running the lattice quench, detecting vortex lines, and extracting the
loop/junction statistics that feed the synthesis pipeline.

## Quick Start

```bash
pip install numpy

python simulations/lattice/lattice_soliton_simulation.py \
  --size 64 \
  --lambda 0.05 \
  --steps 5000 \
  --dt 0.02 \
  --t-initial 1.0 \
  --t-final 0.0 \
  --seed 12345 \
  --output lattice_simulation_output.json
```

The command produces a JSON file with the schema specified in the plan, including
counts of loop components, junction clusters, and the loop length distribution.

## Implementation Notes

- **Initialization** uses Rayleigh-distributed amplitudes (`sigma=0.1`) and
  uniformly random phases on `[0, 2π)`.
- **Dynamics** follow the Euler-Maruyama update rule with a linear temperature
  quench from `T_initial` to `T_final`.
- **Vortex detection** computes plaquette phase winding in the `xy`, `yz`, and
  `zx` planes, and maps nonzero winding numbers to edges on the dual lattice.
- **Topology classification** (Updated):
    - Uses **Proximity Clustering** to group loops within distance $R \approx 2\ell$.
    - Isolated loops -> Dark Matter.
    - Clustered loops -> Baryonic Matter.

## Output

The `outputs` field includes:

- `num_isolated_loops`
- `num_clustered_groups`
- `ratio_dm_baryon` (Length ratio)
- `loop_lengths`
- `total_vortex_length`

## Parameter Guidance

The core size $\ell \approx 1/\sqrt{2\lambda}$.
- `lambda=0.02` -> `ell=5.0` (Large, overlapping cores -> Baryon dominated).
- `lambda=0.12` -> `ell=2.0` (Smaller, resolved cores).
- `lambda=1.0` -> `ell=0.7` (Too small for lattice, leads to pinning/aliasing).

Recommended sweep range: `[0.02, 0.20]`.

## Preliminary Calibration (Feb 2026)

- **Lambda=0.02**: Ratio(DM/B) ~ 0.01. (Dense soup, highly clustered).
- **Lambda=0.20**: Ratio(DM/B) -> 100% Isolated (if density drops enough).
- **Target**: Find lambda where Ratio ~ 5.3.
