# Simulation Task: Lattice Formation & Statistics of Phase-Dark Solitons

---

## 0. Parameter Interface (Canonical Microphysics)

This simulation is part of a unified constraint program. All simulations share a common set of **Canonical Microphysics Parameters** derived from the underlying complex scalar field theory.

### 0.1 Canonical Parameters (Shared Across All Simulations)

| Symbol | Name | Description |
|--------|------|-------------|
| $\alpha_{\rm grad}$ | Gradient Stiffness | Coefficient of $\|\nabla\Psi\|^2$ in the free energy. Sets the energy cost of spatial variation. |
| $\beta_{\rm pot}$ | Potential Curvature | Coefficient of $-w^2$ in the Mexican-hat potential. Controls symmetry breaking depth. |
| $\gamma$ | Quartic Coupling | Coefficient of $w^4$. Stabilizes the potential at large amplitude. |
| $\kappa$ | Phase Stiffness | Coefficient of $w^2\|\nabla\phi\|^2$. Sets the phase-wave speed. |

### 0.2 Derived Physical Scales

| Symbol | Formula | Interpretation |
|--------|---------|----------------|
| $w_*$ | $\sqrt{\beta_{\rm pot} / 2\gamma}$ | Vacuum amplitude (spontaneous symmetry breaking scale). |
| $m_\xi$ | $\sqrt{2\beta_{\rm pot}}$ | Amplitude gap / mass (energy cost to excite amplitude mode). |
| $\ell$ | $\sqrt{\alpha_{\rm grad} / 2\beta_{\rm pot}}$ | Correlation length / vortex core width. |
| $c_s$ | $\sqrt{\kappa} w_*$ | Phase wave speed (effective "speed of light"). |

### 0.3 Mapping to This Simulation's Parameters

This simulation uses a **rescaled lattice Hamiltonian** where:
- Lattice spacing $a = 1$ (all lengths in units of $a$).
- Field rescaled so $w_* = 1$.

The dimensionless **hardness parameter** $\lambda$ used in this simulation maps to the canonical parameters as:
\[
\lambda = \frac{\gamma}{\alpha_{\rm grad}} \cdot a^2
\]
In the natural units where $a = 1$ and $\alpha_{\rm grad} = 1$:
\[
\lambda = \gamma, \qquad \ell = \frac{1}{\sqrt{2\beta_{\rm pot}}}
\]
**Constraint:** The lattice must resolve the core: $\ell \ge 2a$, i.e., $\beta_{\rm pot} \le 0.125$ in lattice units.

---

## 1. Theoretical Context & Motivation

This project investigates a **complex scalar field theory** ($\Psi \in \mathbb{C}$) undergoing a symmetry-breaking phase transition. We posit that in the early universe, the field relaxes from a high-temperature disordered state to a low-temperature ordered vacuum ("mexican hat" potential).

During this process (the Kibble-Zurek mechanism), topological defects form. In a standard $U(1)$ theory, these are cosmic strings (vortex lines). However, our framework modifies the standard picture in two ways:

1.  **Finite Core Size:** The field amplitude $w = |\Psi|$ is dynamic and gapped, giving defects a "soft" core of characteristic width $\ell$.
2.  **Topological Complexity:** We are interested in the ratio of **simple closed loops** (knots/unknots) to **complex junctions** (where strings branch or intersect, effectively forming a network).

**Goal:** Determine the relic abundance ratio $R = N_{\text{loops}} / N_{\text{junctions}}$ and the loop length distribution after the system freezes out. This ratio is a critical predictor for Dark Matter (loops) vs. Baryonic Matter (junctions) in this specific emergent-gravity framework.

---

## 2. The Mathematical Model

We model the complex order parameter $\Psi(x) = w(x) e^{i\phi(x)}$ on a discrete 3D Euclidean lattice.

### 2.1 The Continuum Free Energy

The system evolves to minimize the effective free energy density:
\[
\mathcal{F}[\Psi] = \alpha_{\rm grad} |\nabla \Psi|^2 + V(|\Psi|)
\]
where the local potential $V$ drives symmetry breaking:
\[
V(w) = -\beta_{\rm pot} w^2 + \gamma w^4
\]

### 2.2 Lattice Discretization (Simulation Units)

We discretize on an $N^3$ cubic lattice with spacing $a=1$.
Let $\Psi_{\mathbf{x}}$ be the complex value at site $\mathbf{x} = (i,j,k)$.
The Hamiltonian $H$ (Energy) is calculated using nearest-neighbor differences:
\[
H = \sum_{\mathbf{x}} \sum_{\mu=1}^3 \left| \Psi_{\mathbf{x}+\hat{\mu}} - \Psi_{\mathbf{x}} \right|^2 + \sum_{\mathbf{x}} \lambda (|\Psi_{\mathbf{x}}|^2 - 1)^2
\]

**Boundary Conditions:** Periodic Boundary Conditions (PBC) in all 3 directions. $\Psi_{(N, j, k)} \equiv \Psi_{(0, j, k)}$, etc.

---

## 3. Simulation Protocol

### 3.1 Initialization (High Temperature)

Start with a "hot" randomized state to mimic the pre-geometric era:

-   $\Psi_{\mathbf{x}} = \rho_{\mathbf{x}} e^{i \theta_{\mathbf{x}}}$
-   $\theta_{\mathbf{x}} \sim \text{Uniform}[0, 2\pi)$
-   $\rho_{\mathbf{x}} \sim \text{Rayleigh}(\sigma=0.1)$ (Small fluctuating amplitudes).

### 3.2 Cooling / Annealing (The Quench)

Evolve using **Overdamped Langevin Dynamics** (Gradient Flow with noise).
Update Rule (Euler-Maruyama):
\[
\Psi_{\mathbf{x}}(t+\Delta t) = \Psi_{\mathbf{x}}(t) - \Delta t \frac{\partial H}{\partial \Psi_{\mathbf{x}}^*} + \sqrt{2 T(t) \Delta t} \, \eta_{\mathbf{x}}
\]

-   **Gradient:** The force term is the discrete Laplacian plus the local potential derivative:
    \[ -\frac{\partial H}{\partial \Psi_{\mathbf{x}}^*} = \sum_{\mu} (\Psi_{\mathbf{x}+\hat{\mu}} + \Psi_{\mathbf{x}-\hat{\mu}} - 2\Psi_{\mathbf{x}}) - 2\lambda (|\Psi_{\mathbf{x}}|^2 - 1) \Psi_{\mathbf{x}} \]
-   **Noise:** $\eta_{\mathbf{x}}$ is a complex Gaussian noise $\mathcal{N}(0,1) + i\mathcal{N}(0,1)$.
-   **Schedule:**
    -   Start $T = 2.0$ (High).
    -   Cool to $T = 0.0$ over $N_{\text{steps}}$ (e.g., 5000 steps).

### 3.3 Freeze-out

The simulation ends when the topological network stabilizes. This occurs when $T(t) \ll \Delta E_{\text{core}}$.

---

## 4. Analysis & Observables

### 4.1 Defect Detection (Plaquette Winding)

Detect vortex lines by calculating phase winding on every plaquette (face) of the lattice.
For a plaquette in the $xy$-plane at $(i,j,k)$:

1.  Extract phases: $\theta_{00}, \theta_{10}, \theta_{11}, \theta_{01}$.
2.  Compute phase differences wrapped to $(-\pi, \pi]$:
    \[ \Delta_{\mu} \theta = (\theta_{\text{next}} - \theta_{\text{curr}} + \pi \pmod{2\pi}) - \pi \]
3.  Sum differences: $W = \sum \Delta \theta$.
4.  If $W \approx \pm 2\pi$, a vortex line passes through this face in the $z$-direction.

### 4.2 Topological Tracing (The Dual Lattice)

-   Create a list of "pierced faces" and their orientations.
-   Connect these piercings to form lines on the **dual lattice** (centers of the cubes).
-   **Algorithm:**
    -   Start at a pierced face.
    -   The flux must exit the current dual-cell through another face (Conservation of topological charge, $\nabla \cdot J = 0$).
    -   Follow the flux until it closes (Loop) or hits a messy node (Junction).

### 4.3 Classification

-   **Isolated Loop:** A connected component where every node has degree 2.
-   **Junction/Network:** A connected component containing nodes with degree $\neq 2$.

---

## 5. Parameter Study Grid

1.  **Lattice Size:** $N \in [64, 128]$.
2.  **Hardness:** $\lambda \in [1.0, 5.0]$. (Standard Higgs $\lambda \approx 1$).
3.  **Cooling:** `steps` $\in [0 \text{ (quench)}, 10^4 \text{ (anneal)}]$.

---

## 6. Implementation Specs

-   **Language:** Python + Numba (for the update loop) OR C++.
-   **Libraries:** `numpy`, `scipy.ndimage` (for potential cluster labeling).
-   **Vectorization:** Use `np.roll` to compute neighbor differences efficiently in Python.
    ```python
    laplacian = np.roll(psi, 1, 0) + np.roll(psi, -1, 0) + ... - 6*psi
    ```

---

## 7. Output Interface (Structured Results)

All runs must output a **JSON file** with the following schema:

```json
{
  "simulation": "LatticeFormation",
  "version": "1.0",
  "params": {
    "N": 64,
    "lambda": 1.0,
    "T_initial": 2.0,
    "T_final": 0.0,
    "steps": 5000,
    "seed": 12345
  },
  "canonical_params": {
    "alpha_grad": 1.0,
    "beta_pot": 0.5,
    "gamma": 1.0,
    "ell_over_a": 1.0
  },
  "outputs": {
    "num_loops": 42,
    "num_junction_clusters": 8,
    "loop_lengths": [12, 24, 8, ...],
    "total_vortex_length": 350,
    "ratio_R_count": 5.25,
    "ratio_R_length": 0.72,
    "mean_loop_length": 15.2,
    "std_loop_length": 8.1
  },
  "diagnostics": {
    "final_energy": 123.4,
    "reconnection_events": 0
  }
}
```

### 7.1 Key Outputs Explained

| Field | Meaning | Physical Interpretation |
|-------|---------|-------------------------|
| `num_loops` | Count of isolated closed loops | Proxy for number of DM candidates |
| `num_junction_clusters` | Count of networked/junction components | Proxy for baryon-like structures |
| `ratio_R_count` | `num_loops / num_junction_clusters` | Abundance ratio (count-based) |
| `ratio_R_length` | Total loop length / Total vortex length | Mass fraction in loops |
| `loop_lengths` | List of individual loop sizes | Needed for mass spectrum modeling |

---

## 8. Synthesis Role (How This Constrains the Model)

This simulation constrains the **formation channel**: what defect populations emerge from a phase transition given the microphysics.

### 8.1 Target Observable

The primary target is the **DM-to-Baryon mass ratio**:
\[
\frac{\Omega_{\rm DM}}{\Omega_b} \approx 5.3 \quad (\text{observed})
\]
If we assume loops become DM and junctions become baryons, we need:
\[
R_{\rm length} \approx \frac{5.3}{1 + 5.3} \approx 0.84
\]
(i.e., ~84% of total vortex length should be in isolated loops).

### 8.2 Explicit Constraint Function

**For synthesis, report the following constraint:**

$$
f_{\rm Lattice}(\lambda, q) = R_{\rm length}(\lambda, q)
$$

where:
- $\lambda = \gamma / \alpha_{\rm grad}$ is the hardness parameter (from `canonical_params`)
- $q$ is the quench schedule (nuisance parameter)

**Target band:**
$$
R_{\rm length} \in [0.75, 0.90]
$$

**Viable region:** The set of $\lambda$ values such that $R_{\rm length} \in [0.75, 0.90]$ for *at least one* reasonable quench schedule $q$.

### 8.3 Interpretation

- **If $R_{\rm length} \ll 0.8$:** Too few loops. Either the model is wrong, or the freeze-out is too late (loops annihilate). Try faster quench or lower $\lambda$.
- **If $R_{\rm length} \approx 0.8$:** Viable region. Record the parameter band $(\lambda, \text{quench rate})$.
- **If $R_{\rm length} \gg 0.8$:** More DM than expected. Either acceptable (DM abundance is an upper bound from gravity) or suggests junctions are rare (tension with visible matter).

### 8.4 Interface with Other Simulations

- **Knot Simulation:** Takes the typical loop size from `mean_loop_length` and checks if loops of that size are *stable* (metastable radius exists).
- **Poisson Simulation:** Independent; validates the gravity law.
- **Geometry Simulation:** Provides the substrate (3D lattice) on which this field simulation runs. If Geometry fails to produce 3D, this simulation is moot.

---

## 9. Failure Modes & Mitigation

### 9.1 Defect Annihilation

-   **Issue:** Small loops may annihilate too quickly before statistics are gathered, leaving an empty lattice.
-   **Mitigation:** Introduce **Hubble Friction**. Add a damping term $-\gamma \frac{\partial \Psi}{\partial t}$ (implicit in Langevin, but can be boosted). Or, stop the simulation immediately after the phase transition ($T \approx T_c$) rather than waiting for $T=0$.

### 9.2 Resolution Artifacts

-   **Issue:** If the core size $\xi \approx a$ (lattice spacing), the lattice "locks" defects in place (Peierls-Nabarro barrier), preventing physical motion.
-   **Mitigation:** Ensure $\xi \ge 2.0$. Adjust parameters so $\lambda$ is not too large. Soft cores move more physically.

### 9.3 Topological Ambiguity

-   **Issue:** "Junctions" might just be two loops passing very close to each other (within 1 cell).
-   **Mitigation:** Use **Topological Thinning** or Skeletonization on the vortex voxel map before graph tracing. Treat "crossing" lines as junctions only if they persist over multiple time steps.

### 9.4 Periodic Boundary Artifacts

-   **Issue:** Loops may wrap around the box, making "length" ambiguous.
-   **Mitigation:** Use large enough $N$ that most loops fit within $N/2$. Flag/exclude loops with winding number $\neq 0$ around the torus.
