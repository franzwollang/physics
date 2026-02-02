# Analysis Task: Stability & Energetics of Phase-Dark Solitons

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
| $w_*$ | $\sqrt{\beta_{\rm pot} / 2\gamma}$ | Vacuum amplitude. |
| $m_\xi$ | $\sqrt{2\beta_{\rm pot}}$ | Amplitude gap / mass. |
| $\ell$ | $\sqrt{\alpha_{\rm grad} / 2\beta_{\rm pot}}$ | Correlation length / vortex core width. |
| $c_s$ | $\sqrt{\kappa} w_*$ | Phase wave speed. |

### 0.3 Mapping to This Simulation's Parameters

This simulation uses an **Effective String Model** where the vortex tube is approximated as a 1D elastic curve with:

| Simulation Param | Formula from Canonical | Interpretation |
|------------------|------------------------|----------------|
| $a$ (core radius) | $a = \ell = \sqrt{\alpha_{\rm grad} / 2\beta_{\rm pot}}$ | Tube thickness / excluded volume radius. |
| $T$ (tension) | $T \approx \pi \ell^2 \cdot \beta_{\rm pot} w_*^2 = \frac{\pi \alpha_{\rm grad} w_*^2}{2}$ | Energy per unit length of the string. |
| $B$ (bending) | $B \approx T \cdot \ell^2 = \frac{\pi \alpha_{\rm grad}^2 w_*^2}{4 \beta_{\rm pot}}$ | Bending modulus (resists curvature). |
| $\Delta E_{\rm cut}$ (barrier) | $\Delta E \approx \ell^3 \cdot \beta_{\rm pot} w_*^2$ | Energy to "melt" a core-sized region (reconnection barrier). |

**In simulation units** where $\ell = 1$ and $T = 1$:
- $a = 1$
- $B = 1$ (assuming the $T \ell^2$ scaling)
- $\Delta E_{\rm cut} \approx 1$ (order unity in these units)

---

## 1. Theoretical Context

We are evaluating a specific Dark Matter candidate: a **"Phase-Dark" Soliton**.
In our underlying field theory, this corresponds to a **closed, knotted flux tube** of a complex scalar field $\Psi$. Unlike cosmic strings which typically decay or radiate, we propose that **knotted** configurations are topologically protected and energetically metastable in the late-time universe.

**Goal:** Use an "Effective String" model (Polymer Physics approach) to calculate the energy landscape of these knots. We must determine if a **metastable radius $R_*$** exists that prevents the knot from shrinking to zero/dissolving.

---

## 2. The Effective String Model

We approximate the soliton not as a 3D field, but as a discretized chain of $N$ "beads" at positions $\mathbf{r}_i$.

### 2.1 The Discrete Energy Functional

The total energy $E = E_T + E_B + E_C + E_{LJ}$.

1.  **Tension ($E_T$):** Linear spring potential (approximating arc length cost).
    \[ E_T = \frac{k_T}{2} \sum_{i=0}^{N-1} (|\mathbf{r}_{i+1} - \mathbf{r}_i| - l_0)^2 \]
    -   $k_T$: Spring constant (large $k_T$ enforces near-constant segment length).
    -   $l_0$: Equilibrium segment length.

2.  **Bending Stiffness ($E_B$):** Cosine potential on angles.
    \[ E_B = k_B \sum_{i=0}^{N-1} (1 - \cos \theta_i) \]
    -   $\cos \theta_i = \hat{\mathbf{t}}_i \cdot \hat{\mathbf{t}}_{i-1}$, where $\hat{\mathbf{t}}_i$ is the unit tangent vector of segment $i$.
    -   $k_B$: Bending modulus (maps to $B$ above).

3.  **Twist Stiffness ($E_C$):** (Optional for first pass). Cost for discrete torsion.
    -   Requires defining a "material frame" (normal vector $\mathbf{n}_i$) along the curve.
    -   $E_C = k_C \sum (\alpha_i - \tau_0)^2$.

4.  **Self-Contact (Excluded Volume $E_{LJ}$):** Lennard-Jones repulsion to prevent self-intersection (topology violation).
    \[ E_{LJ} = \sum_{|i-j| > 2} 4\epsilon_{\rm LJ} \left[ \left(\frac{\sigma}{r_{ij}}\right)^{12} \right] \]
    -   $\sigma = a$: Effective tube diameter (core width).
    -   Only applied to non-adjacent beads. Cutoff at $r_{ij} = 2.5\sigma$.

### 2.2 Parameter Defaults (Simulation Units)

| Parameter | Default Value | Notes |
|-----------|---------------|-------|
| $a$ | 1.0 | Core width (unit of length). |
| $l_0$ | 1.0 | Segment length (should be $\approx a$). |
| $k_T$ | 100.0 | Large to enforce near-inextensibility. |
| $k_B$ | 1.0 | Bending modulus (scan this). |
| $\sigma$ | 1.0 | LJ diameter = core width. |
| $\epsilon_{\rm LJ}$ | 10.0 | Large to enforce hard-core repulsion. |

---

## 3. Analysis Tasks

### 3.1 Analytical Minimization (The "Unknot")

Derive the condition for stability of a torus of radius $R$.
Energy: $E(R) \approx 2\pi R T + \frac{B}{2} (2\pi R) \frac{1}{R^2} = 2\pi T R + \frac{\pi B}{R}$.

-   Find $R_*$ such that $dE/dR = 0$: $R_* = \sqrt{B / 2T}$.
-   Check if $d^2E/dR^2 > 0$ (Stability): Yes, $d^2E/dR^2 = 2\pi B / R^3 > 0$.

**Key Result:** A circular unknot has a stable radius $R_* = \sqrt{B/2T} = \ell / \sqrt{2}$ in canonical units.

### 3.2 Numerical Minimization (The "Trefoil")

1.  **Initialization:** Generate coordinates for a Trefoil knot ($3_1$) using the parametric equations:
    \[ x = \sin t + 2 \sin 2t, \quad y = \cos t - 2 \cos 2t, \quad z = -\sin 3t \]
    Rescale so loop length $\approx 50a$. Discretize into $N=100$ beads.
2.  **Optimization:** Use `scipy.optimize.minimize` with method `L-BFGS-B` or a custom damped molecular dynamics loop (Verlet integration with friction).
    -   Minimize $E(\mathbf{r}_0, \dots, \mathbf{r}_{N-1})$.
3.  **Result Check:**
    -   Did it collapse to a point? (Unstable).
    -   Did it settle to a finite size geometry? (Metastable).
    -   Calculate the final Radius of Gyration $R_g$.

### 3.3 The Decay Barrier

Estimate the energy required to "cut" the knot.

-   The barrier is the cost to locally melt the core: $\Delta E \approx a^3 \times V_{\max} \approx a^3 \beta_{\rm pot} w_*^2$.
-   Compare the final potential energy of the knot $E_{\text{knot}}(R_*)$ to this barrier.
-   **Stability Criterion:** If $\Delta E_{\rm cut} \gg k_B T_{\rm cosmic}$ (late-time thermal fluctuations), the knot is frozen.

---

## 4. Parameter Study

Run the Trefoil minimization for a grid of bending stiffnesses:

-   $k_B \in [0.1, 1.0, 10.0]$ (Soft vs Stiff tubes).
-   Report $R_g$ vs $k_B$.

---

## 5. Implementation Specs

-   **Language:** Python.
-   **Libraries:** `numpy`, `scipy.optimize`.
-   **Geometry:** Use specific Trefoil coordinates to ensure the correct topology is initialized.
-   **Topology Check:** (Bonus) Calculate the discrete **writhe** or use `pyknotid` library to confirm the knot type hasn't changed during minimization.

---

## 6. Output Interface (Structured Results)

All runs must output a **JSON file** with the following schema:

```json
{
  "simulation": "KnotStability",
  "version": "1.0",
  "params": {
    "knot_type": "Trefoil",
    "N_beads": 100,
    "k_T": 100.0,
    "k_B": 1.0,
    "sigma": 1.0,
    "epsilon_LJ": 10.0
  },
  "canonical_params": {
    "alpha_grad": 1.0,
    "beta_pot": 0.5,
    "gamma": 1.0,
    "ell": 1.0,
    "T_string": 1.57,
    "B_string": 1.57
  },
  "outputs": {
    "converged": true,
    "final_energy": 45.2,
    "radius_of_gyration": 3.1,
    "metastable_radius_R_star": 2.8,
    "topology_preserved": true,
    "final_writhe": -3.2
  },
  "stability_analysis": {
    "barrier_estimate": 0.5,
    "barrier_to_energy_ratio": 0.011,
    "frozen": true
  }
}
```

### 6.1 Key Outputs Explained

| Field | Meaning | Physical Interpretation |
|-------|---------|-------------------------|
| `converged` | Did the optimizer find a minimum? | Basic sanity check. |
| `radius_of_gyration` | RMS distance from center of mass | Size of the relaxed knot. |
| `metastable_radius_R_star` | Theoretical $\sqrt{B/2T}$ for unknot | Comparison scale. |
| `topology_preserved` | Did the knot type survive minimization? | If false, result is invalid (self-intersection occurred). |
| `barrier_estimate` | $\Delta E_{\rm cut}$ in simulation units | Energy to unknot via reconnection. |
| `frozen` | Is $\Delta E_{\rm cut} \gg 1$ (thermal scale)? | Whether the knot survives cosmologically. |

---

## 7. Synthesis Role (How This Constrains the Model)

This simulation constrains the **late-time viability**: are the loops produced by the Lattice simulation *stable* against decay?

### 7.1 Target Observable

The primary target is **collisionless dark matter**: DM particles must have low self-interaction cross-section.
\[
\frac{\sigma_{\rm int}}{m} < 1 \text{ cm}^2/\text{g} \quad (\text{Bullet Cluster bound})
\]
If the knot has radius $R_g$ and mass $M \sim T \cdot L$ (tension × length):
\[
\frac{\sigma_{\rm int}}{m} \sim \frac{\pi R_g^2}{T \cdot L}
\]

### 7.2 Explicit Constraint Functions

**For synthesis, report the following constraints:**

**Constraint A (Stability):**
$$
f_{\rm Knot}^{(a)}(k_B) = R_g(k_B) \quad \Rightarrow \quad R_g > 0
$$

**Constraint B (Collisionless):**
$$
f_{\rm Knot}^{(b)}(k_B) = \frac{\sigma}{m}(k_B) \quad \Rightarrow \quad \frac{\sigma}{m} < 1 \text{ cm}^2/\text{g}
$$

where:
- $k_B$ is the bending modulus, related to canonical parameters via: $k_B \propto B \propto \alpha_{\rm grad}^2 / \beta_{\rm pot}$

**Viable region:** The set of $k_B$ (equivalently $\alpha_{\rm grad}/\beta_{\rm pot}$) where both constraints hold simultaneously.

**Critical threshold:** Define $k_B^{\rm crit}$ as the minimum bending modulus for which $R_g > 0$.

### 7.3 Interpretation

- **If $R_g \to 0$ (collapse):** No stable knot exists. DM candidate ruled out for this $k_B$.
- **If $R_g > 0$ and `frozen=true`:** Knot is cosmologically stable. Compute $\sigma/m$ and check against bound.
- **If $R_g$ is too large:** Knot is "fluffy" → $\sigma/m$ exceeds bound → ruled out.

### 7.4 Interface with Other Simulations

- **Lattice Simulation:** Provides typical loop length $L$ from `mean_loop_length`. This simulation checks if loops of that length are stable and estimates their mass.
- **Poisson Simulation:** Independent; validates the gravity law.
- **Geometry Simulation:** Provides the substrate. If spectral dimension $\neq 3$, the effective string model assumptions may break.

---

## 8. Failure Modes & Mitigation

### 8.1 Knot Collapse

-   **Issue:** The discretized knot collapses to a single point because the discrete Bending Energy isn't strong enough to oppose Tension at the scale of 1 segment length.
-   **Mitigation:** Add a "Hard Core" constraint $r_{ij} > \sigma$ for *adjacent* beads too (or stiff springs). Or, enforce a fixed total length constraint (inextensible string) and only minimize Bending energy.

### 8.2 Local Minima

-   **Issue:** The optimizer gets stuck in a "crumpled" local minimum that isn't the open trefoil.
-   **Mitigation:** Use **Simulated Annealing** (Metropolis-Hastings) to jiggle the chain out of local traps before doing the final gradient descent. Start with a "inflated" knot (large radius) and let it shrink slowly.

### 8.3 Crossing (Topology Violation)

-   **Issue:** The "Excluded Volume" potential is too soft, and the knot passes through itself, becoming an Unknot.
-   **Mitigation:** Use a very hard LJ potential (power 12 or higher). Check the **Crossing Number** invariant periodically. If it changes, the run is invalid. Use the `pyknotid` library for robust knot identification.

### 8.4 Discretization Error

-   **Issue:** Too few beads ($N$) leads to inaccurate curvature estimates.
-   **Mitigation:** Run convergence test: increase $N$ until $R_g$ stabilizes (varies by $<5\%$).
