# Simulation Task: Primordial Yield via The Inverse Problem

---

## 0. Parameter Interface (Canonical Microphysics)

This simulation suite solves the **Inverse Problem** of Dark Matter genesis. Instead of running a massive cosmological simulation and hoping for the right result, we break the chain into three testable links to determine *what kind of primordial noise* is required to produce the observed universe.

### 0.1 Canonical Parameters

| Symbol | Name | Description |
|--------|------|-------------|
| $\alpha_{\rm grad}$ | Gradient Stiffness | Energetic cost of domain walls. |
| $\beta_{\rm pot}$ | Potential Curvature | Depth of the symmetry-breaking well. |
| $\kappa$ | Phase Stiffness | Emergent from local graph connectivity. |
| $\vec{v}_{phase}$ | Phase Wind | The local gradient of the phase field ($\nabla \phi$) representing noise/flows. |

---

## 1. Theoretical Context: The Three-Stage Chain

We hypothesize that Baryonic Matter ($l=1$) forms from scalar lumps ($l=0$) that are "spun up" by the ambient phase noise via the Hairy Ball Theorem. Dark Matter is simply the population of lumps that failed to spin up.

**The Chain of Causality:**
1.  **Induction (Micro-Physics):** A scalar lump sits in a noisy phase bath ("Phase Wind"). If the wind is strong/coherent enough, it forces the lump to acquire a topological defect ($l=0 \to l=1$) to satisfy the Hairy Ball boundary condition.
2.  **Statistics (The Spectrum):** The universe is filled with a spectrum of lump sizes and noise intensities (e.g., Rogue Wave statistics).
3.  **Origin (Topology):** This spectrum originates from the quench on the Infinite-Clique Graph.

**The Strategy:**
Work backwards. Determine the micro-physics (Stage 1), find the required statistical spectrum to get the 5:1 ratio (Stage 2), and then check if the Graph produces that spectrum (Stage 3).

---

## 2. Stage 1: The Induction Test (Micro-Physics)

**Goal:** Determine the **Selection Function** $P(l=1 \mid R, \sigma_\phi)$.
*   Does a scalar lump of radius $R$ in a phase wind of strength $\sigma_\phi$ turn into a Baryon ($l=1$) or remain Dark Matter ($l=0$)?

### 2.0 Theoretical Requirement: The Energy Barrier
The transition from a Scalar Breather ($l=0$) to a Vector Rotor ($l=1$) is a **topological phase transition**. It is not a smooth rotation; it requires the field amplitude $|\Psi|$ to pass through zero to create the defect core.
*   **The Barrier:** To create a defect, the system must "punch a hole" in the condensate. This costs potential energy proportional to the well depth $\Delta V \sim \beta_{\rm pot}^2 / (4\gamma)$.
*   **The Driver:** The "Phase Wind" (gradient noise) provides the kinetic energy density $\mathcal{E}_K \sim \frac{1}{2} \kappa w^2 (\nabla \phi)^2$.
*   **The Collapse Condition:** A transition is only possible if the local kinetic stress exceeds the potential barrier:
    $$ \frac{1}{2} \kappa w^2 (\nabla \phi)^2 > V_{\text{barrier}} $$
    This implies a critical "Phase Reynolds Number" or noise threshold below which the lump simply stretches but never snaps into a vortex.

### 2.1 Protocol
1.  **Setup:** Initialize a stable Scalar Breather ($l=0$) of radius $R$ in a 2D/3D box.
2.  **Forcing:** Apply a "Phase Wind" $\phi_{ext}(x) = \vec{k} \cdot \vec{x} + \eta(x)$ at the boundaries, or add a driving term to the phase evolution.
    *   $\vec{k}$: Coherent flow (bulk rotation).
    *   $\eta$: Random noise (variance $\sigma_\phi$).
3.  **Evolution:** Evolve SGLE.
4.  **Detection:** Measure the winding number $n$ of the lump after time $T$.
    *   If $n \neq 0$: **Transition to Baryon**.
    *   If $n = 0$: **Remains Dark Matter**.

### 2.2 Output: The Phase Diagram
Construct a 2D map in $(R, \sigma_\phi)$ space:
*   **Region A (Baryonic):** Large $R$, Strong Wind $\to$ Defect Formation.
*   **Region B (Dark Matter):** Small $R$, Weak Wind $\to$ Phase Slip / Breather.

---

## 3. Stage 2: The Rogue Wave Fit (Statistics)

**Goal:** Find the **Power Spectrum** $P(k)$ that yields $\Omega_{DM} / \Omega_B \approx 5$.

### 3.1 Protocol
1.  **Assumption:** The primordial field follows "Rogue Wave" statistics (e.g., Gaussian Random Field or Log-Normal).
2.  **Sampling:**
    *   Generate a synthetic ensemble of lumps based on a trial Power Spectrum $P(k) \propto k^{-n_s}$.
    *   This gives a distribution of Lump Sizes $P(R)$ and Local Phase Gradients $P(\sigma_\phi)$.
3.  **Convolution:** Apply the **Selection Function** from Stage 1.
    *   $N_{B} = \int P(R, \sigma_\phi) \cdot P(l=1 \mid R, \sigma_\phi) \, dR \, d\sigma_\phi$.
    *   $N_{DM} = \int P(R, \sigma_\phi) \cdot P(l=0 \mid R, \sigma_\phi) \, dR \, d\sigma_\phi$.
4.  **Optimization:** Vary the spectral index $n_s$ until $E_{DM} / E_B \approx 5.3$.

### 3.2 Output: The Required Spectrum
*   "To explain the universe, the primordial noise must have had a spectral index of $n_s = X$."

---

## 4. Stage 3: The Graph Origin (Topology)

**Goal:** Determine if the **Infinite-Clique Graph** naturally produces the required spectrum $n_s$.

### 4.1 Protocol
1.  **Graph Synthesis:** Generate a scale-free graph with exponent $\gamma$ (Barabasi-Albert).
2.  **Field Mapping:** Compute the Laplacian Spectrum or the local stiffness statistics of this graph.
3.  **Comparison:** Does the graph's spectral index match the $n_s$ found in Stage 2?

### 4.2 Synthesis Result
*   **If Match:** The ICG naturally explains the Dark Matter ratio.
*   **If Mismatch:** The model requires modification (e.g., different graph topology).

---

## 5. Implementation Specs

### 5.1 Code Structure
*   `stage1_induction.py`: Single-soliton SGLE solver. (High Resolution, small domain).
*   `stage2_statistics.py`: Monte Carlo integrator. (No PDE solving, just probability distributions).
*   `stage3_graph.py`: Graph generation and spectral analysis. (NetworkX / SciPy).

### 5.2 Key Metrics
*   **Yield Ratio:** $R = E_{DM} / E_B$.
*   **Selection Threshold:** The critical "Reynolds Number" for phase winding $\text{Re}_\phi \sim R \cdot \nabla \phi$.

---

## 6. Output Interface

### 6.1 Stage 1 JSON (Selection Map)
```json
{
  "stage": 1,
  "map_resolution": [50, 50],
  "param_R": [0.1, 5.0],
  "param_sigma": [0.0, 2.0],
  "transition_matrix": "[[0,0,1...], ...]"
}
```

### 6.2 Stage 2 JSON (The Fit)
```json
{
  "stage": 2,
  "best_fit_ns": 2.4,
  "yield_ratio": 5.28,
  "confidence": 0.95
}
```

---

## 7. Synthesis Role

This 3-stage pipeline provides a **rigorous, falsifiable path** to cosmology. It separates the **Universal Physics** (how solitons react to noise) from the **Contingent History** (what noise distribution existed).
