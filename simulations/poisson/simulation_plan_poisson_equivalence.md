# Simulation Task: Calibrating the "Noise Field" Gravity Law

---

## 0. Parameter Interface (Canonical Microphysics)

This simulation is part of a unified constraint program. All simulations share a common set of **Canonical Microphysics Parameters** derived from the underlying complex scalar field theory.

### 0.1 Canonical Parameters (Shared Across All Simulations)

| Symbol | Name | Description |
|--------|------|-------------|
| $\alpha_{\rm grad}$ | Gradient Stiffness | Coefficient of $\|\nabla\Psi\|^2$ in the free energy. |
| $\beta_{\rm pot}$ | Potential Curvature | Coefficient of $-w^2$ in the Mexican-hat potential. |
| $\gamma$ | Quartic Coupling | Coefficient of $w^4$. |
| $\kappa$ | Phase Stiffness | Coefficient of $w^2\|\nabla\phi\|^2$. Sets the phase-wave speed. |

### 0.2 Derived Physical Scales

| Symbol | Formula | Interpretation |
|--------|---------|----------------|
| $w_*$ | $\sqrt{\beta_{\rm pot} / 2\gamma}$ | Vacuum amplitude. |
| $m_\xi$ | $\sqrt{2\beta_{\rm pot}}$ | Amplitude gap / mass. |
| $\ell$ | $\sqrt{\alpha_{\rm grad} / 2\beta_{\rm pot}}$ | Correlation length / core width. |
| $c_s$ | $\sqrt{\kappa} w_*$ | Phase wave speed (effective "speed of light"). |

### 0.3 Mapping to This Simulation's Parameters

This simulation models the **phase sector** as a diffusion/conductivity problem. A matter defect (amplitude perturbation) changes the local stiffness, sourcing a $1/r$ noise-variance field.

| Simulation Param | Formula from Canonical | Interpretation |
|------------------|------------------------|----------------|
| $\kappa_0$ (vacuum conductivity) | $\kappa_0 = \kappa w_*^2$ | Phase stiffness in the uniform vacuum. |
| $\epsilon$ (coupling strength) | $\epsilon \propto M / \kappa_0$ | How much a mass $M$ perturbs the stiffness. |

**In simulation units** where $\kappa_0 = 1$:
- $\epsilon$ directly represents the integrated perturbation strength.
- The output $C_{\rm eff}$ relates to the Newtonian constant $G$ via:
\[
G_{\rm eff} = \frac{C_{\rm eff}}{c_s^2} = \frac{C_{\rm eff}}{\kappa w_*^2}
\]

---

## 1. Theoretical Context: Emergent Gravity

In this framework, the gravitational potential $\Phi$ is not fundamental. Instead, it is identified with the local variance of quantum phase fluctuations, denoted $\tau^2$.
\[ \Phi(x) \propto \tau^2(x) \equiv \langle |\nabla \phi|^2 \rangle_{\text{local window}} \]

Matter (defects in the field amplitude) couples to the phase field by altering the local **stiffness** (conductivity) of the vacuum. This perturbation changes the propagation of phase fluctuations, creating a "shadow" of increased noise variance around the object.

**Goal:** Verify that for a localized matter source, the induced noise profile $\delta \tau^2(r)$ falls off as $1/r$, effectively recovering the Poisson equation of Newtonian gravity:
\[ \nabla^2 (\delta \tau^2) = -C_{\text{eff}} \cdot \rho_m \]

---

## 2. The Mathematical Model

We treat the vacuum as a resistor network with variable conductivity $\kappa(\mathbf{x})$.

### 2.1 The Discrete Laplacian

On a 3D lattice, the operator $L = -\nabla \cdot (\kappa \nabla)$ acts on a vector $\phi$ as:
\[ (L\phi)_{\mathbf{x}} = \sum_{\mu=1}^3 \left[ \kappa_{\mathbf{x}, \mathbf{x}+\hat{\mu}} (\phi_{\mathbf{x}} - \phi_{\mathbf{x}+\hat{\mu}}) + \kappa_{\mathbf{x}, \mathbf{x}-\hat{\mu}} (\phi_{\mathbf{x}} - \phi_{\mathbf{x}-\hat{\mu}}) \right] \]

-   **Conductivity:** The edge stiffness $\kappa_{\mathbf{x}, \mathbf{y}}$ is the harmonic mean of the site stiffnesses:
    \[ \kappa_{\mathbf{x}, \mathbf{y}} = \frac{2 \kappa(\mathbf{x}) \kappa(\mathbf{y})}{\kappa(\mathbf{x}) + \kappa(\mathbf{y})} \]

### 2.2 The Perturbation (Matter Source)

-   **Vacuum:** $\kappa(\mathbf{x}) = 1.0$ everywhere.
-   **Matter:** A spherical defect of radius $R_{\rm source}=3$ at center.
    \[ \kappa(\mathbf{x}) = 1.0 + \epsilon \] for $|\mathbf{x} - \mathbf{x}_c| < R_{\rm source}$.
    $\epsilon$: Coupling strength (e.g., 0.5).

### 2.3 The Observable: Noise Variance

The noise variance at site $\mathbf{x}$ is the diagonal of the inverse operator:
\[ \tau^2(\mathbf{x}) = (L^{-1})_{\mathbf{x}\mathbf{x}} \]
Signal: $\delta \tau^2(\mathbf{x}) = \tau^2(\mathbf{x}) - \tau^2_{\text{vacuum}}$.

---

## 3. Simulation Protocol

### 3.1 Matrix Construction

Build $L$ as a `scipy.sparse.csr_matrix`.

-   Size: $N^3 \times N^3$.
-   **Boundary Conditions:** Use **Dirichlet** boundaries ($\phi=0$ at edges). Do **NOT** use Periodic boundaries, as they enforce charge neutrality.
    -   *Correction:* Place the source far from boundaries ($N \ge 30$).

### 3.2 Solving (Stochastic Estimation)

For $N > 20$, exact inversion is too slow. Use the **Hutchinson Trace Estimator** to find the diagonal.

**Algorithm:**

1.  Initialize accumulator vector $D = \mathbf{0}$.
2.  Loop $K=100$ times:
    a. Generate random probe vector $\mathbf{v}$ with entries $\pm 1$ (Rademacher distribution).
    b. Solve the linear system $L \mathbf{u} = \mathbf{v}$ for $\mathbf{u}$.
        -   Use `scipy.sparse.linalg.cg` (Conjugate Gradient).
        -   Preconditioner: `scipy.sparse.linalg.spilu` (Incomplete LU).
    c. Update accumulator: $D_i \leftarrow D_i + v_i u_i$.
3.  Average: $\tau^2_i \approx D_i / K$.

### 3.3 Analysis

-   Compute radial profile $\delta \tau^2(r)$ by binning lattice sites based on distance from center.
-   Fit model: $f(r) = \frac{A}{r} + B$ (constant offset for finite-size effects).
-   Extract $A$ and check fit residuals ($R^2 > 0.95$).

---

## 4. Parameter Study

1.  **Lattice Size:** $N \in [20, 40]$ (Direct Solver) and $N=60$ (Stochastic).
2.  **Coupling:** $\epsilon \in [0.1, 1.0, 10.0]$. Verify linearity $A \propto \epsilon$.
3.  **Source Size:** $R_{\rm source} \in [2, 5]$. Check that $A$ depends only on total $\int \epsilon \, dV$, not shape.

---

## 5. Implementation Specs

-   **Language:** Python.
-   **Solver:** `scipy.sparse.linalg.cg`.
-   **Fitting:** `scipy.optimize.curve_fit`.
-   **Validation:** For small $N$, compare Stochastic result vs `inv(L.todense()).diagonal()`.

---

## 6. Output Interface (Structured Results)

All runs must output a **JSON file** with the following schema:

```json
{
  "simulation": "PoissonEquivalence",
  "version": "1.0",
  "params": {
    "N": 40,
    "epsilon": 0.5,
    "R_source": 3,
    "K_samples": 100,
    "boundary": "Dirichlet"
  },
  "canonical_params": {
    "kappa_0": 1.0,
    "alpha_grad": 1.0,
    "beta_pot": 0.5
  },
  "outputs": {
    "fit_A": 0.0234,
    "fit_B": 0.0012,
    "fit_R_squared": 0.987,
    "C_eff": 0.0468,
    "linearity_check": true
  },
  "diagnostics": {
    "cg_iterations_avg": 45,
    "hutchinson_variance": 0.0003
  }
}
```

### 6.1 Key Outputs Explained

| Field | Meaning | Physical Interpretation |
|-------|---------|-------------------------|
| `fit_A` | Coefficient of $1/r$ in the fit | Strength of the monopole potential. |
| `C_eff` | $2 \times \text{fit\_A} / (\epsilon \cdot V_{\rm source})$ | The emergent "Newton constant" in lattice units. |
| `linearity_check` | Is $A \propto \epsilon$ across runs? | Validates the linear response regime (weak field). |
| `fit_R_squared` | Goodness of fit | If $<0.9$, the $1/r$ law may not hold. |

---

## 7. Synthesis Role (How This Constrains the Model)

This simulation constrains the **gravity sector**: does the emergent potential obey Newtonian $1/r$?

### 7.1 Target Observable

The primary target is **Newton's Constant** $G$:
\[
G = 6.674 \times 10^{-11} \text{ m}^3 \text{kg}^{-1} \text{s}^{-2}
\]
In our framework:
\[
G_{\rm eff} = \frac{C_{\rm eff}}{c_s^2}
\]
where $c_s = \sqrt{\kappa} w_*$ is the phase speed.

### 7.2 Explicit Constraint Function

**For synthesis, report the following constraint:**

$$
f_{\rm Poisson}(\kappa, \alpha_{\rm grad}, \beta_{\rm pot}) = C_{\rm eff}
$$

**Target value (calibration condition):**
$$
C_{\rm eff} = G \cdot c_s^2 = G \cdot \kappa \cdot w_*^2 = G \cdot \kappa \cdot \frac{\beta_{\rm pot}}{2\gamma}
$$

where:
- $\kappa_0 = \kappa w_*^2$ is the vacuum phase stiffness used in the simulation
- $C_{\rm eff} = 2 \cdot \text{fit\_A} / (\epsilon \cdot V_{\rm source})$ is extracted from the $1/r$ fit

**This is a calibration constraint, not a band:** It fixes one combination of $(\kappa, \beta_{\rm pot}, \gamma)$ given the measured $G$.

**Practical use:** Given the anchor $\ell = 1$ fm (which fixes $\alpha_{\rm grad}/\beta_{\rm pot}$), this constraint determines the ratio $\kappa/\gamma$.

### 7.3 Interpretation

- **If $1/r$ holds and $C_{\rm eff}$ is constant:** Gravity is Newtonian in the weak-field limit. Extract the relation $C_{\rm eff}(\alpha_{\rm grad}, \beta_{\rm pot}, \kappa)$.
- **If $1/r$ fails (e.g., Yukawa $e^{-mr}/r$):** The model predicts deviations from Newtonian gravity at some scale. Compare to fifth-force experiments.
- **Calibration:** The numerical value of $C_{\rm eff}$ sets the scale at which $\alpha_{\rm grad}, \kappa$ become physical. Use Solar System tests to fix one combination of parameters.

### 7.4 Interface with Other Simulations

- **Lattice Simulation:** The "matter source" in this simulation corresponds to a vortex defect from the Lattice simulation. The stiffness perturbation $\epsilon$ is related to the defect's energy density.
- **Knot Simulation:** If knots are DM, this simulation validates that they gravitate normally (via the same $C_{\rm eff}$).
- **Geometry Simulation:** If the graph doesn't produce a 3D lattice, this simulation's Laplacian assumptions are invalid.

---

## 8. Failure Modes & Mitigation

### 8.1 Noise Floor

-   **Issue:** For small coupling $\epsilon$, the signal $\delta \tau^2$ is smaller than the stochastic noise of the Hutchinson estimator.
-   **Mitigation:** Use **Control Variates**. Solve for *both* $L$ and $L_0$ (vacuum) with the *same* random vectors $\mathbf{v}$. The difference $\delta \tau^2 = \mathbf{v}^T(L^{-1} - L_0^{-1})\mathbf{v}$ has correlated noise that cancels.

### 8.2 Boundary Effects

-   **Issue:** Dirichlet boundaries create "image charges" that distort the potential slope.
-   **Mitigation:** Fit only the inner region ($R_{\rm source} < r < N/4$). Or, use a "Screened Poisson" solver (add a tiny mass term $m^2 \phi$) to make the Green's function decay exponentially, removing boundary sensitivity.

### 8.3 Anisotropy

-   **Issue:** On a coarse grid, the "sphere" looks like a cube/diamond, introducing quadrupole moments.
-   **Mitigation:** Use a "soft sphere" or Gaussian profile for $\kappa(x)$ instead of a hard step function. Look at the monopole term of the multipole expansion only.

### 8.4 Non-Linearity

-   **Issue:** For large $\epsilon$, the response may become nonlinear ($A \not\propto \epsilon$).
-   **Mitigation:** Verify linearity by plotting $A$ vs $\epsilon$. If nonlinear, restrict to small $\epsilon$ or model the nonlinear corrections.
