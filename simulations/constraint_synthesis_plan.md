# Constraint Synthesis Plan: Combining Simulation Results

This document defines the unified framework for combining results from the four simulation tasks into constraints on the fundamental parameters of the Infinite-Clique Graph (ICG) model.

---

## 1. The Parameter Space

### 1.1 Canonical Microphysics Parameters

| Symbol | Name | Role |
|--------|------|------|
| $\alpha_{\rm grad}$ | Gradient Stiffness | Energy cost of spatial variation in $\Psi$ |
| $\beta_{\rm pot}$ | Potential Curvature | Depth of symmetry-breaking well |
| $\gamma$ | Quartic Coupling | Stabilizes the potential at large amplitude |
| $\kappa$ | Phase Stiffness | Sets the phase-wave speed (emergent $c$) |

### 1.2 Derived Scales

| Symbol | Formula | Physical Meaning |
|--------|---------|------------------|
| $w_*$ | $\sqrt{\beta_{\rm pot} / 2\gamma}$ | Vacuum amplitude |
| $m_\xi$ | $\sqrt{2\beta_{\rm pot}}$ | Amplitude gap (mass scale) |
| $\ell$ | $\sqrt{\alpha_{\rm grad} / 2\beta_{\rm pot}}$ | Correlation length / core width |
| $c_s$ | $\sqrt{\kappa} w_*$ | Phase speed (emergent "speed of light") |
| $T_{\rm string}$ | $\pi \ell^2 \beta_{\rm pot} w_*^2$ | Vortex string tension |

### 1.3 Degrees of Freedom

The 4 canonical parameters have **4 degrees of freedom**, but physical observables depend on:

- **3 dimensionless ratios** (e.g., $\lambda = \gamma/\alpha_{\rm grad}$, $\ell/a$, $\kappa/\alpha_{\rm grad}$)
- **1 overall scale** (set by a physical anchor)

We can parameterize as:

| Reduced Parameter | Definition | Range |
|-------------------|------------|-------|
| $\lambda$ | $\gamma / \alpha_{\rm grad}$ | $[0.1, 10]$ |
| $\xi$ | $\ell / \ell_{\rm QCD}$ | $[0.5, 2]$ (near unity if $\ell \sim 1$ fm) |
| $\eta$ | $\kappa / \alpha_{\rm grad}$ | $[0.1, 10]$ |
| $\ell_{\rm QCD}$ | **Anchor** | Fixed at $1$ fm |

---

## 2. Physical Anchors

To convert dimensionless simulation outputs to physical predictions, we need **anchor quantities** that tie the model to measured values.

### 2.1 Primary Anchor: QCD/Nuclear Scale

| Quantity | Value | Sets |
|----------|-------|------|
| $\ell = 1$ fm | $10^{-15}$ m | The correlation length / core width |

This choice is motivated by the expectation that the amplitude gap is near the QCD scale ($\sim 200$ MeV $\sim 1$ fm$^{-1}$ in natural units).

From this anchor:
$$
\alpha_{\rm grad} = 2 \beta_{\rm pot} \cdot \ell^2 = 2 \beta_{\rm pot} \cdot (1 \text{ fm})^2
$$

### 2.2 Secondary Anchors (Consistency Checks)

| Quantity | Value | Constrains |
|----------|-------|------------|
| $G$ (Newton's constant) | $6.674 \times 10^{-11}$ m³/kg/s² | $C_{\rm eff} / (\kappa w_*^2)$ |
| $m_p$ (proton mass) | $0.938$ GeV | $T_{\rm string} \cdot L_{\rm junction}$ |
| $c$ (speed of light) | $3 \times 10^8$ m/s | $c_s = \sqrt{\kappa} w_*$ |

### 2.3 Anchor Consistency Condition

If the model is correct, these anchors must be mutually consistent:
$$
G \cdot m_p^2 / (\hbar c) \sim 10^{-38} \quad \text{(observed)}
$$
The model must reproduce this hierarchy from the ratios $(\lambda, \xi, \eta)$.

---

## 3. Constraint Functions (From Each Simulation)

Each simulation produces a constraint of the form: $f_i(\theta) \in [\text{target}_i^{\rm lo}, \text{target}_i^{\rm hi}]$.

### 3.1 Geometry Simulation

**Output:** Spectral dimension $d_s(\alpha, T, \beta)$.

**Constraint Function:**
$$
f_1(\theta) = d_s \quad \Rightarrow \quad d_s \in [2.9, 3.1]
$$

**What it constrains:** The vacuum functional coefficients $(\alpha, T, \beta)$ in $F[J]$. This is a **model validation**—if $d_s \neq 3$, the entire framework fails.

**Parameter mapping:** The coefficients $(\alpha, T, \beta)$ are related to $(\alpha_{\rm grad}, \kappa, T_{\rm eff})$ via the vacuum marginalization derivation (see ICG SI). For now, treat this simulation as a binary pass/fail.

### 3.2 Lattice Simulation

**Output:** Loop-to-junction mass ratio $R_{\rm length}(\lambda, \text{quench})$.

**Constraint Function:**
$$
f_2(\lambda, q) = R_{\rm length} \quad \Rightarrow \quad R_{\rm length} \in [0.75, 0.90]
$$
where $q$ parameterizes the quench schedule (a nuisance parameter).

**What it constrains:** The dimensionless ratio $\lambda = \gamma / \alpha_{\rm grad}$.

**Target:** $R_{\rm length} \approx 0.84$ corresponds to $\Omega_{\rm DM} / \Omega_b \approx 5.3$.

**Nuisance handling:** Run multiple quench schedules. Report the **envelope** of $\lambda$ values that achieve $R_{\rm length} \in [0.75, 0.90]$ for *any* reasonable quench.

### 3.3 Knot Simulation

**Output:** Metastable radius $R_g(k_B)$ and self-interaction cross-section $\sigma/m$.

**Constraint Functions:**
$$
f_3^{(a)}(k_B) = R_g \quad \Rightarrow \quad R_g > 0 \text{ (stability)}
$$
$$
f_3^{(b)}(k_B, R_g) = \frac{\sigma}{m} \quad \Rightarrow \quad \frac{\sigma}{m} < 1 \text{ cm}^2/\text{g (SIDM bound)}
$$

**What it constrains:** The bending modulus $k_B \propto \alpha_{\rm grad}^2 / \beta_{\rm pot}$, which maps to a threshold on the ratio $\alpha_{\rm grad} / \beta_{\rm pot}$.

**Combined:** Define $k_B^{\rm crit}$ as the minimum bending for which both conditions hold.

### 3.4 Poisson Simulation

**Output:** Effective gravitational coupling $C_{\rm eff}(\kappa_0, \epsilon)$.

**Constraint Function:**
$$
f_4(\kappa, \alpha_{\rm grad}) = C_{\rm eff} \quad \Rightarrow \quad C_{\rm eff} = G \cdot c_s^2 = G \cdot \kappa w_*^2
$$

**What it constrains:** The combination $\kappa \cdot w_*^2 = \kappa \cdot \beta_{\rm pot} / (2\gamma)$.

**Calibration:** This is the strongest quantitative constraint. Given $G$ is measured to high precision, the model must hit this value.

---

## 4. Combination Procedure

### 4.1 Phase 1: Constraint Intersection (Current Approach)

**Method:** Grid search + set intersection.

**Steps:**

1. **Define a grid** over the reduced parameter space:
   - $\lambda \in [0.1, 10]$ (log-spaced, 20 points)
   - $\xi \in [0.5, 2]$ (linear, 10 points)
   - $\eta \in [0.1, 10]$ (log-spaced, 20 points)

2. **For each grid point**, check if all constraints are satisfied:
   - $f_1 \in [2.9, 3.1]$ (Geometry) — or skip if treating as pass/fail
   - $f_2 \in [0.75, 0.90]$ (Lattice)
   - $f_3^{(a)} > 0$ and $f_3^{(b)} < 1$ (Knot)
   - $|f_4 - G \cdot c_s^2| / (G \cdot c_s^2) < 0.1$ (Poisson)

3. **Mark viable points.** The union of viable points is the **viable region** $V$.

4. **Report:**
   - Is $V$ non-empty? (Model viability)
   - What are the bounds on each reduced parameter?
   - Are the bounds tight or loose?

### 4.2 Phase 2: Bayesian Upgrade (Future)

**Method:** Surrogate model + MCMC.

**Steps:**

1. **Train a Gaussian Process (GP)** on simulation outputs:
   - Input: $(\lambda, \xi, \eta)$
   - Output: $(R_{\rm length}, R_g, \sigma/m, C_{\rm eff})$
   - Use $\sim 100$ simulation runs to train.

2. **Define likelihoods:**
   $$
   \mathcal{L}_{\rm Lattice} = \exp\left( -\frac{(R_{\rm length} - 0.84)^2}{2 \sigma_R^2} \right)
   $$
   $$
   \mathcal{L}_{\rm Knot} = \text{sigmoid}\left( \frac{R_g}{\delta_R} \right) \cdot \text{sigmoid}\left( \frac{1 - \sigma/m}{\delta_\sigma} \right)
   $$
   $$
   \mathcal{L}_{\rm Poisson} = \exp\left( -\frac{(C_{\rm eff} - G c_s^2)^2}{2 \sigma_C^2} \right)
   $$

3. **Sample the posterior:**
   $$
   P(\theta | \text{data}) \propto \mathcal{L}_{\rm Lattice} \cdot \mathcal{L}_{\rm Knot} \cdot \mathcal{L}_{\rm Poisson} \cdot P(\theta)
   $$
   Use `emcee` or `dynesty` for sampling.

4. **Report:** Posterior distributions, credible intervals, corner plots.

---

## 5. Success Criteria

### 5.1 Model Viability

The model passes if the viable region $V$ is **non-empty** and contains parameter values that are:

- **Natural:** Dimensionless ratios are $O(1)$, not fine-tuned.
- **Consistent:** The anchors $(G, m_p, c)$ are mutually compatible.

### 5.2 Predictive Power

The model is interesting if it makes **testable predictions** beyond the anchors, e.g.:

- DM self-interaction cross-section (within SIDM bounds but potentially detectable).
- Gravitational wave signatures from defect annihilation.
- Deviations from Newtonian gravity at short range (fifth-force experiments).

### 5.3 Failure Modes

| Outcome | Interpretation |
|---------|----------------|
| $V = \emptyset$ | Model is ruled out by current constraints |
| $V$ requires $\lambda \gg 1$ or $\ll 1$ | Fine-tuning problem |
| $C_{\rm eff}$ cannot match $G$ | Gravity sector is broken |
| $R_g = 0$ always | DM candidate is unstable |

---

## 6. Data Products and Deliverables

### 6.1 From Each Simulation

| Simulation | Primary Output File | Key Fields |
|------------|---------------------|------------|
| Geometry | `geometry_results.json` | `final_d_s`, `converged` |
| Lattice | `lattice_results.json` | `ratio_R_length`, `canonical_params.lambda` |
| Knot | `knot_results.json` | `radius_of_gyration`, `sigma_over_m`, `stable` |
| Poisson | `poisson_results.json` | `C_eff`, `fit_R_squared` |

### 6.2 Synthesis Output

| File | Contents |
|------|----------|
| `viable_region.json` | List of $(\lambda, \xi, \eta)$ points that pass all constraints |
| `constraint_summary.md` | Human-readable summary of bounds |
| `corner_plot.png` | (Phase 2) Posterior distribution visualization |

---

## 7. Workflow Summary

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         SIMULATION PHASE                                │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│   ┌──────────────┐   ┌──────────────┐   ┌──────────────┐   ┌──────────┐│
│   │   Geometry   │   │   Lattice    │   │    Knot      │   │ Poisson  ││
│   │   (d_s → 3)  │   │  (R_length)  │   │   (R_g, σ/m) │   │  (C_eff) ││
│   └──────┬───────┘   └──────┬───────┘   └──────┬───────┘   └────┬─────┘│
│          │                  │                  │                │      │
│          ▼                  ▼                  ▼                ▼      │
│   geometry_results   lattice_results   knot_results   poisson_results  │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                         SYNTHESIS PHASE                                 │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│   1. Load all *_results.json files                                      │
│   2. Map simulation params → canonical params (using Section 0 of each) │
│   3. Evaluate constraint functions f₁, f₂, f₃, f₄                       │
│   4. Identify viable region V = {θ : all constraints satisfied}         │
│   5. Report bounds on (λ, ξ, η) and physical predictions                │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                         OUTPUT                                          │
├─────────────────────────────────────────────────────────────────────────┤
│   viable_region.json   constraint_summary.md   (optional: corner_plot)  │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 8. Instructions for Subagents

Each simulation subagent should:

1. **Read its own simulation plan** for task details.
2. **Implement the simulation** and produce JSON output matching the schema.
3. **Include `canonical_params`** in the output so results can be mapped to the shared space.
4. **Do not assume knowledge of other simulations**—each task is independent.

The **synthesis step** (this document) will be executed after all simulations complete.
