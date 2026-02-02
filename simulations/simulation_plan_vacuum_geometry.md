# Simulation Task: Graph Geometrogenesis & Spectral Dimension

---

## 0. Parameter Interface (Canonical Microphysics)

This simulation is part of a unified constraint program. It operates at a **foundational level**: it validates the emergence of a 3D lattice substrate on which the other simulations run.

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
| $\ell$ | $\sqrt{\alpha_{\rm grad} / 2\beta_{\rm pot}}$ | Correlation length / lattice spacing proxy. |
| $c_s$ | $\sqrt{\kappa} w_*$ | Phase wave speed. |

### 0.3 Mapping to This Simulation's Parameters

This simulation models the **vacuum marginalization** process: integrating out fast phase fluctuations to derive an effective free energy for the slow graph weights $J_{ij}$.

| Simulation Param | Origin | Interpretation |
|------------------|--------|----------------|
| $\alpha$ (Log-Det weight) | $\propto T_{\rm eff}$ (effective temperature of phase fluctuations) | Strength of the "expander" drive (favors connectivity). |
| $T$ (Entropy weight) | Related to the entropic drive on links | Favors democratic weight spreading. |
| $\beta$ (Dimension bias) | Phenomenological | Strength of the penalty for $d_s \neq 3$. |

**Connection to other simulations:**
- If this simulation succeeds ($d_s \to 3$, local geometry emerges), the resulting graph *is* the lattice on which the Lattice, Knot, and Poisson simulations are defined.
- The emergent lattice spacing $a$ should be identified with the correlation length $\ell$ from the field theory.

---

## 1. Theoretical Context: Spacetime from Graphs

This project explores **Geometrogenesis**: the hypothesis that continuous 3D spacetime emerges from a disordered, high-dimensional graph structure.

We posit a "Vacuum Free Energy" $F[J]$ that governs the evolution of the graph's link weights $J_{ij}$. This functional contains competing terms:

1.  **Entropic Drive:** Favors random, "small-world" connectivity (infinite dimension).
2.  **Dimensional Bias:** Explicitly penalizes configurations where the "Spectral Dimension" $d_s \neq 3$.

**Goal:** Demonstrate that minimizing this functional naturally drives a random graph to "condense" into a structure that looks like a 3D lattice.

---

## 2. The Mathematical Model

We operate on a graph of $N$ nodes. The state is defined by the symmetric adjacency matrix $J_{ij} \ge 0$.

### 2.1 The Vacuum Functional

\[
F[J] = \alpha E_{\text{phase}}[J] - T S[J] + \beta (d_s[J] - d^*)^2
\]

-   **Phase Energy:** $E_{\text{phase}} = \text{Tr}(\ln (L + \epsilon I))$. $\epsilon=10^{-6}$ is a regulator to handle the zero mode.
-   **Entropy:** $S[J] = - \sum_{i,j} P_{ij} \ln P_{ij}$, where $P_{ij} = J_{ij} / \sum J$.
-   **Dimension Bias:** Penalizes deviation from target dimension $d^* = 3$.

### 2.2 The Graph Laplacian

\[ L_{ij} = \begin{cases} -J_{ij} & i \neq j \\ \sum_k J_{ik} & i = j \end{cases} \]

### 2.3 Measuring Spectral Dimension

\[ d_s(\tau) = -2 \frac{\partial \ln P(\tau)}{\partial \ln \tau} \]
where $P(\tau) = \frac{1}{N} \text{Tr}(e^{-L\tau}) = \frac{1}{N} \sum_{k} e^{-\lambda_k \tau}$.

-   $\lambda_k$ are the eigenvalues of the Laplacian $L$.
-   We target a specific scale $\tau_* \approx 1.0$ (or average over a "diffusive window").

---

## 3. Simulation Protocol

### 3.1 Initialization

Start with a **Random Geometric Graph** or Erdős-Rényi graph to ensure initial connectivity but high dimension.

-   $N = 100$ nodes.
-   Random weights $J_{ij} \in [0, 1]$ for connected pairs.

### 3.2 Optimization (Gradient Descent)

**Update Loop:**

1.  Compute Laplacian $L$ and its eigenvalues/vectors via `numpy.linalg.eigh(L)`.
2.  Compute current $d_s(\tau_*)$.
3.  Compute gradients:
    -   **Analytic** for Entropy: $\partial S / \partial J_{ij} = -\ln(J_{ij}/Z) - 1$.
    -   **Analytic** for Log-Det: $\partial E_{\rm phase} / \partial J_{ij} = (L^{-1})_{ii} + (L^{-1})_{jj} - 2(L^{-1})_{ij}$ (Effective Resistance).
    -   **Analytic** for Dimension Bias: Use eigenvalue perturbation theory (see Section 5.1).
4.  Update: $J_{ij} \leftarrow J_{ij} - \eta \nabla F$.
5.  **Project:** Enforce $J_{ij} \ge 0$ (set negatives to 0).
6.  **Normalize:** Rescale so $\sum J = \text{const}$ (prevents runaway weights).

### 3.3 Analysis Observables

1.  **Dimension Flow:** Plot $d_s$ vs. Step.
2.  **Embedding (Isomap):**
    -   Use `sklearn.manifold.Isomap` with `n_neighbors=10`.
    -   Embed the final graph into 3D.
    -   Calculate **Residual Variance** to see if it fits well in 3D compared to higher dimensions.
3.  **Visualization:** Plot the 3D embedding. Does it look like a blob, a line, or a grid?
4.  **Regularity:** Compute the degree distribution. A 3D lattice should have near-constant degree $\approx 6$.

---

## 4. Parameter Study

-   **Bias Strength:** $\beta \in [0, 10, 100]$.
-   **Target Dimension:** Try driving to $d^*=1$ (line), $d^*=2$ (sheet), and $d^*=3$.
-   **Initial Connectivity:** Sparse (Erdős-Rényi $p=0.1$) vs Dense ($p=0.5$).

---

## 5. Implementation Specs

-   **Language:** Python.
-   **Libraries:** `numpy`, `scipy`, `sklearn` (for Isomap).
-   **Computational Cost:** Eigenvalue decomposition is $O(N^3)$. Keep $N \le 200$ for tractability.

### 5.1 Analytic Eigenvalue Gradient (Recommended)

To avoid slow finite differences:
\[ \frac{\partial d_s}{\partial J_{uv}} = \sum_k \frac{\partial d_s}{\partial \lambda_k} \frac{\partial \lambda_k}{\partial J_{uv}} \]

-   $\frac{\partial d_s}{\partial \lambda_k} = \frac{2\tau e^{-\lambda_k \tau}}{P(\tau)} \left( d_s(\tau) / 2 - 1 \right) / N$ (chain rule on the log-derivative).
-   $\frac{\partial \lambda_k}{\partial J_{uv}} = (v_k^{(u)} - v_k^{(v)})^2$ where $v_k$ is the $k$-th eigenvector.

---

## 6. Output Interface (Structured Results)

All runs must output a **JSON file** with the following schema:

```json
{
  "simulation": "VacuumGeometry",
  "version": "1.0",
  "params": {
    "N": 100,
    "alpha": 1.0,
    "T": 0.1,
    "beta": 10.0,
    "d_star": 3,
    "tau_star": 1.0,
    "initial_graph": "ErdosRenyi_p0.3",
    "steps": 1000,
    "eta": 0.01
  },
  "outputs": {
    "converged": true,
    "final_d_s": 3.02,
    "d_s_std": 0.15,
    "isomap_residual_3D": 0.05,
    "isomap_residual_2D": 0.45,
    "mean_degree": 5.8,
    "std_degree": 1.2,
    "fiedler_eigenvalue": 0.12
  },
  "embedding": {
    "coordinates_file": "embedding_3D.npy"
  }
}
```

### 6.1 Key Outputs Explained

| Field | Meaning | Physical Interpretation |
|-------|---------|-------------------------|
| `final_d_s` | Spectral dimension at convergence | Should be $\approx 3$ if geometrogenesis succeeds. |
| `isomap_residual_3D` | Reconstruction error for 3D embedding | Low residual = graph is intrinsically 3D. |
| `mean_degree` | Average node degree | For a 3D cubic lattice, expect $\approx 6$. |
| `fiedler_eigenvalue` | Second smallest eigenvalue of $L$ | Measures connectivity; $>0$ means connected. |

---

## 7. Synthesis Role (How This Constrains the Model)

This simulation is the **foundation check**: it validates that the vacuum marginalization process produces a 3D spatial substrate.

### 7.1 Target Observable

The primary target is **Spectral Dimension = 3**:
\[
d_s(\tau_*) = 3.0 \pm 0.1
\]
and **Local Geometry**: the graph can be embedded in $\mathbb{R}^3$ with low distortion.

### 7.2 Explicit Constraint Function

**For synthesis, report the following constraint:**

$$
f_{\rm Geometry}(\alpha, T, \beta) = d_s(\tau_*)
$$

**Target band:**
$$
d_s \in [2.9, 3.1]
$$

**This is a model validation constraint (binary pass/fail):**
- If $d_s \in [2.9, 3.1]$: **PASS** — The vacuum functional produces 3D geometry. Proceed.
- If $d_s \notin [2.9, 3.1]$: **FAIL** — The model does not produce 3D spacetime. Other simulations are invalid.

**Parameter mapping:** The functional coefficients $(\alpha, T, \beta)$ are phenomenological in this simulation. They relate to the canonical microphysics via the vacuum marginalization derivation (integrating out fast phase modes). For synthesis purposes:
- $\alpha \propto T_{\rm eff}$ (effective temperature of phase fluctuations)
- $\beta$ (dimension bias) is a free knob that must be "large enough" to enforce $d_s = 3$

**What this constrains:** This simulation validates the *existence* of a 3D substrate, not the numerical values of $(\alpha_{\rm grad}, \beta_{\rm pot}, \gamma, \kappa)$. It is a prerequisite for the other simulations.

### 7.3 Interpretation

- **If $d_s \to 3$ and `isomap_residual_3D` is low:** The vacuum drive successfully produces 3D geometry. Proceed with other simulations on this substrate.
- **If $d_s \neq 3$:** The model fails to produce 3D spacetime. Either adjust $\beta$ (bias strength) or conclude the vacuum functional is incomplete.
- **If graph is not regular (high `std_degree`):** The emergent geometry is "glassy" rather than crystalline. This may still be acceptable if it's locally 3D.

### 7.4 Interface with Other Simulations

- **Lattice Simulation:** The output graph (or a 3D cubic lattice approximating it) is the substrate for the field simulation.
- **Knot Simulation:** The effective string model assumes a 3D ambient space. If $d_s \neq 3$, the model may need modification.
- **Poisson Simulation:** The discrete Laplacian is defined on the 3D lattice. If the lattice doesn't exist, the simulation is moot.

---

## 8. Failure Modes & Mitigation

### 8.1 Fragmentation

-   **Issue:** The optimizer sets most $J_{ij} \to 0$, breaking the graph into disconnected islands.
-   **Mitigation:** Add a penalty for the **Fiedler Eigenvalue** (algebraic connectivity): add $+ \gamma_{\rm conn} (\lambda_2 - \lambda_{\rm target})^2$ to $F$ to force $\lambda_2 > 0$. Or, start with a Minimum Spanning Tree and only relax the *extra* edges.

### 8.2 Topological Trap (Star Graph)

-   **Issue:** The graph optimizes to a "Star Graph" (Hub-and-spoke) which has low diameter but isn't a lattice.
-   **Mitigation:** Enforce a **Maximum Degree Constraint** (e.g., $\deg(v) \le 10$). Or add a repulsion term for high-degree nodes: $-\sum_i (\deg_i - 6)^2$.

### 8.3 Oscillation

-   **Issue:** Gradient descent oscillates around the target $d_s=3$ without settling.
-   **Mitigation:** Use **Momentum** (Adam optimizer) or a decaying learning rate schedule. Or use a soft constraint (anneal $\beta$ from 0 to target).

### 8.4 Local Minima

-   **Issue:** The graph gets stuck in a local minimum (e.g., $d_s = 2.5$) and won't move to $d_s = 3$.
-   **Mitigation:** Use **Simulated Annealing**: add noise to the updates early, then cool down. Or, use multiple random restarts and select the best final state.

### 8.5 Eigenvalue Degeneracy

-   **Issue:** Degenerate eigenvalues cause numerical instability in the gradient formula.
-   **Mitigation:** Add a small random perturbation to $L$ before computing eigenvalues, or use a regularized pseudo-inverse for nearly-degenerate subspaces.
