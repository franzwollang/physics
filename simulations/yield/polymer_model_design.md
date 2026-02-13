# Design Document: Polymer Physics Adaptation for Primordial Yield
**Status:** Draft / Iteration 1
**Context:** Enhancing Stage 1 of the Yield Simulation (Selection Function)

## 1. Core Concept: The "Spaghetti" Mapping

We model the formation of topological defects (leptons/baryons) inside a scalar lump as a polymer statistics problem. The "noodles" are the integral curves of the phase gradient field $\nabla\phi$ within the lump's volume.

### The Dictionary

| Physics Framework | Polymer Model | Physical Meaning |
| :--- | :--- | :--- |
| **Phase Gradient Path** | **Polymer Chain** | Integral curve of $\nabla\phi$. |
| **Phase Stiffness** ($\kappa w_*^2$) | **Persistence Length** ($\ell_p$) | Resistance to bending; stiffer phase = straighter paths. |
| **Lump Volume** ($R$) | **Confinement** (Bowl) | The spatial domain where the field is excited. |
| **Phase Correlation** | **Chain Thickness** ($a$) | Transverse coherence width; acts as excluded volume. |
| **Noise Intensity** ($\rho_N$) | **Temperature** ($T$) | Drives the "wiggliness" of the paths. |
| **Vortex Loop** ($l=1$) | **Ring Polymer** | A closed path with $2\pi$ winding (Lepton-like). |
| **Y-Junction** ($l=1$ composite) | **Star Polymer** (3-arm) | Three paths meeting at a core (Baryon-like). |

---

## 2. Theoretical Models

We treat the phase paths as **Worm-Like Chains (WLC)** under spherical confinement.

### A. The Chain Model (WLC)
Unlike an ideal random walk, the WLC includes bending stiffness.
- **Key Parameter:** Persistence length $\ell_p \sim \frac{\kappa w_*^2}{\text{Noise Power}}$.
- **Regimes:**
    - $\ell_p \gg R$: Stiff rods. Defects are rare/suppressed.
    - $\ell_p \ll R$: Flexible coils. High entropy, many defects possible.

### B. The Observables (What counts as a "Particle"?)

#### 1. Single Loop Formation (Lepton Channel)
*Polymer equivalent:* End-to-end closure probability (Cyclization).
*Physics:* Probability that a phase path wraps around and reconnects to form a stable vortex ring.
*Scaling:*
$$ P_{\text{loop}}(L) \sim \left(\frac{a}{L}\right)^{3/2} \exp\left(-\frac{\ell_p}{L}\right) $$
*Interpretation:* Short loops are stiff-suppressed; long loops are entropy-suppressed. Confinement ($R$) enhances closure probability by reducing the available volume.

#### 2. Three-Leg "Tripod" Formation (Baryon Channel)
*Polymer equivalent:* Formation of a 3-arm Star Polymer (branched junction).
*Physics:* Probability that three distinct phase paths meet at a common "core" with compatible orientations (color neutrality).
*Scaling:*
$$ P_{\text{3-arm}} \sim R^{-\gamma_3} $$
where $\gamma_3$ is the star-polymer vertex exponent.
*Interpretation:* This is entropically suppressed relative to single loops, but **confinement** in a dense lump increases the density of near-contacts, boosting the branching fraction.

---

## 3. Simulation Strategy: The "Hybrid" Approach

Instead of running expensive full-field PDE simulations for every point in the $(R, \sigma_\phi)$ parameter space, we use the polymer scaling laws as an **analytical interpolator**.

### Step 1: Analytical Skeleton (The Prior)
Define the theoretical selection function based on WLC statistics:
$$ P(\text{Type} \mid R, \sigma_\phi) \approx A \cdot f_{\text{conf}}(R/\ell_p) \cdot g_{\text{topo}}(\text{Type}) $$
- $f_{\text{conf}}$: Confinement function (sigmoid-like).
- $g_{\text{topo}}$: Topological weight (1 for loop, $\epsilon_{\text{branch}}$ for tripod).

### Step 2: Calibration (The Likelihood)
Run a sparse grid of full PDE simulations (Stage 1 of existing plan) to measure the actual formation rates at key points.
- **Fit the prefactors:** Determine $A$ and $\epsilon_{\text{branch}}$ from the data.
- **Calibrate $\ell_p$:** Map the simulation's $\sigma_\phi$ (noise) to the effective persistence length $\ell_p$.

### Step 3: Synthesis (The Posterior)
Use the calibrated polymer model to generate the high-resolution selection function needed for the Stage 2 cosmological yield calculation.

---

## 4. Open Questions & Iteration Goals

1.  **Chiral Bias:** The polymer model is naturally symmetric. How do we insert the "Chiral Tilt" ($\theta$-term)?
    *   *Idea:* Add a "torsion" term to the WLC energy that favors one winding direction (helical bias).
2.  **Junction Stability:** Polymer stats tell us a junction *forms*, but not if it *survives*.
    *   *Idea:* Multiply the formation probability by a stability factor derived from the "Bulk Resistance" potential (from Particles paper).
3.  **Dense Packing:** What happens when the "bowl" is packed with noodles (high noise)?
    *   *Idea:* Move from single-chain statistics to Flory mean-field theory for polymer melts.

## 5. Next Steps
- [ ] Implement a simple 3D WLC Monte Carlo script to generate the baseline $P(R)$ curves.
- [ ] Define the exact "contact criteria" for a Y-junction in the MC model.
- [ ] Run the calibration against the first batch of PDE results.
