# Open Issues / Research Tasks

This file tracks **open conceptual + technical issues** that are not yet fully pinned down,
but which (if resolved) would substantially tighten the framework, reduce free parameters,
and create cross-paper consistency checks.

Keep entries short, actionable, and tied to concrete places in the papers/simulations.

---

## Issue 1 — Make cross-window agreement explicit: exact re-windowing covariance (RG matching of observables)

### What we are already doing implicitly
Across the framework we repeatedly use a hidden EFT principle:

- **Overlapping windows must agree** on the same physical picture in their overlap domain (same invariants; same sign/ordering of trends; no contradictory runaways), even if they use different effective variables and renormalised couplings.

This is “windowing” as more than epistemics: it is **consistency of overlapping effective descriptions**.

### The tighter restriction to explore
Promote the implicit overlap-agreement to an explicit, stronger requirement:

- **Exact re-windowing covariance:** under coarse-graining / re-windowing, the *functional form* of the effective equations is preserved; only parameters/couplings run.

Equivalently: the effective theory lies in a universality class stable under re-windowing (an RG-style statement).

### Why this matters (payoff)
If correct, this turns many “soft” claims into sharp constraints:

- **Over-determination:** the *same* small set of running parameters must satisfy Gravity + Forces + Particles + BH constraints simultaneously.
- **Predictive cross-links:** infer running in one sector (or simulation) → predict consequences in others.
- **Parameter bounding:** combine independent simulation/observation bounds into a single intersecting posterior on the underlying parameters.

### What is not yet clear (the actual open question)
Exact covariance can mean different strengths. We need to decide (and justify) which one holds:

- **Weak form (default):** only *dimensionless observables* are invariant; dimensionful ones co-scale approximately. This is the baseline expectation unless the stronger claim can be derived for the specific regime from the graph + free-energy dynamics.
- **Strong form (conditional):** the **effective action / equations** are form-invariant under re-windowing; couplings run by beta functions. This should be treated as a hypothesis to be checked in concrete scenarios (it may hold in some windows/eras and fail in others).
- **Near-pivot involution (optional):** an additional “self-duality-like” involution under window exchange that could fix a near-unique crossover shape (hyperbolic bounce), rather than just a smooth class.

The issue is to determine which strength is warranted by the framework’s mechanics (window saturation + free-energy minimization + nesting) and which is an extra axiom.

Practical rule of thumb: adopt the weak form as the “always on” consistency constraint; invoke the strong form only after verifying the dynamical prerequisites in the target regime (adiabatic window drift, clear scale separation, stable fixed-point behaviour of the coarse-grained functional, etc.).

### Concrete tasks / checks (high value)

#### A) Identify the minimal parameter set and their “running targets”
Make an explicit list of the few parameters that appear everywhere and may run with window / background noise proxy:

- **phase-sector**: \(\kappa\), \(w_*\) (with \(c_s^2=\kappa w_*^2\))
- **amplitude-sector**: \(\alpha_{\rm grad}\), \(\beta_{\rm pot}\), \(\eta_0\), (and any curvature/volume stiffnesses in Particles)
- **window functionals**: \(k_{\min},k_{\max}\), \(\Lambda_{\rm eff}\), \(A_{\rm win}\), etc.
- **noise proxy**: \(\rho_N\) / \(\tau^2\) and the mapping to \(\Phi\) in the Newtonian window

Goal: express the candidate covariance principle as beta-function constraints for these quantities.

#### B) Gravity: induced action matching as an RG constraint
In `exploratory/Gravity/SI.tex`, windowing enters explicitly via \(\Lambda_{\rm eff}\) and \(A_{\rm win}\).

Check whether “form-invariant induced action” implies a specific running relation such as:

- \(C_1\Lambda_{\rm eff}^2 \sim (16\pi G)^{-1}\) with fixed \(C_1\), hence \(G \propto \Lambda_{\rm eff}^{-2}\) (up to dimensionless constants),

and whether that is compatible with the operational definition of \(\Lambda_{\rm eff}\) in terms of the observation window and spectrum.

Deliverable: one page derivation sketch (“if exact covariance, then…”), and a list of admissible \(S(k)\)/window definitions consistent with it.

#### C) Forces: exact cone invariance and constitutive relations
Forces already asserts:

- \(\varepsilon_0^{\rm eff}\mu_0^{\rm eff}=1/c_s^2\) with mild window dependence.

Exact covariance would sharpen:

- **cone invariance** (\(\alpha_c = 0\) at the fixed-form level), with finite-window deviations as controlled corrections,
- correlated running of screening/coherence lengths (Yukawa \(\ell\), phase coherence \(\lambda_*\)) through \(\alpha_{\rm grad},\beta_{\rm pot},\kappa,w_*\).

Deliverable: explicit “RG consistency identities” tying the runnings so that Coulomb \(1/r\), Yukawa screening, and confinement \(+\sigma r\) retain their *form* across windows.

#### D) Particles: window-running of mass matrices becomes computable
Particles already says changing experimental scale is changing the coarse-graining window.

Exact covariance would require:

- the mass-matrix overlap form stays fixed; \(\mathcal W_f\) coefficients run in a controlled way from the same underlying beta functions,
- mass ratios / mixing angles as RG invariants (up to finite-window corrections) are not just “likely stable” but *forced*.

Deliverable: rewrite the “mild running” statement as an explicit RG-covariance statement; derive \(\delta_f(\rho_N)\) from parameter runnings rather than leaving it phenomenological.

#### E) Black holes / pivot: does exact covariance + optional involution fix the crossover shape?
BH work already shows that window saturation + minimization implies a smooth U-turn class.

Question: if we adopt exact re-windowing covariance (and optionally an involutive UV↔IR exchange), do we get:

- a more unique crossover shape (e.g. harmonic-mean / quadrature “hyperbolic bounce”),
- and does that choice propagate to other sectors (e.g. particle mass-radius scaling, screening-length scaling)?

Deliverable: a short “closure taxonomy” (generic smooth class vs harmonic mean vs involutive duality), with what each assumption buys and costs.

### Outputs we want (so simulations can constrain it)
- A **small set of RG-consistency equations** linking the running of \(\kappa,w_*,\alpha_{\rm grad},\beta_{\rm pot},\eta_0,\dots\)
- A mapping from measurable/simulatable quantities to those parameters:
  - filament width / void profile ↔ screening length ↔ \(\alpha_{\rm grad}/\beta_{\rm pot}\)
  - EM constitutive parameters ↔ \(\kappa,w_*,\beta_{\rm pot}\)
  - mass ratios / mixing drifts ↔ stiffness ratios / \(\mathcal W_f\) running
  - (BH side) pivot crossover sharpness ↔ saturation functional form

### Failure modes (things to watch for)
- Exact covariance may only hold **in an intermediate window** (not all the way to extremes) due to non-adiabatic drift of window edges.
- Different sectors may sit in different universality classes unless additional assumptions enforce unification.
- The operational definition of \(\Lambda_{\rm eff}\) (and/or \(A_{\rm win}\)) may be too ambiguous: we may need a sharper micro-to-macro prescription for the window functionals.

