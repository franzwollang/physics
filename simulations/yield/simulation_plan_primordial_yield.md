# Simulation Task: Primordial Yield via The Inverse Problem

---

## 0. Parameter Interface (Canonical Microphysics)

This simulation suite solves the **Inverse Problem** of Dark Matter genesis. Instead of running a massive cosmological simulation and hoping for the right result, we break the chain into three testable links to determine *what kind of primordial noise* is required to produce the observed universe.

## 0.A Scope, Assumptions, and “Knobs” (Read This First)

This plan is written to be **directly implementable by a subagent**. It is also explicit about where modeling choices enter.

### 0.A.1 Scope (what we are / aren’t simulating)

- **We are simulating:** defect induction probability of a localized lump under controlled phase-gradient forcing (Stage 1), and the inverse statistical mapping from a primordial spectrum family to an energy-weighted baryon/DM ratio (Stage 2), and a graph-spectrum plausibility check using vacuum-relaxed Laplacians (Stage 3).
- **We are not simulating:** full cosmological expansion, plasma microphysics, baryogenesis chemistry, or detailed hadron structure. “Baryon” here means **\(l=1\) topological rotor** vs “DM” as **\(l=0\) non-rotor** in the model’s soliton taxonomy.

### 0.A.2 The minimal invariants we care about (must be robust)

- **Stage 1:** a stable empirical map \(P(l=1\mid R,\sigma_\phi)\) with confidence bands.
- **Stage 2:** an **energy-weighted** ratio \(E_{DM}/E_B\) matching \(\approx 5.3\) for some spectrum family; sensitivity of the fit to cutoffs and intermittency must be reported.
- **Stage 3:** existence (or non-existence) of vacuum-relaxed graphs whose mid-band Laplacian spectrum can reproduce the Stage-2 implied slope *in a declared fit window*.

### 0.A.3 “Knobs” that must be fixed (defaults provided; subagent must not guess silently)

**Stage 1 knobs (PDE + forcing):**

- Spatial dimension: **2D only** for Stage 1. (3D is forbidden in this task; do not “skip ahead”.)
- Boundary condition type: default **Dirichlet phase + pinned amplitude**.
- Forcing decomposition: coherent \(\vec{k}\) (default \(\vec{k}=0\) for pure noise scan; later add \(|\vec{k}|\neq 0\)) + boundary noise \(\eta\).
- Boundary-noise spectrum: \(\sigma_m\propto m^{-q}\) with default \(q=1\), and cutoff \(m_{\max}\) default 64.
- Induction success criteria: persistence time \(T_{\rm persist}\), core threshold \(\epsilon_{\rm core}\), and loop radius for winding measurement.

**Stage 2 knobs (statistics):**

- Spectrum family: fixed in Section 3.1 (no alternatives).
- Lump candidate definition: fixed in Section 3.1 (no alternatives).
- Energy proxy: fixed \(E(R)=R^3\).
- Intermittency / lognormal / “rogue waves”: forbidden in this implementation.

**Stage 3 knobs (graph→spectrum mapping):**

- Vacuum graph generator: fixed in Section 4.1 (must use `simulations/vacuum/vacuum_geometry_simulation.py`).
- Mapping eigenvalues to “wavenumbers”: fixed \(k=\sqrt{\lambda}\).
- Phase model: fixed OU model in Section 4.1 (no alternatives).
- Fit band: fixed quantile band in Section 4.1 (no alternatives).

### 0.A.5 Anti-gaming rules (non-negotiable)
If you implement this pipeline, you must follow the exact numbers in Sections 2.3 / 3.1 / 4.1. The following are explicitly forbidden:

- Reducing any grid size or number of realizations/trials.
- Using fewer bootstrap resamples.
- Replacing the grid search with an “optimizer” to reduce compute.
- Replacing the OU graph model with a different heuristic.
- Skipping the boundary calibration step in Stage 1.
- “Testing on fewer points” and then extrapolating.
- Reporting only the best-fit without saving the raw aggregates needed to verify it.

If your computer cannot run the full spec, you must report “INSUFFICIENT COMPUTE TO RUN SPEC” and stop. Do not silently downscale.

### 0.A.6 Audit artifacts (must be produced so results are verifiable)
Each stage must write:

- A JSON result file (schemas in Section 6).
- A `run_manifest.json` containing: git hash (if available), timestamp, OS, Python version, CPU model (if available), and the full CLI arguments used.
- A `raw_aggregates.npz` file containing the minimum raw arrays needed to recompute the reported metrics (details per stage below).

### 0.A.4 Deliverables (ungamable pass/fail checklist)

Your implementation is accepted ONLY if all items below are true:

- **Stage 1 output files**
  - `stage1_selection.json` exists and matches Section 6.1.
  - `run_manifest.json` exists.
  - `raw_aggregates.npz` exists and contains the required arrays with required shapes (Section 2.4).

- **Stage 2 output files**
  - `stage2_fit.json` exists and matches Section 6.2.
  - `run_manifest.json` exists.
  - `raw_aggregates.npz` exists and contains the required arrays with required shapes (Section 3.2).

- **Stage 3 output files**
  - `stage3_graph_check.json` exists and matches Section 6.3.
  - `run_manifest.json` exists.
  - `raw_aggregates.npz` exists and contains the required arrays with required shapes (Section 4.1.5).

- **No forbidden shortcuts**
  - Any deviation from the fixed numbers in Sections 2.3 / 3.1 / 4.1 is an automatic failure.
  - Any reduction of trials / realizations / grid size / bootstrap resamples is an automatic failure.

- **Acceptance tests met**
  - Stage 1 passes Section 2.5 acceptance tests.
  - Stage 2 passes Section 3.2 acceptance test.
  - Stage 3 is evaluated against Stage 2 CI exactly per Section 4.1.5.

### 0.1 Canonical Parameters

| Symbol | Name | Description |
| ------ | ---- | ----------- |
| $\alpha_{\rm grad}$ | Gradient Stiffness | Energetic cost of domain walls. |
| $\beta_{\rm pot}$ | Potential Curvature | Depth of the symmetry-breaking well. |
| $\gamma$ | Quartic Coupling | Stabilizes the Mexican-hat potential; sets \(w_*=\sqrt{\beta_{\rm pot}/(2\gamma)}\). |
| $\kappa$ | Phase Stiffness | Emergent from local graph connectivity. |
| $\vec{v}_{phase}$ | Phase Wind | The local gradient of the phase field ($\nabla \phi$) representing noise/flows. |

**Derived scales (used to nondimensionalize Stage 1):**
- \(w_*=\sqrt{\beta_{\rm pot}/(2\gamma)}\) (vacuum amplitude)
- \(\ell=\sqrt{\alpha_{\rm grad}/(2\beta_{\rm pot})}\) (amplitude correlation length / core-width proxy)
- \(c_s=\sqrt{\kappa}\,w_*\) (phase-wave speed)

**Conventions for this simulation suite (mandatory):**
- Stage 1 is run in **dimensionless units** with \(w_*=1\), \(\ell=1\), \(c_s=1\), and must report the mapping to \((\alpha_{\rm grad},\beta_{\rm pot},\gamma,\kappa)\).
- \(\sigma_\phi\) in Stage 1 is the **RMS phase-gradient magnitude** at the soliton boundary, i.e. \(\sigma_\phi := \sqrt{\langle|\nabla\phi|^2\rangle_{\text{shell}}}\) with “shell” defined below.

### 0.2 Stage Interfaces (what passes between stages)

- **Stage 1 → Stage 2:** \(P(l=1\mid R,\sigma_\phi)\) map + its confidence bands + the *exact definition of* \(R\) and \(\sigma_\phi\) (shell definition, measurement method).
- **Stage 2 → Stage 3:** target spectrum family and best-fit parameters with credible intervals, plus the fit band/cutoffs.

### 0.3 Failure Modes (STOP / REPORT rules; no discretion)

- **Stage 1 never nucleates defects (fails Section 2.5 acceptance tests):**
  - STOP immediately.
  - DO NOT tweak parameters, reduce thresholds, or change the algorithm.
  - REPORT exactly:
    - `stage1_selection.json`
    - `run_manifest.json`
    - `raw_aggregates.npz`
    - and a short text note: “STAGE1_FAIL_NO_NUCLEATION”

- **Stage 1 nucleates defects everywhere (also fails Section 2.5 acceptance tests):**
  - STOP immediately.
  - DO NOT tweak parameters or weaken success criteria.
  - REPORT exactly the same artifacts as above and the note: “STAGE1_FAIL_ALWAYS_NUCLEATES”

- **Stage 2 is wildly sensitive or cannot hit the target band (fails Section 3.2 acceptance test):**
  - STOP immediately.
  - DO NOT add intermittency/lognormal terms, do not change peak-finding, do not change the spectrum family, and do not change the grids.
  - REPORT:
    - `stage2_fit.json`
    - `run_manifest.json`
    - `raw_aggregates.npz`
    - and the note: “STAGE2_FAIL_FIXED_MODEL_INCONSISTENT”

- **Stage 3 mismatches (no kept vacuum run matches Stage-2 CI):**
  - STOP immediately.
  - DO NOT widen fit bands, change the OU model, change the histogram bins, or change the sweep.
  - REPORT:
    - `stage3_graph_check.json`
    - `run_manifest.json`
    - `raw_aggregates.npz`
    - and the note: “STAGE3_FAIL_GRAPH_SPECTRUM_MISMATCH”

---

## 1. Theoretical Context: The Three-Stage Chain

We hypothesize that Baryonic Matter ($l=1$) forms from scalar lumps ($l=0$) that are "spun up" by the ambient phase noise via the Hairy Ball Theorem. Dark Matter is simply the population of lumps that failed to spin up.

**The Chain of Causality:**
1. **Induction (Micro-Physics):** A scalar lump sits in a noisy phase bath ("Phase Wind"). If the wind is strong/coherent enough, it forces the lump to acquire a topological defect ($l=0 \to l=1$) to satisfy the Hairy Ball boundary condition.
2. **Statistics (The Spectrum):** The universe is filled with a spectrum of lump sizes and noise intensities (e.g., Rogue Wave statistics).
3. **Origin (Topology):** This spectrum originates from the quench on the Infinite-Clique Graph.

**The Strategy:**
Work backwards. Determine the micro-physics (Stage 1), find the required statistical spectrum to get the 5:1 ratio (Stage 2), and then check if the Graph produces that spectrum (Stage 3).

---

## 2. Stage 1: The Induction Test (Micro-Physics)

**Goal:** Determine the **Selection Function** \(P(l=1 \mid R, \sigma_\phi)\).

- Does a scalar lump of radius $R$ in a phase wind of strength $\sigma_\phi$ turn into a Baryon ($l=1$) or remain Dark Matter ($l=0$)?

**Important clarification (implementation):**
- Stage 1 does *not* assume the barrier criterion is exact; it uses it only to motivate a phase diagram.
- The output must therefore be empirical: for each \((R,\sigma_\phi)\), estimate \(p=\Pr[n\neq 0]\) over \(N_{\rm trials}\) randomized initial noise realizations.

### 2.1 Model: SGLE / dissipative GPE (define the PDE precisely)
Use the split fields \(w(x,t)\) and \(\phi(x,t)\) so \(\alpha_{\rm grad}\) and \(\kappa\) can be set independently. (The shortcut \(\partial_t\Psi=\alpha\nabla^2\Psi+\dots\) implicitly ties phase stiffness to amplitude stiffness and is forbidden for this simulation suite.)

**Fixed dimensionless parameters (no tuning in Stage 1):**
- \(\alpha_{\rm grad}=2\)
- \(\beta_{\rm pot}=1\)
- \(\gamma=1/2\)
- \(\kappa=1\)

This enforces \(w_*=1\), \(\ell=1\), \(c_s=1\).

**Unknowns:**
- \(w(x,t)\ge 0\)
- \(\phi(x,t)\in\mathbb R\)
- \(\Psi(x,t)=w(x,t)\,e^{i\phi(x,t)}\)

**Evolution equations (implement exactly):**
\[
\partial_t w \;=\; \alpha_{\rm grad}\left(\nabla^2 w - w|\nabla\phi|^2\right) + \beta_{\rm pot} w - 2\gamma w^3
\]
\[
\partial_t \phi \;=\; \kappa\left(\nabla^2\phi + 2\,\frac{\nabla w}{w+\epsilon_w}\cdot \nabla\phi\right)
\]
with \(\epsilon_w=10^{-6}\).

**Dimensionality:** Stage 1 is **2D only**. Do not implement 3D until Stage 1 passes the acceptance tests in Section 2.5.

### 2.2 Theoretical Requirement: The Energy Barrier (motivation only)
The transition from a Scalar Breather ($l=0$) to a Vector Rotor ($l=1$) is a **topological phase transition**. It is not a smooth rotation; it requires the field amplitude $|\Psi|$ to pass through zero to create the defect core.

- **The Barrier:** To create a defect, the system must "punch a hole" in the condensate. This costs potential energy proportional to the well depth $\Delta V \sim \beta_{\rm pot}^2 / (4\gamma)$.
- **The Driver:** The "Phase Wind" (gradient noise) provides the kinetic energy density $\mathcal{E}_K \sim \frac{1}{2} \kappa w^2 (\nabla \phi)^2$.
- **The Collapse Condition:** A transition is only possible if the local kinetic stress exceeds the potential barrier:

  $$ \frac{1}{2} \kappa w^2 (\nabla \phi)^2 > V_{\text{barrier}} $$

  This implies a critical "Phase Reynolds Number" or noise threshold below which the lump simply stretches but never snaps into a vortex.

### 2.3 Protocol (idiot-proof, no choices)
Follow these steps exactly.

#### 2.3.1 Fixed domain and discretization
- Dimension: \(d=2\).
- Domain: \(\Omega=[-40,40]\times[-40,40]\) (so \(L=80\)).
- Grid: \(N=512\) per axis (so \(\Delta x = 80/512 \approx 0.15625\)).
- Time step: \(\Delta t = 5\times 10^{-4}\).
- Finite differences: 2nd-order central differences for gradients and Laplacians.

#### 2.3.2 Fixed time schedule
- Relax time: `T_relax = 20.0` (no forcing).
- Drive time: `T_drive = 80.0` (forcing on).
- Measurement cadence: compute winding every `dt_sample = 0.2`.
- Persistence requirement: `T_persist = 10.0` (must be continuously nonzero for this long).

#### 2.3.3 Initial condition (single \(l=0\) lump)
For each radius parameter \(R\) (defined below), initialize:
- Center at \(x_0=(0,0)\).
- Phase: \(\phi(x,0)=0\).
- Amplitude:
  \[
  w(r,0)=1 - 0.6\exp\!\left(-\frac{r^2}{2R^2}\right),\quad r=\|x\|
  \]

#### 2.3.4 Boundary forcing (this defines the “phase wind”)
Force only through boundary conditions; no bulk forcing terms.

- Amplitude boundary condition: \(w|_{\partial\Omega}=1\).
- Phase boundary condition:
  \[
  \phi|_{\partial\Omega}(s) = \eta(s)
  \]
  where \(s\) is arclength along the boundary (total boundary length \(L_b=4L=320\)).

Construct \(\eta(s)\) as:
\[
\eta(s)=\sum_{m=1}^{64}\Big(a_m\cos(2\pi m s/L_b)+b_m\sin(2\pi m s/L_b)\Big),
\]
with \(a_m,b_m \sim \mathcal N(0,\sigma_m^2)\) and \(\sigma_m = A_\eta\, m^{-1}\).

**Calibration rule (do not skip):** do not choose \(A_\eta\). Calibrate to hit \(\sigma_\phi\) target using this deterministic 2-pass rescaling:
1. Set \(A_\eta=1\), run only the phase equation for `T_cal = 5.0` with \(w\equiv 1\) fixed, compute the achieved \(\sigma_{\phi,\rm meas}\) (definition below).
2. Rescale all \(a_m,b_m\) by factor \(g = \sigma_{\phi,\rm target}/(\sigma_{\phi,\rm meas}+10^{-12})\).
3. Proceed with the full coupled run (relax + drive) using the rescaled boundary \(\eta\).

#### 2.3.5 Definition of control parameter \(\sigma_\phi\) (must match Stage 2)
Given a run with a lump radius parameter \(R\), define:
- Shell annulus: \(r\in [R, R+2]\) (since \(\ell=1\)).
- Compute \(|\nabla\phi|\) on the grid by finite differences.
- Define \(\sigma_\phi\) as the RMS over all grid points in the shell:
  \[
  \sigma_\phi := \sqrt{\langle |\nabla\phi|^2\rangle_{\text{shell}}}.
  \]
Measure and store \(\sigma_\phi\) immediately after calibration and at the start of the driven coupled run.

#### 2.3.6 Evolution schedule (exact)
For each trial:
1. **Relax (no forcing):** integrate coupled \(w,\phi\) with boundary conditions \(w|_{\partial\Omega}=1\), \(\phi|_{\partial\Omega}=0\) for `T_relax`.
2. **Calibrate boundary forcing:** run the phase-only calibration described above to set \(\eta(s)\) amplitude for the target \(\sigma_\phi\).
3. **Drive:** integrate coupled \(w,\phi\) for `T_drive` with boundary conditions \(w|_{\partial\Omega}=1\), \(\phi|_{\partial\Omega}=\eta(s)\).

#### 2.3.7 Defect detection (exact algorithm)
Define a fixed circular loop \(C\) at radius:
- `r_loop = R + 1.0`

Algorithm to compute winding at a given time:
1. Sample `M_loop = 512` points uniformly on the circle.
2. At each sample point, bilinearly interpolate \(\Psi=w e^{i\phi}\) from the grid.
3. Compute phase \(\phi=\mathrm{atan2}(\Im\Psi,\Re\Psi)\).
4. Compute successive differences \(\Delta\phi_j\), wrap each into \((-\pi,\pi]\), and sum: \(\Delta\Phi=\sum_j \Delta\phi_j\).
5. Winding estimate: \(n=\mathrm{round}(\Delta\Phi/(2\pi))\).

**Core confirmation:** declare a vortex “real” only if, at the same time sample,
- `min_absPsi_in_disk(R/2) < epsilon_core` with `epsilon_core = 0.2`.

**Persistence rule:** declare “success” if there exists a contiguous time interval of length `T_persist` during the driven phase where:
- \(n\neq 0\) at every sample time, AND
- core confirmation holds at least once inside that interval.

#### 2.3.8 Parameter grid and trials (fixed)
Compute the selection function on a fixed grid:
- Radii: `R_values = 30` log-spaced values in \([0.5, 8.0]\).
- Targets: `sigma_phi_values = 30` linear values in \([0.0, 3.0]\).
- Trials per grid point: `N_trials = 200`.
- RNG seeding: seed each trial deterministically as `seed = base_seed + hash(R_idx, sigma_idx, trial_idx)`.

For each grid point, compute:
- \(\hat p = \#\text{success}/N_{\rm trials}\)
- Wilson 95% CI \([p_{\rm low}, p_{\rm high}]\)

### 2.4 Output: The Phase Diagram (what to write to disk)
Write a single JSON file exactly matching the schema in Section 6.1, with:
- the full \(\hat p(R,\sigma_\phi)\) matrix,
- CI low/high matrices,
- all numeric defaults in this section embedded in the JSON under `solver` and `protocol`,
- and a small set of representative debug artifacts:
  - save `n_time_series.npy` and `min_absPsi_time_series.npy` for **one** representative run in each qualitative regime (low/med/high \(\sigma_\phi\)).

**Stage-1 raw artifacts (must write, ungamable):**
- `raw_aggregates.npz` must include:
  - `R_values` (shape `(30,)`)
  - `sigma_phi_values` (shape `(30,)`)
  - `success_counts` (shape `(30,30)`, integer)
  - `trial_counts` (shape `(30,30)`, integer; must be exactly 200 everywhere)
  - `sigma_phi_measured_calibration` (shape `(30,30,200)`; the achieved \(\sigma_\phi\) after calibration per trial)
  - `sigma_phi_measured_start_drive` (shape `(30,30,200)`; \(\sigma_\phi\) at start of drive per trial)
  - `n_final` (shape `(30,30,200)`; winding at end of drive)
  - `min_absPsi_min_over_time` (shape `(30,30,200)`; minimum \(|\Psi|\) observed during drive)

If any of these arrays are missing or have different shapes, the run is invalid.

### 2.5 Stage-1 Implementation Checklist (must pass)

- **Numerics**
  - [ ] stable time integrator (explicit Euler acceptable initially; semi-implicit preferred later)
  - [ ] phase unwrapping tested on synthetic vortices
  - [ ] unit tests for winding computation on known fields
- **Diagnostics**
  - [ ] save a few “movie frames” (or snapshots) for representative points in each regime
  - [ ] log whether failures are “no vortex,” “transient vortex,” or “persistent vortex”
  - [ ] verify conservation/relaxation sanity: without forcing, lump stays \(n=0\)
- **Sampling**
  - [ ] per-point trials = 200
  - [ ] report CI low/high matrices

**Acceptance tests (hard, no discretion):**
- With \(\sigma_\phi=0\), success rate must be \(\hat p < 0.05\) for all \(R\).
- For at least one \(\sigma_\phi \ge 2.0\) and some \(R\), success rate must satisfy \(\hat p > 0.30\).
- If either condition fails, you must fix the bug in forcing calibration or winding detection before proceeding to Stage 2.

---

## 3. Stage 2: The Rogue Wave Fit (Statistics)

**Goal:** Find a **primordial phase-spectrum family** that yields \(\Omega_{DM} / \Omega_B \approx 5\) after applying the Stage-1 selection function.

**Critical detail (conceptual + implementation):** the observed ratio is in **energy density**, not number counts. Stage 2 must therefore weight outcomes by an energy proxy per lump:
- Energy proxy is fixed for this simulation: \(E(R)=R^3\). (We are fitting an energy-density ratio; we use a 3D volume proxy even though sampling is done on a 2D field.)
- You must report both “count ratio” and “energy-weighted ratio”; match the energy-weighted one to \(\sim 5.3\).

### 3.1 Protocol (idiot-proof, no choices)
Follow these steps exactly.

#### 3.1.1 Fixed domain and grids
- Domain: periodic square \([0, 80)\times[0, 80)\). (Matches Stage-1 length scale \(L=80\).)
- Grid: `N = 1024` per axis.
- Wavenumbers: use FFT frequencies \(k_x,k_y = 2\pi\cdot\text{fftfreq}(N, d=80/N)\) and \(k=\sqrt{k_x^2+k_y^2}\).

#### 3.1.2 Fixed spectrum family and fixed cutoffs
Stage 2 fits exactly two parameters:
- spectral index \(n_s\)
- gradient-scale factor \(s_\sigma\) (a global multiplier that rescales \(\phi\) so typical \(\sigma_\phi\) falls into Stage-1’s \([0,3]\) range)

The raw GRF spectrum is:
\[
P_\phi(k) = k^{-n_s}\,\mathbf 1_{[k_{\min},k_{\max}]}(k)
\]
with fixed cutoffs:
- `k_min = 2π / 80` (the fundamental mode)
- `k_max = π / (80/N)` (the Nyquist mode)

#### 3.1.3 Exact GRF synthesis for \(\phi(x)\)
For each realization:
1. Draw complex Fourier coefficients \(\hat\phi(\mathbf k)\) with Hermitian symmetry so \(\phi(x)\) is real.
2. For each \(\mathbf k\neq 0\), set variance:
   \[
   \mathbb E|\hat\phi(\mathbf k)|^2 = P_\phi(k).
   \]
3. Set \(\hat\phi(\mathbf 0)=0\).
4. Inverse FFT to obtain \(\phi(x)\).
5. Compute \(\nabla\phi\) by spectral differentiation (\(\partial_x\phi\leftrightarrow i k_x\hat\phi\), \(\partial_y\phi\leftrightarrow i k_y\hat\phi\)).

**Mandatory normalization step (removes the arbitrary amplitude \(A\)):**
- Compute the global RMS \(g := \sqrt{\langle |\nabla\phi|^2\rangle_{\Omega}}\).
- Replace \(\phi \leftarrow \phi/(g+10^{-12})\) so that global RMS gradient is exactly 1.
- Then apply the fit scale: \(\phi \leftarrow s_\sigma \phi\).

#### 3.1.4 Exact seed field \(u(x)\) for lump candidates
Generate an independent GRF \(u(x)\) with *fixed* white spectrum:
\[
P_u(k) = \mathbf 1_{[k_{\min},k_{\max}]}(k),
\]
same cutoffs and grid as \(\phi\). Set \(\hat u(\mathbf 0)=0\). Inverse FFT to real \(u(x)\).

#### 3.1.5 Exact peak finding (lump candidates)
Find candidate peaks as local maxima on the grid:
- A pixel \((i,j)\) is a peak if \(u_{i,j}\) is strictly greater than its 8 neighbors.
- Threshold: keep only peaks with \(u_{i,j} > \mu_u + 2\sigma_u\) where \(\mu_u,\sigma_u\) are the global mean/std of \(u\).
- If the number of peaks is < 2000, re-run threshold at \(\mu_u+1.5\sigma_u\).
- If still < 2000, take the top 2000 points by \(u\) value and then discard any that are within Chebyshev distance ≤ 3 of a higher-ranked point (deterministic non-maximum suppression).
- If more than 5000 peaks remain, keep the top 5000 by \(u\).

These surviving peaks are the lump centers \(x_i\).

#### 3.1.6 Exact radius estimator \(R_i\)
At each peak \(x_i\), estimate radius from curvature:
- Compute discrete Laplacian \(\Delta u(x_i)\) by 5-point stencil.
- Set
  \[
  \kappa_i = \max\!\left(10^{-6},\; -\frac{\Delta u(x_i)}{u(x_i)+10^{-12}}\right),
  \qquad
  R_i = \mathrm{clip}\!\left(\kappa_i^{-1/2},\; 0.5,\; 8.0\right).
  \]

#### 3.1.7 Exact \(\sigma_{\phi,i}\) estimator (must match Stage 1)
For each peak \(x_i\) with radius \(R_i\):
- Shell annulus is \(r\in [R_i, R_i+2]\).
- Compute \(|\nabla\phi|\) on the grid.
- Compute RMS in the shell:
  \[
  \sigma_{\phi,i} := \sqrt{\langle|\nabla\phi|^2\rangle_{\text{shell around }x_i}}.
  \]

#### 3.1.8 Apply Stage-1 selection map (exact interpolation)
Load Stage-1 JSON and interpolate \(p_i=P(l=1\mid R_i,\sigma_{\phi,i})\) by bilinear interpolation over the grid.
Clamp \(p_i\) into \([0,1]\).

Compute yields:
- Count-weighted:
  \[
  N_B=\sum_i p_i,\qquad N_{DM}=\sum_i (1-p_i).
  \]
- Energy-weighted (fixed \(E(R)=R^3\)):
  \[
  E_B=\sum_i R_i^3\,p_i,\qquad E_{DM}=\sum_i R_i^3\,(1-p_i).
  \]

#### 3.1.9 Fit procedure (fixed grid search, no optimizer)
Run a brute-force grid search over:
- `n_s_grid = [0.5, 0.55, ..., 4.0]` (step 0.05)
- `s_sigma_grid = logspace(log10(0.25), log10(4.0), 25)`

For each \((n_s,s_\sigma)\):
- Generate `N_realizations = 200` independent pairs \((\phi,u)\).
- Aggregate \(E_B,E_{DM},N_B,N_{DM}\) across realizations.
- Compute the objective:
  \[
  \mathcal L(n_s,s_\sigma)=\left(\log\frac{E_{DM}}{E_B}-\log 5.3\right)^2.
  \]

Pick the minimizer \((n_s^\*,s_\sigma^\*)\).

#### 3.1.10 Uncertainty (fixed bootstrap)
Bootstrap over realizations:
- draw 500 bootstrap resamples of the 200 realizations (with replacement),
- recompute \((n_s^\*,s_\sigma^\*)\) by minimizing \(\mathcal L\) on the same grids,
- report 16–84% and 2.5–97.5% intervals for \(n_s^\*\) and \(s_\sigma^\*\), and for the final \(E_{DM}/E_B\).

### 3.2 Output: The Required Spectrum (what to write to disk)
Write a single JSON file exactly matching the schema in Section 6.2, including:
- `best_fit.n_s` and `best_fit.s_sigma`
- fixed `k_min` and `k_max`
- `yield_ratio_energy = E_DM/E_B` and `yield_ratio_count = N_DM/N_B`
- bootstrap intervals for both ratios and both fitted parameters

**Hard acceptance test (no discretion):**
- The best-fit run must achieve \(E_{DM}/E_B \in [4.5, 6.5]\).
- If it cannot, Stage 2 is considered FAILED and you must stop and report “Stage 1 selection map is inconsistent with simple GRF spectrum under fixed sampling model.”

**Stage-2 raw artifacts (must write, ungamable):**
- `raw_aggregates.npz` must include:
  - `n_s_grid` (shape `(71,)` for 0.5..4.0 step 0.05)
  - `s_sigma_grid` (shape `(25,)`)
  - `objective` (shape `(71,25)`)
  - `E_ratio_energy` (shape `(71,25)`; aggregated \(E_{DM}/E_B\))
  - `E_ratio_count` (shape `(71,25)`; aggregated \(N_{DM}/N_B\))
  - `E_B` and `E_DM` (shape `(71,25)`)
  - `N_B` and `N_DM` (shape `(71,25)`)
  - `bootstrap_best_ns` (shape `(500,)`)
  - `bootstrap_best_s_sigma` (shape `(500,)`)
  - `bootstrap_yield_ratio_energy` (shape `(500,)`)

If any of these arrays are missing or have different shapes, the run is invalid.

---

## 4. Stage 3: The Graph Origin (Topology)

**Goal:** Determine if the **ICG / vacuum-relaxation** mechanism naturally produces the spectrum family found in Stage 2.

Stage 3 is not a “conceptual plausibility” step. It is a deterministic measurement on vacuum-relaxed Laplacians.

### 4.1 Protocol (idiot-proof, no choices)
Follow these steps exactly.

#### 4.1.1 Generate vacuum-relaxed graphs (fixed sweep)
Use `simulations/vacuum/vacuum_geometry_simulation.py` as the generator.

Run the following parameter sweep:
- `N = 200`
- `steps = 2000`
- `eta = 0.01`
- `tau_star = 1.0`
- `d_star = 3.0`
- `initial_graph = "ErdosRenyi_p0.3"`
- `temperature = 0.1`
- `alpha ∈ {0.5, 1.0, 2.0}`
- `beta ∈ {10.0, 30.0, 100.0}`
- seeds: `seed ∈ {0,1,2,3,4,5,6,7,8,9}`

Keep only runs that satisfy:
- `final_d_s ∈ [2.9, 3.1]`
- `fiedler_eigenvalue > 1e-6`

#### 4.1.2 Extract eigenvalues and define wavenumbers (fixed)
For each kept run:
- Compute Laplacian eigenvalues \(\lambda_0\le\lambda_1\le\cdots\le\lambda_{N-1}\).
- Discard the zero mode \(\lambda_0\).
- Define \(k_i=\sqrt{\lambda_i}\).

#### 4.1.3 Define the graph-implied phase spectrum (fixed OU model)
Model the phase as an Ornstein–Uhlenbeck process on the graph:
\[
\frac{d\phi}{dt} = -L\phi + \sqrt{2}\,\xi(t),
\]
which yields stationary mode variance \(\mathrm{Var}(\phi_i)\propto 1/\lambda_i\).

Define the effective spectrum for comparison to Stage 2 as:
\[
\rho(k) := \text{histogram density of }k_i,\qquad
P_{\mathrm{eff}}(k) := \frac{\rho(k)}{k^2}.
\]

#### 4.1.4 Fixed fit band and fixed slope estimator
Let \(k\) be the sorted list of \(k_i\). Define:
- `i_lo = floor(0.2 * len(k))`
- `i_hi = floor(0.8 * len(k))`
- Fit band is \([k[i_lo], k[i_hi]]\).

Compute \(\rho(k)\) by histogramming \(k_i\) into 30 equal-width bins on that band.
Compute \(P_{\mathrm{eff}}(k_b)=\rho(k_b)/k_b^2\) at bin centers.
Fit a line to \(\log P_{\mathrm{eff}}\) vs \(\log k\) by ordinary least squares. The fitted slope is \(-n_{s,\mathrm{graph}}\).

#### 4.1.5 Pass/fail criterion (fixed)
Let Stage 2 report \(n_s^\*\) and its 95% CI \([n_{s,\rm lo}, n_{s,\rm hi}]\).

- PASS if there exists at least one kept vacuum run with \(n_{s,\mathrm{graph}}\in[n_{s,\rm lo}, n_{s,\rm hi}]\).
- Otherwise FAIL.

**Stage-3 raw artifacts (must write, ungamable):**
- `raw_aggregates.npz` must include:
  - `sweep_alpha_grid` (shape `(3,)`)
  - `sweep_beta_grid` (shape `(3,)`)
  - `sweep_seeds` (shape `(10,)`)
  - `kept_mask` (shape `(3,3,10)`; boolean)
  - `final_d_s` (shape `(3,3,10)`)
  - `fiedler_eigenvalue` (shape `(3,3,10)`)
  - `n_s_graph` (shape `(3,3,10)`; NaN where not kept)
  - `k_fit_band_lo` and `k_fit_band_hi` (shape `(3,3,10)`; NaN where not kept)

If any of these arrays are missing or have different shapes, the run is invalid.

### 4.2 Stage-3 Implementation Checklist (for subagent)

- [ ] run the fixed sweep and keep only runs passing the filters
- [ ] compute \(n_{s,\mathrm{graph}}\) for each kept run using the OU model and fixed fit band
- [ ] write Stage 3 JSON exactly matching Section 6.3

### 4.3 Synthesis Result

- **If Match:** The ICG naturally explains the Dark Matter ratio.
- **If Mismatch:** The model requires modification (e.g., different graph topology).

---

## 5. Implementation Specs

### 5.1 Code Structure

- `stage1_induction.py`: Single-soliton SGLE solver. (High Resolution, small domain).
- `stage2_statistics.py`: Monte Carlo integrator. (No PDE solving, just probability distributions).
- `stage3_graph.py`: Vacuum-graph spectrum extraction + mapping to \(n_s\). (NumPy/SciPy; optionally reuse `vacuum_geometry_simulation.py`).

**Shared utilities (mandatory):**
- `common_grid.py`: FFT-based GRF generator + spectral differentiation.
- `common_topology.py`: winding/vortex detection utilities.
- `common_io.py`: JSON schema validation + reproducible RNG seeding.

**Change:** These are no longer “recommended.” They are mandatory. If you do not create shared utilities, you must still implement the exact same functionality in each stage and keep it byte-for-byte identical where applicable (e.g. GRF synthesis).

### 5.2 Key Metrics

- **Yield Ratio:** $R = E_{DM} / E_B$.
- **Selection Threshold:** The critical "Reynolds Number" for phase winding $\text{Re}_\phi \sim R \cdot \nabla \phi$.

**Reproducibility requirements (for subagent):**
- Every stage must accept `--seed` and write it into JSON output.
- Every stage must write the *exact* parameterization of spectra/forcing used (cutoffs, slopes, intermittency flags).

---

## 6. Output Interface

### 6.1 Stage 1 JSON (Selection Map)

```json
{
  "stage": 1,
  "map_resolution": [30, 30],
  "param_R": [0.5, 8.0],
  "param_sigma": [0.0, 3.0],
  "transition_matrix": "[[0,0,1...], ...]",
  "transition_matrix_ci_low": "[[...]]",
  "transition_matrix_ci_high": "[[...]]",
  "solver": {
    "dimension": 2,
    "L": 80.0,
    "N": 512,
    "dt": 0.0005,
    "t_relax": 20.0,
    "t_drive": 80.0,
    "forcing": {
      "k_vec": [0.0, 0.0],
      "boundary_noise_powerlaw_q": 1.0,
      "boundary_noise_mmax": 64,
      "boundary_noise_calibration": "two-pass: phase-only with w=1 for T_cal=5, then rescale coefficients to hit sigma_phi_target"
    }
  },
  "trials_per_point": 200
}
```

### 6.2 Stage 2 JSON (The Fit)

```json
{
  "stage": 2,
  "spectrum_family": "P_phi(k)=k^{-n_s}*1_[kmin,kmax], normalized so RMS(|∇phi|)=1 then scaled by s_sigma",
  "best_fit": {
    "n_s": 2.4,
    "s_sigma": 1.0,
    "k_min": 0.07853981633974483,
    "k_max": 40.21238596594935
  },
  "yield_ratio_energy": 5.28,
  "yield_ratio_count": 4.7,
  "bootstrap": {
    "resamples": 500,
    "n_s_ci_95": [2.2, 2.6],
    "s_sigma_ci_95": [0.8, 1.2],
    "yield_ratio_energy_ci_95": [4.8, 6.1]
  },
  "protocol": {
    "domain_L": 80.0,
    "grid_N": 1024,
    "num_realizations": 200,
    "peak_threshold_sigma": 2.0,
    "min_peaks": 2000,
    "max_peaks": 5000,
    "radius_clip": [0.5, 8.0],
    "shell_width": 2.0,
    "energy_weight": "R^3",
    "n_s_grid": [0.5, 4.0, 0.05],
    "s_sigma_grid": [0.25, 4.0, 25, "logspace"]
  }
}
```

### 6.3 Stage 3 JSON (Graph Origin Check)

```json
{
  "stage": 3,
  "vacuum_graph_source": "vacuum_geometry_simulation.py (fixed sweep)",
  "graph_phase_model": "OU: dphi/dt = -L phi + sqrt(2) xi, so P_eff(k)=rho(k)/k^2",
  "sweep": {
    "N": 200,
    "steps": 2000,
    "eta": 0.01,
    "tau_star": 1.0,
    "d_star": 3.0,
    "initial_graph": "ErdosRenyi_p0.3",
    "temperature": 0.1,
    "alpha_grid": [0.5, 1.0, 2.0],
    "beta_grid": [10.0, 30.0, 100.0],
    "seeds": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    "keep_filters": {
      "final_d_s": [2.9, 3.1],
      "fiedler_eigenvalue_min": 1e-6
    }
  },
  "vacuum_runs": [
    {
      "alpha": 1.0,
      "T": 0.1,
      "beta": 10.0,
      "d_s": 3.02,
      "fiedler_eigenvalue": 0.12,
      "n_s_graph": 2.35,
      "k_fit_band": [0.2, 2.0]
    }
  ],
  "stage2_target": {
    "n_s_ci_95": [2.2, 2.6]
  },
  "match": true
}
```

---

## 7. Synthesis Role

This 3-stage pipeline provides a **rigorous, falsifiable path** to cosmology. It separates the **Universal Physics** (how solitons react to noise) from the **Contingent History** (what noise distribution existed).
