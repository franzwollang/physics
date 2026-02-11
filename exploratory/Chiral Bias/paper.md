# Chiral Bias as a Local Window Phenomenon with Global Net Zero

## 0. Abstract

We argue that chirality in the soliton–noise framework is a window‑local effect built atop globally neutral foundations. Because rotor‑ladder degrees of freedom exist empirically, a chiral (pseudoscalar) order parameter is available. Globally, the net parity‑odd content across the infinite fractal hierarchy vanishes by (i) conservation/topological telescoping of total‑derivative densities, (ii) convergence via scale dilution and axion‑like relaxation, and (iii) isotropy (no cosmic handedness). Locally (per window/horizon), sliding‑window quenches generically trigger spontaneous symmetry breaking, Kibble–Zurek seeding, and rapid annealing to a single sign, with a tiny θ/μ5 tilt selecting the sign. Bias patches perturb neighbors through interfaces/defects; interface energy minimization favors alternation (often opposite signs across boundaries), while same‑sign matches occur transiently and leave predictable, small anomalies. Repeated window changes over cosmic history imply such patches recur with nonzero frequency, yet the global sum remains zero.

## Introduction

- Premise: rotor‑ladder degrees of freedom (U(1)→SU(2)→SU(3)) are realized in nature, so a chiral/pseudoscalar order parameter exists at low energy.
- Global principle: net parity‑odd content over the full hierarchy is zero (telescoping boundary terms, convergence of scale‑diluted contributions, and isotropy).
- Local mechanism: sliding‑window cooling drives SSB → KZ seeding → Allen–Cahn annealing to a single sign per patch, with a tiny θ/μ5 tilt selecting the sign.
- Patch interactions: biased regions perturb neighbors via boundary energies and defect networks; alternation (opposite signs) is favored but not guaranteed.
- Recurrence: because window changes repeat across epochs/scales, bias patches should arise with nonzero frequency, while global neutrality persists.

## 1. Global Principle: Net Zero Over All Scales

- Topological telescoping: the parity-odd densities are total derivatives (e.g., `F·F̃ = ∂_μ K^μ`, `Tr(W·W̃)=∂_μ K_W^μ`, `Tr(G·G̃)=∂_μ K_G^μ`, `R·R̃=∂_μ K_GW^μ`). Integrated over nested domains they telescope to boundary terms; with no preferred orientation at infinity, the outer boundary vanishes ⇒ net integral = 0.
- Convergence: axion-like relaxation and scale dilution make contributions decay with scale; the infinite sum converges. Random sign selection across domains drives the measure-weighted mean to zero.
- Isotropy: a nonzero global chirality would pick a cosmic handedness; net zero preserves isotropy while allowing strong uniform signs inside horizon patches.

## 2. Local Mechanism: Statistical Genesis and Topological Lock-in

The emergence of a local chiral bias arises inevitably from the statistics of a finite observational window acting on a globally neutral background.

### 2.1 The Order Parameter: Topological Helicity
We identify the chiral order parameter $m(x)$ with the coarse-grained topological density of the internal rotor sector (e.g., the SU(2) winding density):
$$ m(x) \;\sim\; \langle \mathrm{Tr}(F_{\mu\nu}\tilde{F}^{\mu\nu}) \rangle_{\text{window}} $$
This field is a Lorentz-invariant pseudoscalar. Globally, the vacuum is parity-symmetric ($\langle m \rangle_{\text{universe}} = 0$).

### 2.2 The "Tilt" from Sample Variance
In the early high-noise epoch, $m(x)$ fluctuates with mean zero. However, for any finite causal horizon with volume $V$, the spatial average possesses a non-zero statistical variance (Central Limit Theorem):
$$ \Theta_{\text{local}} = \frac{1}{V} \int_V m(x)\,d^3x \;\sim\; \pm \frac{1}{\sqrt{N_{\text{domains}}}} $$
This fluctuation $\Theta_{\text{local}}$ acts as a scalar chemical potential. While the averaging occurs in the cosmic rest frame, the resulting condensate is a spacetime scalar, preserving Lorentz invariance while spontaneously breaking Parity.

### 2.3 Spontaneous Symmetry Breaking and Domain Collapse
As the noise temperature drops, the potential for $m$ develops degenerate minima. The statistical tilt $\Theta_{\text{local}}$ breaks this degeneracy, lowering the energy of one sign relative to the other.
*   **Without Tilt:** The universe forms a persistent 50/50 patchwork of domains.
*   **With Tilt:** The energy difference creates a pressure on domain walls, driving the rapid expansion of the favored vacuum. This "False Vacuum" collapse anneals the entire horizon to a single chiral sign ($+m_*$ or $-m_*$) exponentially faster than random coarsening.

## 3. Local Effective Theory (Window Level)

### 3.1 Why SU(2) is the Chiral Receiver
While the vacuum bias $\Theta_{\text{local}}$ is a universal pseudoscalar, it does not affect all forces equally. We observe maximal parity violation only in the Weak (SU(2)) sector. The framework explains this via topological susceptibility:

*   **U(1) (Vector-Like):** The topological density $\mathbf{E}\cdot\mathbf{B}$ is a total derivative. Unless monopoles are active, the vacuum angle $\theta_{\text{EM}}$ has no effect on the equations of motion. U(1) remains parity-symmetric.
*   **SU(3) (Rigid/Screened):** The Strong sector forms "volume-lock" baryons (tripods) that mechanically resist chiral twisting. Furthermore, relaxation mechanisms (effective axions) in the high-stiffness amplitude sector efficiently screen any $\theta_{\text{QCD}}$, restoring Parity to the strong force (solving the Strong CP problem).
*   **SU(2) (The Goldilocks Zone):** The Weak sector governs surface orientations (rotors). It is topologically distinct enough to support chiral textures (unlike U(1)) but "soft" enough not to mechanically screen the bias (unlike SU(3)). The vacuum tilt $\Theta_{\text{local}}$ thus latches exclusively onto the SU(2) sector, energetically filtering out "Right-Handed" textures and leaving the observed "Left-Handed" weak force.

### 3.2 Effective Lagrangian Terms
At the window level, the bias manifests as tiny, dimensionless θ parameters:
- Weak Sector (Active): `(θ_W/4)Tr(W·W̃)` — The dominant source of PV.
- EM Sector (Inactive): `(θ_EM/4)F·F̃` — Physical effects suppressed/vanish.
- Strong Sector (Screened): `(θ_s/4)Tr(G·G̃)` — Relaxed to zero dynamically.

### 3.3 Axial Background
An effective axial background `J_5^μ B_{5μ}` (with `B_{50}≡μ_5`) may be present transiently during domain alignment; it damps with the bath (`T_eff(τ)`). Relaxation requires an axion-like rotor `a/f` to relax `θ_s` (EDM bounds); EM/weak θ can remain tiny. θ’s are dimensionless local invariants under conformal rescaling.

## 4. Putting It Together: Local Bias with Global Neutrality

- Patch uniformity: within a window/horizon, SSB + annealing selects one chirality; observables (e.g., parity-odd correlations, birefringence) reflect the local θ.
- Inter-patch cancellation: across the hierarchy, signs flip between domains/levels; parity-odd integrals telescope to zero; magnitude decays with scale ⇒ convergence.
- Conformal consistency: θ’s are dimensionless and persist under window rescalings while global neutrality is maintained.

### 4.1 Black-Hole Shell Anisotropy and Detectability

- Constant-power preserved: parity matching/mismatching at the inner surface redistributes energy angulary and temporally, but the integrated radiative flux remains ≈ constant (window-exhaustion argument unchanged). Global R(t) shrinkage is monotonic; anisotropy modulates only sectoral dR/dt(θ,φ).
- Transient anisotropy: same-parity inflow suppresses annihilation at the boundary ⇒ temporary “thickening”/hot spots (larger |m| ordered regions). Opposite-parity inflow annihilates efficiently, keeping the wall thin. Timescales: t_micro ≪ t_smooth ≪ t_evap ⇒ quasi-spherical long-term evolution.
- Chiral relaxation: for strong-sector components, axion-like relaxation further damps sustained mismatches; EM/weak θ’s remain tiny.

Detectable signatures (if same-parity feeding occurs)

- Horizon-scale polarization (EHT-class):
  - Enhanced circular polarization (CP) near the photon ring; azimuthal CP/EVPA asymmetries persisting longer than GRMHD baselines.
  - Parity-odd mode mixing/birefringence patterns localized to hot sectors; multi-epoch variability correlated with accretion episodes.
- Jet/inner-disk correlations:
  - Temporary shifts in jet launching efficiency/handedness; asymmetric Faraday conversion and RM near the jet base.
  - Fewer CP sign flips along the first few parsecs during same-parity episodes.
- Light-curve/variability:
  - Short-lived flares or QPO phase/timing tweaks as excess ordered energy is dissipated/unwound; relaxation back to baseline once inflow parity decorrelates.
- Population-level test:
  - Bimodal, patchy helicity distributions with enhanced sign flips across BH interfaces (bulk vs core), but rare systems showing aligned signs and the anomalies above.

## 5. Open Derivations

1. Rotor Landau coefficients: derive `α(T_eff)` from integrating out fast rotor modes (dependence on `(κ, β_pot, γ, η_0)`), fix `m_*`.
2. KZ scaling with windows: compute freeze-out length for sliding-window quenches; compare to defect statistics.
3. Telescoping/convergence: quantify decay and cancellations of parity-odd integrals across a synthetic hierarchical network.
4. Axion-like rotor: build `V(a)` from rotor bundle energetics; match strong-CP and EDM bounds.

## 6. Summary

Global net-zero chirality is enforced by conservation, convergence, and isotropy. Local chiral bias arises naturally from sliding-window dynamics that trigger SSB → seeding → annealing, yielding a single sign per horizon patch with tiny θ and/or transient `μ_5`. Local asymmetries and global neutrality coexist, predicting patchy, scale-local signatures with vanishing global handedness.
