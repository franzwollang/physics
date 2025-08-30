# Infinite‑Clique Graph — One‑Page Claims Sheet (internal)

## Core claims (what is established)

- Common origin of entanglement and projective‑latching measurement: the same phase‑locking free‑energy (inertia vs coupling + weak stabilizer) yields (i) a symmetric droplet–droplet double‑well (entanglement) and (ii) an asymmetric system–apparatus double‑well with decoherence and barrier‑induced latching (measurement).
- Entanglement as phase‑link hysteresis: metastable double‑well between droplets, Kramers‑type lifetime, nonlocal correlations without signalling; GR‑compatible (no stress–energy flux; waves propagate at `c_s`).
- Single‑mechanism continuity of measurement: tuning the same landscape parameters `(B, K, λ_A, ζ, T_eff)` smoothly deforms continuous→POVM‑like→projective+latching; primary switch at `B = K/4`; remaining thresholds packaged as conformal‑invariant ratios `(α, β, Ξ1–Ξ3)`.
- Effective collapse without a postulate: Born statistics `p ∝ |ψ|^2` arise from drift–diffusion guidance `p_* ∝ e^{Φ/Θ} ≈ I(x)`; einselection (dephasing) + Kramers‑latching yield durable records and Lüders‑like updates, so projective outcomes are dynamical, not axiomatic.
- Emergent complex state from a real field: analytic‑signal construction justifies `Ψ = A + i(1/ω_0)∂_t A` at grain scale; amplitude/phase `(w, φ)` well‑defined under narrowband/coarse‑grain assumptions.
- Volume‑normalized couplings tame infinities: grain‑internal weights `J_{ij}=η_0/V_K`, `J'_{ij}=J'_0/V_K` → extensive energies and finite continuum limits.
- On‑site stiffness emerges from pairwise cohesion: `γ_K = (η_0/2) V_K` (discrete); continuum pair density analogue `γ ∝ V_K^2` (appendix explains discrete vs continuum scaling).
- Lorentzian phase sector: hyperbolic quadratic action with single cone `c_s^2 = κ w_*^2` (with `κ ≡ J'_0`), no preferred frame on near‑regular subgraphs.
- Newtonian 1/r window: clustered screened sources with Yukawa tails produce a monopole far‑field `Φ(r) ∝ 1/r` for `R_cl ≪ r ≪ ℓ`.
- Spectral‑dimension fixed point: minimizing `F(d)=b_E e^{−α_dim d}+β_dim d^2` yields `d_* = (1/α_dim) W(α_dim^2 b_E/2β_dim)` (small‑integer bias).
  - Robustness: replacing the exponential cost by power‑law/log forms still yields a unique stable minimum (analytic scalings); the exponential ansatz gives the compact Lambert‑W expression.

## Assumptions and scope

- Substrate: uncountable nodes on an infinite clique; grains are countable coarse‑grains with spectral‑volume `V_K ∝ R_K^{D_s}`.
- Narrowband/separation of scales at grain level; coarse‑grain windows define numerical couplings, but invariants (e.g., `ħ_eff`) persist.
- Two channels: amplitude (cohesion) depends on `|Ψ|^2` at both ends; phase (correlation) composition‑independent.

## Key relations

- Equilibrium amplitude: `w_* = √(β_pot/(2γ))`.
- Signal speed: `c_s^2 = κ w_*^2`, `κ ≡ J'_0 > 0`.
- Phase inertia: `M_p ∝ V_K R_K^2`.
- Effective action (invariant): `ħ_eff ≃ E_cell · τ_cell`.

## Where established (paper refs)

- Complex field emergence: §3 (Coarse‑Graining, Analytic Signal).
- Volume‑normalized couplings and `γ`: §4, §7; Appendix (discrete vs continuum scaling).
- Lorentzian action and cone: §8 (Hyperbolic Action and Lorentz Cone).
- Entanglement as phase‑link hysteresis: §10.
- Entanglement → Measurement: §11 (From Entanglement to Measurement) with derivation in Appendix (Measurement Double‑Well and Pointer States) and Appendix (Threshold Criteria for Measurement Regimes).
- Newtonian window: §7 (Helmholtz/Yukawa → 1/r limit).
- Spectral‑dimension fixed point: §6 (vacuum relaxation) and Appendix D.

## Deferred/explicitly not done here

- Curvature sector (dynamical metric) derived from phase operator; left to future work.
- Full CHSH/no‑signalling with stochastic measurement back‑action; framework outlined, proofs deferred.
- Broad simulation calibration of `(α_dim, β_dim)` and dispersion/isotropy maps.

## Validation/diagnostics (minimal targets)

- Dispersion: measure `ω^2 ≈ c_s^2 k^2` on near‑regular subgraphs; isotropy of `c_s`.
- Newtonian window: embed compact sources, recover `Φ∝1/r` for `R_cl ≪ r ≪ ℓ`.
- Spectral dimension: heat‑kernel `K(t)` slope → `D_s`; fit `F(d)` near minimum.
- Local‑unit invariance test: search for tiny high‑resolution deviations/drifts in “constants” relative to external standards (thermodynamic convergence of co‑variation).
- Entanglement lifetime: verify Kramers scaling `τ ≈ (ω_0 ω_b/2π ζ) exp(ΔE_b/(k_B T_eff))` vs coupling `C`, inertia `M_p`, and environment `(ζ, T_eff)`.
- Regime mapping (single‑mechanism continuity): vary `α=B/K` and environment to cross `α=1/4`, `Ξ1≈1`, `Ξ2≈1`, `Ξ3≈1` boundaries; observe continuous→POVM→projective+latching deformation.
- Bell/no‑signalling checks: CHSH violation with invariant single‑wing statistics; all phase excitations limited by `c_s`.
