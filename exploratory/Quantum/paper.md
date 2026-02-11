# Interference, Measurement, and Entanglement in the Soliton‑Noise Framework: Simple Tests and Small Deviations (Draft)

## 0. Abstract

We develop simple, parameter‑light tests of quantum phenomena within the soliton‑noise framework and identify small, falsifiable deviations from standard predictions. Using only the phase‑channel dynamics of the foundations paper and a narrowband envelope approximation, we show:

- Schrödinger evolution for the complex envelope in local units (scope and conditions stated).
- Born statistics and effective “collapse” arise dynamically from drift–diffusion guidance (einselection + Kramers‑latching), without a collapse postulate.
- Interference: a concrete multi‑slit deviation—excess normalized variance beyond Poisson scales ≈ 1/k at bright fringes due to soliton radiative phase noise; in‑slit detectors yield a pattern not identical to (k−m) slits because disturbed paths coherently perturb the guidance field.
- Entanglement: a practical test using gravity/acceleration as a “noise knob” predicts symmetric environment dependence of lifetime and hysteresis near threshold (latching/unlatching), distinct from smooth Markov dephasing.
  We provide minimal protocols, micro‑to‑macro links for the phase‑noise parameter σ*ε in terms of (C, M_p, S*ξ(0)), and controls. The goal is not a full reconstruction of quantum theory, but focused, testable bridges: interference, measurement (continuous→projective+latching), entanglement, and small deviations suitable for near‑term experiments.

## 1. Motivation and scope

The foundational paper established that an emergent complex field Ψ arises from a pre‑geometric substrate and that its phase sector realises a `U(1)` symmetry—the complex phase rotor—providing a first‑principles basis for electromagnetic phenomena and a universal, composition‑independent phase channel. With this `U(1)` structure already on solid footing, the present work focuses on accessible quantum behaviours that directly involve EM‑like phase dynamics and can be tested with minimal additional assumptions.

Concretely, we isolate simple, high‑leverage phenomena and turn them into parameter‑light protocols with falsifiable deviations from standard expectations:

- Interference and measurement: show how self‑guidance by the phase channel produces interference, while the same system–apparatus potential in an asymmetric setup yields a continuous deformation from weak (POVM‑like) to projective, latched measurement. Born statistics and effective “collapse” emerge dynamically from drift–diffusion guidance plus einselection and Kramers‑latching, with no external collapse postulate.
- Entanglement as a physical link: model entanglement as a hysteretic, bistable phase‑link whose lifetime depends symmetrically on the combined environments of both parties, distinguishing it from purely local Markov dephasing.
- Testable deviations: propose two concrete, parameter‑tight tests: (1) a multi‑slit excess‑variance law that falls ≈ 1/k at bright fringes; (2) an entanglement experiment using gravity/acceleration as a practical “noise knob” to probe symmetric lifetime dependence and hysteresis near threshold.
- Connection to standard dynamics: state the Schrödinger envelope mapping and its validity conditions in local units, clarifying when familiar nonrelativistic dynamics are recovered between interactions.

What this is not: a full axiomatization or reconstruction. We focus on mechanisms already present in the phase channel that (i) admit compact derivations (envelope/SVEA, Fokker–Planck, Kramers), and (ii) lead to accessible experiments with few new parameters. Throughout we provide:

- Minimal translations and validity conditions (narrowband, weak inhomogeneity, paraxial/slow envelopes, local units).
- Micro‑to‑macro links for σ_ε via (C, M_p, \(S_{\xi}(0)\)) with conformal‑covariance caveats (dimensionless thresholds).
- Protocols and controls designed to distinguish framework‑specific predictions from standard expectations.

Scope note. This work focuses on U(1) phase dynamics and accessible quantum behavior. Non‑Abelian internal frames, covariant transport, and induced Yang–Mills curvature (SU(2), SU(3)) are deferred to follow‑up QFT papers.

## 2. Minimal model (phase rotor with weak stochasticity)

- Field at the detection screen (scalar, single polarization/window):
  I(x) = |Σ_{j=1}^k A_j(x) e^{i(φ_j(x) + ε_j)}|^2, with independent, zero‑mean phase perturbations ε_j from soliton radiative noise.
- Assumptions (leading order):
  - Equalized slits: |A_j(x)| = A(x)/√k (constant total flux across k).
  - ε\*j are small, stationary, mutually uncorrelated at the combining plane (correlations discussed in §6) with Var[ε_j] = σ_ε^2 (frequency‑band limited).
  - Detection integrates over a short time bin Δt within which ε statistics are approximately stationary.

Expand to O(ε^2). Let S*k(x) := Σ*{j=1}^k e^{iφ_j(x)} and C_k(x) := Σ*{j≠m} e^{i(φ_j(x) − φ_m(x))}.

- Mean intensity (to O(ε^2)):
  E[I] ≈ |A|^2 (|S*k|^2/k) · (1 − σ*ε^2/2) + |A|^2 · O(σ_ε^2/k).
- Variance (dominant stochastic term at fixed photon number mean):
  Var[I] ≈ |A|^4 · (σ_ε^2/k) · F_k(x) + higher orders,
  where F_k(x) is an order‑unity geometric factor built from {φ_j(x)} (maximal near fringe maxima where |S_k| is large).

Resulting normalized variance (beyond counting statistics):
Var[I]/E[I]^2 ≈ (σ*ε^2/k) · G_k(x) + O(σ*ε^4),
with G_k(x) ≈ F_k(x)/(|S_k|^4/k^2). For symmetric k‑slit arrays and well‑collimated sources, G_k(x) is O(1) near high‑visibility fringes.

Interpretation: the soliton radiative noise contribution averages down as ~1/k for fixed total flux, yielding a distinct k‑dependence not present in standard ideal QM (which predicts only Poisson 1/N_count plus technical noise).

## 2.4 Schrödinger envelope limit (translation and scope)

Under the narrowband/slowly varying envelope approximation (SVEA) used throughout, the complex envelope $\psi(x,t)=e^{-i\Omega_0 t}\,\Psi(x,t)$ of the phase field obeys an approximate Schrödinger equation in local units,
\begin{equation}
i\,\hbar*{\mathrm{eff}}\,\partial_t \psi(x,t) \;=\; -\frac{\hbar*{\mathrm{eff}}^2}{2\,m*{\mathrm{eff}}}\,\nabla^2 \psi(x,t) \; +\; V*{\mathrm{eff}}(x,t)\,\psi(x,t),
\end{equation}
with $m_{\mathrm{eff}}=\hbar_{\mathrm{eff}}\,\Omega_0/c_s^2$ and $V_{\mathrm{eff}}=\hbar_{\mathrm{eff}}\,\delta\omega$ as induced by weak index/potential variations. This follows from the phase-sector wave equation and the dispersion expansion around the carrier. For massive (subluminal) soliton cores, a multi-scale expansion around the internal mode yields the same nonrelativistic form up to renormalised prefactors. Validity requires narrowband excitation, weak inhomogeneity, amplitudes near $w_*$, and paraxial/small-angle propagation. Standard canonical commutation and the uncertainty relation for the envelope follow in this limit; see Appendix B.5.

## 2.5 Born statistics and effective collapse (no postulate)

The drift–diffusion guidance derived in §2.2 leads to Born statistics at leading order and to an effective collapse dynamics without an external postulate (see Appendix B.1–B.2 for the corresponding CPTP/Kraus map and Lindblad dephasing form):

- Stationary paths obey $p_*(x) \propto e^{\Phi(x)/\Theta} \approx I_k(x)$ when $\Phi/\Theta$ is small, with $I_k=|\sum_j A_j e^{i\phi_j}|^2\equiv|\psi|^2$. Calibrating $\Theta$ on the two-slit baseline fixes the proportionality.
- In the asymmetric system–apparatus potential, tuning $(\alpha,\beta,\Xi_1\!:\!\Xi_3)$ continuously deforms weak/POVM-like readout into projective+latching. Einselection ($\Gamma_{\rm dephase} T_{\rm int}\gtrsim1$) and Kramers-latching ($\Delta E_b/(k_B T_{\mathrm{eff}})\gg1$) yield durable records; conditionalisation on the selected pointer minimum reproduces Lüders-like updates.

Thus, in this framework, Schrödinger evolution governs the envelope between interactions, and the measurement transition is a dynamical outcome of the same phase-channel mechanism, not an added collapse axiom.

## 2.6 No‑signalling and delayed‑choice consistency (concise)

- No‑signalling (screen marginal): with joint state $|\Psi_{SA}\rangle$ and system POVM $M_x$, apparatus POVM $E_j$ (possibly chosen later), the screen marginal is $P(x)=\mathrm{Tr}[(M_x\otimes I)|\Psi_{SA}\rangle\langle\Psi_{SA}|]$, independent of $\{E_j\}$ and their settings or time of choice. Conditional fringes arise only after post‑selection on $j$.
- Delayed‑choice eraser: unconditional pattern remains classical; conditional patterns $I_\pm(x)=|c_1\psi_1\pm c_2\psi_2|^2$ are recovered by measuring the apparatus in superposition bases. Imperfect eraser yields partial recovery with visibility set by the basis rotation.
- Empirical assumption: pointer overlap $|\langle A_n|A_m\rangle|\approx\exp(-\mathcal A_{mn}/\hbar_{\mathrm{eff}})$ from action separation; in practice $\ll 1$ in projective regimes; finite overlap maps to partial visibility.

## 2.7 Casimir effect and zero‑temperature residuals (framework account)

Casimir forces arise here from how boundaries reshape the spectrum of phase‑sector normal modes and, consequently, the windowed ground‑state energy. The calculation is a difference‑of‑spectra in local units and needs no absolute vacuum energy.

- Mode picture (local units). In any bounded region the phase sector decomposes into harmonic modes with frequencies $\{\omega_n\}$ set by boundary conditions. The coarse‑grain invariant $E_{\mathrm{cell}}\,\tau_{\mathrm{cell}}\simeq\hbar_{\mathrm{eff}}$ and the harmonic normal form together give a zero‑point energy $E_{0,n}=\tfrac12\hbar_{\mathrm{eff}}\,\omega_n$ per mode. Only differences in spectral sums are observable; we therefore compute the change induced by inserting boundaries.
- Pressure as a spectral difference. For two parallel, high‑admittance plates a distance $a$ apart, the interior spectrum differs from the exterior continuum. The renormalised free energy $\mathcal F(a)$ is the difference of the two spectral sums evaluated with the same sliding window; the Casimir pressure is $P(a)=-\partial\mathcal F/\partial a$. For ideal reflectors this reproduces the standard magnitude with the replacements $\hbar\to\hbar_{\mathrm{eff}}$ and $c\to c_s$:
  \[ P*{\rm ideal}(a) \;=\; -\,\frac{\pi^2}{240}\,\frac{\hbar*{\mathrm{eff}}\,c_s}{a^4}. \]
- Real materials (Lifshitz mapping). Finite‑conductivity and dispersion enter through the same admittance map used for screens in §2.2: $\mathcal A(\mathbf x)$ defines reflection coefficients $r_{\rm TE/TM}(\omega,\mathbf k_\parallel)$. The force follows the Lifshitz‑type expression (Matsubara sum or contour integral) with $\hbar_{\mathrm{eff}}$ and $c_s$, reducing to the formula above in the perfect‑reflector, zero‑temperature limit.
- $T=0$ residuals. The nonzero $T=0$ force is the difference of ground‑state (zero‑point) energies of the phase modes in the two geometries, not an extractable bulk “vacuum reservoir.” Energy bookkeeping is ordinary: mechanical work done while changing $a$ equals the change in the windowed spectral energy; no net energy can be cyclically extracted once dissipation and path dependence are accounted for.
- Consistency with conformal covariance. Because predictions are stated in local units, homogeneous drifts of the window do not alter the dimensionless $a$‑dependence. Geometry, material response (via $\mathcal A$ or $r$), and temperature control the observable deviations; thermal corrections are included by the standard Matsubara sum with $k_B T$ in local units.

Connections to earlier sections: the same admittance picture used for interference with screens (§2.2) supplies the boundary conditions here; the propagation speed is the universal $c_s$ from the phase cone; and the appearance of $\hbar_{\mathrm{eff}}$ reflects the coarse‑grain invariant action scale of the phase sector.

## 2.8 Quantum eraser (spotlight): nonlocal phase‑link, time‑like interactions, and delayed‑choice

We summarise the quantum eraser in this framework, emphasising that all couplings are strictly local and time‑like at the devices, while the entangled pair shares a single nonlocal “phase‑link” degree of freedom that determines coincidence statistics without enabling signalling.

1. Source (pair creation; activation)

- A nonlinear source (e.g., SPDC) creates two phase‑rotor solitons (signal→Alice, idler→Bob) and kicks the shared relative‑phase coordinate into the metastable well of the two‑body potential
  \[ V*{\rm rel}(\theta) \approx -K\cos\theta,\quad K>0, \]
  equivalently realised by the system–apparatus Landau form
  \[ V*{\rm tot}(\theta*A,\theta_S;\theta*{\rm set}) = -K\cos(\theta*A-\theta_S)\; -\; B\cos\big(2(\theta_A-\theta*{\rm set})\big)\; +\; \lambda*A\,\cos^4(\theta_A-\theta*{\rm set}). \]
  The bistability threshold for a two‑outcome pointer is \(B*{\rm crit}=K/4\); for \(B>B*{\rm crit}\) the apparatus has two stable minima (projective+latching), for \(B\ll K\) it behaves as a weak/continuous meter.

1. Alice path optics (slits, free space)

- Refraction/propagation is local; unless a pointer couples and latches, the shared link remains delocalised.

1. Bob pre‑eraser optics (delays, basis setting)

- Local unitaries set the effective axis \(\theta\_{\rm set}\) for any subsequent pointer latching; no record yet.

1. Bob’s which‑way vs eraser module (local coupling, two regimes)

- Which‑way (distinguishable): Pointer states orthogonalise (double‑well formed, latched), so the overlap
  \[ \gamma\_{mn} := \langle A_n|A_m\rangle \approx 0 \]
  and Alice’s off‑diagonals are suppressed in each coincidence class.
- Eraser (indistinguishable): Coupling/optics maps alternatives to nearly the same pointer state so that
  \[ \gamma*{mn} \approx 1 \quad (\text{in the eraser basis}), \]
  recovering fringes within the corresponding post‑selected sub‑ensembles.
  Microscopically, pointer overlaps obey
  \[ |\langle A_n|A_m\rangle| \;\approx\; \exp\!\Big(-\tfrac{\mathcal A*{mn}}{\hbar*{\rm eff}}\Big),\qquad \mathcal A*{mn} \sim \pi\, \tfrac{\Delta p*A^{(mn)}}{\hbar*{\rm eff}} ,\qquad \Delta p_A \approx \tfrac{2C}{M_p}. \]

1. Alice detection (local latching)

- The screen’s detector cell couples via the same phase sector and latches a pointer minimum (projective in the large‑barrier regime). This is a local, time‑like interaction; no influence propagates to Bob outside the light cone.

1. Bob detection (local latching in the chosen basis)

- Bob’s detector latches in the basis set by his module; which‑way vs eraser determines \(\gamma\_{mn}\) and therefore the visibility in coincidence classes.

1. Coincidence sorting and no‑signalling

- The single‑wing marginal is independent of the distant choice:
  \[ P*A(x)\;=\; \mathrm{Tr}\big[(M_x\otimes \mathbb I)\,\rho*{SA}\big], \]
  so delayed‑choice settings on Bob’s side cannot modulate Alice’s marginal distribution. Conditional (post‑selected) patterns follow directly: with two alternatives and amplitudes \(\psi*1,\psi_2\), the eraser‑basis outcomes yield
  \[ I*\pm(x) \;=\; \big|c*1\,\psi_1(x) \pm c_2\,\psi_2(x)\big|^2, \]
  while the unconditional sum is classical. More generally, off‑diagonals are multiplied by \(\gamma*{mn}\):
  \[ \rho*{mn}' \;=\; \gamma*{mn}\,\rho\_{mn}\quad (m\neq n). \]

Frame‑dependent causal narratives, single mechanism

- The nonlocal “phase‑link” is a shared constraint variable (no stress–energy flux), so different inertial frames may order spacelike‑separated local choices differently (Alice‑first vs Bob‑first narratives). Coincidence statistics depend only on local couplings through \(\gamma\_{mn}\) and are frame‑invariant; dynamical excitations still propagate within the light cone at \(c_s\), precluding signalling.

Dimensional fixed point, pre‑geometric origin of the “weirdness”

- Substrate and vacuum relaxation. Fundamentally there is no metric or fixed dimensionality: the substrate is an all‑to‑all (infinite‑clique) connectivity. Locality and an effective dimension emerge when links relax by minimising a coarse functional that balances information/entropy against a small dimensional bias, e.g.
  \[ F(d) = b*E\,e^{-\alpha\,d} + \beta*{\rm dim}\,d^2, \]
  whose stationary point is
  \[ d*\* = \frac{1}{\alpha}\, W\!\Big( \frac{\alpha^2 b_E}{2\beta*{\rm dim}} \Big) \approx 3 \quad (\text{vacuum fixed point}). \]
- Emergent metric and locality. Near this fixed point the phase kernel becomes a Laplace–Beltrami operator on a near‑regular, low‑valence subgraph; waves propagate on a single light cone at \(c_s\). Local devices couple only along time‑like worldlines in this emergent geometry.
- Why nonlocal links persist. Sparse, targeted cross‑links in the phase sector between distant droplets (the phase‑link bundle) are energetically cheap and shift the spectral dimension only at subleading order \(\mathcal O(N\_{\rm link}/N)\). They therefore coexist with the 3D fixed point without upsetting locality. These cross‑links carry no stress–energy and so are nontraversable constraints, not signalling channels.
- Interpretation. The apparent “weirdness” of delayed‑choice then reflects remnants of the pre‑geometric, maximally connected substrate: the theory permits non‑signalling, system‑spanning constraints (the phase‑link) atop an emergent 3D metric background. This preserves Bell‑type nonlocal correlations while keeping all dynamical influences inside the emergent light cone.

## 2.1 Size–wavelength relation from time‑averaged stochastic emission

Consider a localized phase‑rotor soliton with characteristic width (attractor) σ. Let its internal coarse mode be a narrowband oscillator with central frequency ω*0 and correlation time τ_c (set by the internal relaxation and bath). The outward phase field obeys
\[ (\partial_t^2 - c_s^2 \nabla^2)\,\phi*{\rm out} = J*\sigma(t)\,\chi*\sigma(\mathbf r), \]
where χ*σ is the normalized spatial source (support ~ σ) and J*σ(t)=J*0\,[\cos(ω_0 t)+\xi(t)] combines the coherent mode and a zero‑mean stochastic drive ξ with spectrum S*ξ(ω) supported near ω_0 and bandwidth Δω\ll ω_0.

In the far field and over bins \( \Delta t \gg \tau_c \), the emitted spectrum is filtered by the soliton’s spatial form factor \(F_{\sigma}(\mathbf k)=\int \chi_{\sigma}(\mathbf r)e^{-i\mathbf k\cdot\mathbf r}\,d^3r\), selecting dominant wavenumber \( k_{\ast} \sim 1/\sigma \). Time‑averaging suppresses the stochastic sidebands, yielding a narrow line at (ω*0,k*_) with
\[ \langle |\Phi(\mathbf k,\omega)|^2 \rangle \propto |F*σ(\mathbf k)|^2\,\delta(\omega-ω_0) + O(Δ\omega). \]
Thus the observed wavelength scales with the soliton width, \( \lambda*\_ \sim 2\pi/k\_\_ \propto σ \). This holds for electron‑like (sub‑luminal) and photon‑like (luminal) cases; the latter requires the retarded construction in §2.3.

## 2.2 Self‑action guidance with stochasticity (Langevin and Born‑like limit)

Let the soliton center $\mathbf x*s(t)$ experience a weak backreaction from its own outward field scattered by apertures. In the paraxial/far‑field regime, the superposed field at the detection plane defines an intensity pattern $ I_k(\mathbf x) = |\sum*{j=1}^k A*j(\mathbf x) e^{i\phi_j(\mathbf x)}|^2 $. The phase‑channel coupling produces an effective guidance potential $ U*{\rm guide}(\mathbf x*s) = -\gamma_G\,(G*\ell * I*k)(\mathbf x_s) $, a local convolution with kernel range ℓ (set by the near‑field of the droplet). The center obeys an overdamped Langevin equation
\[ \dot{\mathbf x}\_s = \mu\,\nabla (G*\ell*I*k)(\mathbf x_s) + \sqrt{2D}\,\boldsymbol\eta(t), \quad \langle\eta_a(t)\eta_b(t')\rangle=\delta*{ab}\delta(t-t'), \]
where the noise D originates from soliton radiative stochasticity (same source as §2.1). The associated Fokker–Planck equation for the path density p(\mathbf x,t) is
\[ \partial*t p = -\nabla\cdot(\mu p\,\nabla \Phi) + D\,\nabla^2 p, \quad \Phi:=(G*\ell*I\*k). \]

Screen model (imperfect reflector/resonator). Let the slitted screen be characterized by a spatial map of transmissivity $\tau(\mathbf x)\in[0,1]$ (high inside slits, near zero in the bulk) and an additional slot‑resonator enhancement $Q*s(\mathbf x)\ge 0$ (nonzero within the slit cavities). Define an
effective phase‑channel admittance map
\[ \mathcal A(\mathbf x) := \tau(\mathbf x) + Q*s(\mathbf x). \]
The soliton experiences a lower effective potential where $\mathcal A$ is large. A minimal guidance potential that captures this is
\[ U*{\rm guide}^{(\rm scr)}(\mathbf x*s) = -\gamma_G\, (G\*\ell * \mathcal A)(\mathbf x*s). \]
When the outward wavefront is present, the available admittance is further modulated by the retarded field envelope, yielding
\[ U*{\rm guide}(\mathbf x*s,t) \;\approx\; -\gamma_G\, \big(G*\ell \_ [\mathcal A\,\overline I_k^{\rm (ret)}]\big)(\mathbf x*s,t), \]
which reduces to the intensity‑based form above when $\mathcal A$ is approximately constant over the illuminated region. Through‑slit regions ($\tau$ large, $Q_s$ moderate) produce valleys in $U*{\rm guide}$; bulk screen ($\tau\approx 0$) produces ridges, so the soliton drifts preferentially through the slits.

Slit‑width/spacing conditions. For robust guidance: (i) slit width $w$ supports at least one transverse mode, e.g.
\[ w \;\gtrsim\; \frac{\lambda*\*}{2 n*{\rm eff}} \]
(with $n*{\rm eff}$ the effective index in the slot), (ii) slit spacing $d$ is not much larger than a few $\lambda*\*$ so that the admittance valleys overlap into a connected channel network, and (iii) $Q_s$ is moderate (too high traps; too low reduces contrast).

With Einstein relation $ D=\mu\,\Theta $ (Θ an effective local unit “temperature” of the phase bath), the stationary solution is
\[ p*\_(\mathbf x) \propto e^{\Phi(\mathbf x)/\Theta}. \]
For weak convolution (ℓ\to 0) and small $ \Phi/\Theta $, expanding the exponential yields $ p\_* \propto I\*k $ at leading order; more generally, choosing $ \Theta $ via calibration on two‑slit data fixes the map $ p\** = f(I*k) $ and recovers Born‑like statistics in the diffraction regime. The stochastic term ensures exploration and prevents trapping in spurious local extrema; the drift aligns with “paths of least resistance” (intensity channels) provided slit spacing is within the guiding kernel range (d \lesssim few \lambda\*\_).

## 2.3 Photon‑like soliton limit (co‑moving wavefront)

For luminal solitons, the outward wavefront co‑moves at $ v\approx c*s $, so guidance must be retarded. Let $ t_r=t-\tau_r $ with $ \tau_r \approx L/c_s $ the path delay to the nearest scattering plane of distance L. The effective potential becomes
\[ U*{\rm guide}^{\rm (ph)}(\mathbf x*s,t) = -\gamma_G\,(G*\ell \* I*k^{\rm (ret)})(\mathbf x_s,t), \quad I_k^{\rm (ret)}(\cdot,t):=I_k(\cdot; J*σ(t*r)). \]
When the internal correlation time obeys $ \tau_c\ll \tau_r \ll \tau*{\rm dyn} $ (soliton center varies slowly over $ \tau*r $), a WKB/eikonal averaging yields the same Langevin form as §2.2 with $ I_k $ replaced by the time‑averaged retarded intensity $ \overline{I}\_k^{\rm (ret)} $. Hence the transverse drift is governed by contemporaneous gradients of the accumulated phase intensity, ensuring consistent guidance even when the wavefront co‑moves with the soliton. In the limit of negligible delay (near‑field or cavity‑like geometries), the convolution kernel $ G*\ell $ regularizes the self‑interaction and preserves stability.
In practice, the slight effective retardation of the center near slit edges (due to guidance ridges) means the leading wavefront has already sampled both apertures; the resulting retarded interference field sets up a guidance network ahead of the center that steers it along interference-compatible paths toward the detector screen without violating $c_s$.

Lead–lag kinematics (causal field lead near slits). Close to a partially reflecting screen, the soliton center experiences an effective slow‑down due to the guidance potential ridges, with instantaneous center speed $v*{\rm eff}(t)\le c_s$. The radiated wavefront still propagates at $c_s$ along admissible optical paths that do not suffer the center’s barrier detours. Over a local interaction time $\tau*{\rm loc}$, the wavefront acquires a causal lead distance
\[ \ell*{\rm lead} \;=\; (c_s - v*{\rm eff})\,\tau*{\rm loc} \;=\; \epsilon*{\rm lead}\, c*s\,\tau*{\rm loc},\qquad \epsilon*{\rm lead}:=1-\frac{v*{\rm eff}}{c*s}\in[0,1]. \]
If the geometric detour the center would incur by skimming the screen is $\Delta L*{\rm detour}$, then whenever $\ell*{\rm lead} \gtrsim \Delta L*{\rm detour}$ the retarded intensity reaching the center’s vicinity is already biased toward the slit cavities (shorter optical paths with higher admittance), producing a net drift through the slits. This respects $c*s$ strictly (no superluminal influence); the field’s lead arises solely from the center’s effective slow‑down near barriers. In practice, $\Delta L*{\rm detour}\sim O(w)$ for slit width $w$, while $\ell*{\rm lead}$ is set by the local potential ridge strength and $\tau*{\rm loc}\sim L/c_s$, so narrow slits and moderate ridge strength favour guidance.

## 3. Benchmark: standard QM expectation

- Ideal quantum optics with k slits and fixed total flux predicts mean pattern I_QM(x) ∝ |Σ_j A_j e^{iφ_j}|^2 and variance dominated by counting statistics: Var[I]/E[I]^2 ≈ 1/N_count (Poisson) plus detector/technical noise.
- No additional path‑indexed stochastic phase term exists in the ideal theory; hence no systematic 1/k decay in excess variance beyond Poisson is expected solely from increasing k at fixed flux.

## 4. Predicted deviation and scaling summary

- Excess (beyond Poisson) normalized variance: R_ex(x; k) := Var[I]/E[I]^2 − 1/N_count.
- Framework prediction (leading order): R\*ex(x; k) ≈ (σ_ε^2/k) · G_k(x), maximized at bright fringes.
- Two‑slit baseline: R*ex(x; 2) ≈ (σ_ε^2/2) · G*2(x) sets σ_ε; as k increases, R_ex should fall ~1/k if geometry and coherence are maintained.

Scope and comparison. The $1/k$ normalized excess‑variance prediction is a multi‑slit effect under fixed total flux and weak, approximately independent path‑phase noise. Crucially, a configuration with $k$ slits where $m$ slits host detectors is **not** equivalent (in this framework) to a pure $(k-m)$‑slit block: measured paths still contribute a disturbed wavefront that reshapes guidance (prefactor change), whereas blocks remove those paths entirely. Classical optics often treats these as equivalent; here they are measurably different (see §4.1).

## 4.1 In‑Slit Detection: A Second Testable Deviation

The framework provides a further, distinct prediction when detectors are placed at $m$ of the $k$ slits.

\paragraph{Standard QM prediction (collapse).} In the standard interpretation, a which‑way measurement collapses the state. If a particle is detected at one of the $m$ slits, it is absorbed and does not reach the screen. If it is not detected, its state is projected onto a clean superposition of the remaining $(k-m)$ slits. The resulting interference pattern for transmitted particles is therefore **identical** to that of a $(k-m)$‑slit experiment.

\paragraph{Framework prediction (physical interaction).} In this framework, detectors are physical objects that interact with both the soliton core and its extended guiding wavefront.

- **Core interaction:** A core–detector interaction produces a local registration appropriate to the species (e.g., absorption for photons; inelastic scattering/charge transfer/secondary emission for electrons), thereby removing that trial from the transmitted ensemble.
- **Wavefront interaction:** The portion of the guiding wavefront passing through the $m$ instrumented slits interacts with the detectors, which act as imperfect absorbers/scatterers. This part of the wave is attenuated and phase‑shifted, but not perfectly zeroed.

The guidance potential $U\_{\rm guide}$ is therefore shaped by the superposition of $(k-m)$ clean paths and $m$ weak, disturbed paths. The resulting intensity pattern, $I'\_k(x)$, will be subtly different from a pure $(k-m)$‑slit pattern. This constitutes a direct, falsifiable prediction: an experiment that compares the final pattern from a system with in‑slit detectors to one where those same slits are physically blocked should find a measurable difference in fringe visibility, shape, or peak positions. The excess variance scaling would also be modified from a pure $1/k$ or $1/(k-m)$ law.
In particular, the normalized excess‑variance law with $\sim 1/k$ scaling applies to the multi‑slit case with measured (disturbed) paths; replacing measured slits by blocks changes the prefactor and the pattern and is not equivalent in this framework.

## 4.2 Visibility–distinguishability tests without which‑way hardware

Design a tunable, weak apparatus that couples to the phase channel but leaves no macroscopic record (no latching), e.g., a short interaction region with adjustable coupling $K$ and bath temperature $T_{\mathrm{eff}}$. Measure fringe visibility $V(K,T_{\mathrm{eff}})$ and estimate distinguishability via an independent overlap proxy $F\approx|\gamma_{12}(K,T_{\mathrm{eff}})|$ from homodyne readout of the weak pointer mode. Test Englert: $V^2 + D^2 \le 1$, with $D\ge\sqrt{1-F^2}$. Framework mapping: $F\approx\exp(-\mathcal A/\hbar_{\mathrm{eff}})$, $\mathcal A\propto\Delta p_A(K)/\hbar_{\mathrm{eff}}$ with $\Delta p_A\approx 2C/M_p$.

Falsifiable features:

- Saturation for near‑pure pointers (low bath): $V\approx F$ and $V^2+D^2\approx 1$.
- Monotone $V\downarrow$ with increasing $K$ or $T_{\mathrm{eff}}$; quantitative fit to $\exp(-\mathcal A/\hbar_{\mathrm{eff}})$.
- Deviations from saturation diagnose mixedness (finite bath) consistent with fidelity bounds.

Unequal intensities (two paths): if path powers are $p_1, p_2$ (with $p_1+p_2=1$), the observable visibility obeys $V = 2\sqrt{p_1 p_2}\,|\gamma_{12}|$. Experiments should report $p_1,p_2$ and use this formula when paths are not balanced.

## 4.3 Reflectance‑dependent guidance and variance prefactor (falsifiable)

- Guidance with back‑reflection. The screen’s optical response enters the guidance potential through an admittance map that includes reflectance $R(\omega,k_\parallel)$ and its phase $r$ (TE/TM):
  \[ U*{\rm guide}(\mathbf x_s,t) \;=\; -\gamma_G\,\big( G_{\ell} [\,\mathcal A(\mathbf x;\omega,k_\parallel,R)\; \overline I*k^{\rm (ret)}(\mathbf x,t)\,] \big)(\mathbf x*s), \]
  where $\overline I_k^{\rm (ret)}$ includes the near‑field, back‑reflected component. Thus $R$ modifies the spatial gradients of $\Phi := (G*\ell * [\mathcal A\,\overline I_k^{\rm (ret)}])$ that drive drift.

- Prediction for the variance law. The $1/k$ scaling survives, but the geometric factor depends on reflectance and polarization:
  \[ \frac{\mathrm{Var}[I]}{\mathbb E[I]^2} \;\approx\; \frac{\sigma*\varepsilon^2}{k}\; G_k\big(x; R(\omega,k*\parallel),\,\text{pol}\big) \; +\; O(\sigma*\varepsilon^4). \]
  Reflectance changes the prefactor $G_k(\cdot;R)$ (and the calibration factor $\kappa*\phi$ in §9), not the $1/k$ exponent.

- Falsifiable signatures (with geometry fixed):

  - Tune coating reflectivity $R$; map $R_{\rm ex}(x;k)$ at bright fringes. Expect the $\sim 1/k$ slope to persist but a systematic, polarization‑ and wavelength‑dependent change in the prefactor and fringe shape via $G_k(x;R)$.
  - Non‑monotonic contrast vs $R$: increasing $R$ strengthens back‑reflection (steeper $\nabla\Phi$ → contrast $\uparrow$) up to a cavity‑like regime (effective $Q_s$ too high) where throughput drops and patterns distort.
  - Phase‑sensitive shifts: TE/TM phase upon reflection shifts fringe positions in coincidence‑resolved sub‑ensembles without affecting single‑wing marginals.

- Novelty. Classical EM already links reflectance to fringe envelopes, but the reflectance‑dependent prefactor in the normalized excess variance law together with the $1/k$ scaling arises from the soliton noise channel and is specific to this framework.

## 4.4 Massive vs. Massless Slit-Spacing Sensitivity (A Falsifiable Prediction)

The causal separation between a soliton's core and its radiated guidance field creates a mechanical distinction between massive and massless particles, leading to a novel prediction about their sensitivity to the
center-to-center slit spacing \(d\) (not just slit width).

- **Massive Solitons (Integrated Guidance):** A massive soliton at velocity \(v<c*s\) radiates a massless phase field that arrives at the screen with a lead time \(\Delta t*{\rm lead}=L/v - L/c*s\). During this interval the guidance potential \(U*{\rm guide}\) integrates the back-scattered field and forms lateral "valleys" toward the slits. The integrated lateral drift scale obeys, to leading order,
  \[ L*{\rm drift} \;\sim\; \mu\,|\nabla \Phi|\, \Delta t*{\rm lead}, \]
  with mobility \(\mu\) and guidance potential \(\Phi=(G*\ell \* [\mathcal A\,\overline I_k^{\rm (ret)}])\). When \(L*{\rm drift}\gtrsim d/2\) near the screen, the core can be funneled into the nearest slit even for relatively larger spacings.

- **Massless Solitons (Instantaneous Steering):** For massless solitons, \(\Delta t*{\rm lead}=0\). The lateral drift available just at the boundary is effectively limited to the local interaction time; hence \(L*{\rm drift}\) is much smaller. The trajectory is set by the instantaneous gradient at the edge, yielding a sharper loss of transmission/interference as \(d\) increases (valleys no longer overlap).

- **Falsifiable Prediction (spacing):** Holding wavelength \(\lambda\) and slit width \(w\) fixed (with \(w\gtrsim \lambda/(2n*{\rm eff})\) so a transverse mode is supported), the visibility and transmission as functions of spacing \(d\) will show a **softer fall-off** for massive solitons than for massless solitons. Equivalently, the critical spacing \(d*{\rm crit}\) beyond which contrast collapses (valleys no longer overlap) satisfies
  \[ d*{\rm crit}^{\,(\rm massive)} \;>\; d*{\rm crit}^{\,(\rm massless)} \quad \text{at fixed } (\lambda,w), \]
  because massive solitons benefit from an integrated look-ahead guidance (\(\Delta t\_{\rm lead}>0\)). This difference is **not** attributable to wavelength; it follows from distinct guidance dynamics.

- **Width vs spacing (scope):** Width still matters through the modal condition \(w\gtrsim \lambda/(2 n\_{\rm eff})\) (cf. §2.2), and both species will show reduced transmission below that threshold. However, at fixed \(w\) above threshold, the **clean discriminator** is the spacing dependence: massive solitons retain interference/throughput to larger \(d\) than massless solitons.

## 5. Experimental protocol (minimal)

- Geometry: equal‑width slits, center‑to‑center spacing d, far‑field Fraunhofer screen. Fixed source flux and spectrum; fixed integration bin Δt.
- Procedure: measure time‑binned intensity at several fixed pixels x spanning bright/dark fringe regions while ramping k = 2→3→…→K; keep total transmitted power constant (attenuate per slit to A/√k).
- Estimation: at each k and pixel x, estimate mean and variance across time bins; subtract Poisson 1/N_count using calibrated detector QE to form R_ex(x; k).
- Fit: R*ex(x; k) ≈ (σ_ε^2/k) · G*k(x) with a single σ_ε estimated from k=2 and geometry‑computed G_k(x); test 1/k law.

Note: For environment‑driven tests (gravity/acceleration as a practical “noise knob”), see Section 10 for a dedicated treatment in the entanglement context and how the same mapping applies to multi‑slit excess‑variance measurements.

## 6. Controls and confounds

- Partial coherence: increasing k can change effective coherence if slits/sample are not rigorously stabilized. Control by verifying |S_k| visibility and spectral bandwidth.
- Path correlations: common‑mode phase noise across slits (e.g., source jitter) does not average down with k; design to suppress common‑mode (e.g., independent slit‑proximal phase shifters with matched statistics) or model a correlation coefficient ρ, yielding R_ex ∝ (1−ρ)/k + ρ (plateau test).
- Detector nonlinearity/afterpulsing: characterize and subtract with dark and bright calibrations.
- Diffraction envelope changes: maintain fixed per‑slit aperture; adjust per‑slit attenuation to keep total flux constant.

## 7. Parameter extraction and falsifiability

- From k=2 data at bright fringes, estimate σ*ε via R_ex ≈ (σ*ε^2/2) · G_2.
- Predict R_ex at higher k with no new parameters; failure to observe ~1/k scaling (after controls) falsifies the proposed stochastic term.
- Spatial diagnostic: R_ex peaks at bright fringes (large |S_k|) and is suppressed at near‑dark regions; map vs x for additional confirmation.

## 8. Framework link (local units and invariants)

- ε_j arise from low‑frequency phase‑current fluctuations (soliton radiative noise) in the phase channel; their local‑unit variance is controlled by the same environment that fixes ħ_eff.
- The 1/k scaling is expressed in dimensionless, local‑unit quantities (normalized variance); slow thermodynamic drifts in σ_ε are allowed but should be subleading within a run.

## 9. Micro‑to‑macro link: estimating σ*ε from (C, M_p, S*ξ(0))

To close the loop with the entanglement/measurement derivations, we express the path‑phase variance per time bin, σ_ε^2, in terms of microscopic parameters already used in the main paper.

Starting point (main paper, measurement/entanglement sections): the dephasing rate between two macroscopic pointer states is
\[ \Gamma\*{\rm dephase} \;\approx\; \frac{(\Delta p*A)^2}{2\,\hbar*{\mathrm{eff}}\,}\, S*\xi(0), \qquad \Delta p_A \;\approx\; \frac{2C}{M_p}. \]
Here C is the phase‑channel coupling, M_p the phase inertia of the droplet (or effective subregion), and S*ξ(0) the low‑frequency bath power.

In the multi‑slit setting, the phase on a given path j acquires small stochastic kicks ε*j(t). For short integration bins Δt and in the white‑noise limit (bandwidth ≲ 1/Δt), the phase variance grows linearly with time under a Markovian approximation. Introducing a geometric/projection factor κ*φ ∈ O(1) that maps the two‑pointer dephasing channel to a single‑path phase channel, we obtain
\[ \sigma*ε^2(\Delta t) \;\approx\; 2\,\Gamma*φ\,\Delta t, \qquad \Gamma*φ := \kappa*φ\,\Gamma\*{\rm dephase}. \]
Combining gives the practical estimate
\[ \boxed{\; \sigma*ε^2(\Delta t) \;\approx\; \kappa*φ\,\frac{(2C/M*p)^2}{\hbar*{\mathrm{eff}}\,}\, S\_\xi(0)\,\Delta t \; } \]
up to O(1) geometry and filtering constants. This provides the sought micro‑to‑macro link:

- Stronger coupling C, smaller inertia M*p, or a “hotter” bath S*ξ(0) increase σ_ε and thus the excess variance contribution R_ex.
- Larger ħ_eff (in local units) suppresses phase diffusion per unit time.

Remarks and usage:

- κ_φ captures slit geometry, path projection, and finite‑band effects; it can be calibrated from two‑slit data at bright fringes.
- The same (C, M*p, S*ξ(0)) appear in the entanglement lifetime and measurement latching formulas; using those independently estimated values makes the multi‑slit variance prediction parameter‑tight.
- Local‑unit invariance: the combination (C^2/M*p^2)(S*ξ(0)/ħ*eff) is dimensionless in local units, so σ*ε^2 scales with Δt as expected and is robust to conformal co‑variation.

With this relation, the predicted 1/k excess‑variance scaling (§4) can be written without new free parameters beyond κ_φ and known geometric factors G_k(x).

## 10. Entanglement Deviations: Testing the Physical Link Model

### 10.1 Foundational Model Summary and Objective

The soliton-noise framework models entanglement not as an abstract correlation but as a physical, metastable **phase-link** between two solitons (droplets). As detailed in the foundations paper, this link emerges from the same phase-channel dynamics that govern measurement, resulting in a double-well potential in the joint phase space of the two particles. This physical link has two key, testable properties:

1. **Metastability and Hysteresis:** The link is a bistable system with a finite energy barrier, `ΔE_b`, separating the entangled (correlated) state from the uncorrelated ground state. Its lifetime, `τ_link`, follows a Kramers-like escape law, `τ_link ∝ exp[ΔE_b / (k_B T_eff)]`, making it exquisitely sensitive to the effective temperature `T_eff` of its environment near the stability threshold. This implies hysteretic behavior (latching and unlatching) when the relevant dimensionless threshold (e.g., `ΔE_b/(k_B T_eff)`) is driven across unity; under exact conformal co‑variation this threshold remains fixed and no change is observed.
2. **Symmetric Environmental Dependence:** The phase-link is a single, non-local object. Its stability and lifetime depend on the *combined* noise environment of both solitons. The effective temperature, `T_eff,link`, is a function of the baths at both locations, `T_eff,link ≈ f(T_A, T_B)`.

The objective of a minimal experiment is therefore to test these two core predictions: **(1)** look for a sharp, hysteretic drop-off in entanglement lifetime as environmental noise is increased, and **(2)** verify that the lifetime depends symmetrically on the noise injected at *both* locations.

### 10.2 Experimental Test Using Gravity/Acceleration as a "Noise Knob"

Directly engineering and injecting calibrated noise is difficult. However, the framework predicts that gravitational potential and acceleration can serve as practical, albeit subtle, knobs on the effective temperature of the phase bath.

- **Setup (Photonic Entanglement):**

  - Place two entangled photon receivers (A and B) at locations with a controllable difference in gravitational potential, `ΔΦ = Φ_A - Φ_B` (e.g., one at ground level, one on a mountain or on a high-altitude platform).
  - Alternatively, for a stronger effect, mount one station (e.g., B) on a precision centrifuge or rotation stage to induce a controlled acceleration, `a`.

- **Mapping to Effective Temperature:** The framework provides a qualitative mapping from these physical variables to the effective bath temperature in local units:

  - **Gravitational:** `T_eff(Φ) ≈ T_0 [1 + α_g (Φ/c_s^2)]`, where `α_g` is a coupling constant. The link temperature becomes `T_eff,link ≈ T_0 [1 + α_g (Φ_A + Φ_B)/(2c_s^2)]`.
  - **Acceleration:** The Unruh-like effect in the soliton's co-moving frame adds a term `T_U(a) = ħa / (2πc_s k_B)`. The total effective temperature at station B would be `T_B(a) ≈ T_0 + T_U(a)`.

- **Measurements and Predictions:**

  1. **Symmetric Dependence Test:**
      - Measure the entanglement lifetime `τ_link` (via heralded coincidence decay) for different configurations: (A fixed, B at high altitude), (B fixed, A at high altitude), and (both at high altitude).
      - **Prediction:** The change in `τ_link` should depend on the *sum* of the potentials `(Φ_A + Φ_B)`. A change at site A should have the same effect as an identical change at site B. This contrasts with models where decoherence is a purely local process.
  2. **Hysteresis Test:**
      - Starting with a stable link (low altitude/acceleration), slowly ramp up the potential difference `ΔΦ` (or acceleration `a`) until the entanglement visibility drops sharply (unlatching). Then, slowly ramp it back down.
      - **Prediction:** The visibility will not recover along the same path. It will re-latch at a lower value of `ΔΦ` (or `a`), exhibiting a clear hysteresis loop characteristic of a bistable system. Standard Markovian decoherence models predict a single, reversible curve.

Conformal‑covariance caveat (observability). In local units, many quantities co‑vary with the analysis window and background. The experimentally relevant drivers are the **dimensionless** thresholds (e.g., `Ξ_1, Ξ_2, Ξ_3`, particularly `ΔE_b/(k_B T_eff)` and `Γ_dephase T_int`). If conformal co‑variation is exact, these ratios remain constant under changes in `ΔΦ` or acceleration and no change is observed. Observable drifts and hysteresis arise from small **non‑covariant residuals** (thermodynamic convergence), slight **asymmetries** in the co‑variation at A vs B, or **retarded/non‑equilibrium** effects in the link. Accordingly, operation near threshold is recommended to amplify small departures via the exponential Kramers sensitivity.

**Practical Note:** The fractional shifts from gravity are extremely small (`ΔΦ/c_s^2 ≪ 1`), requiring long integration times and operation near the critical threshold where the exponential Kramers rate sensitivity can amplify the effect. The acceleration‑based test is likely more feasible.

This experimental design provides a direct, parameter-light test of the framework's core claims about the physical nature of the entanglement link.

### 10.3 CHSH and the Quantum Potential: From Tanh to Cosine

We model the shared degree of freedom as a bistable phase‑current orientation with barrier $\Delta E_b$ and coupling $K$. Measurement at A/B is a rotation of the apparatus setting axes $(\alpha, \beta)$; outcomes $A,B\in\{\pm1\}$ correspond to the two pointer minima.

#### The "Classical Phase" Deviation ($\tanh\Delta$)

To see why the standard stochastic model fails to reproduce the exact quantum correlation, consider the "classical limit" where the phase link is treated as a simple Langevin rotor in a joint potential, ignoring amplitude backreaction. The effective potential combines the link coupling with the deep latching potentials of the apparatuses:
\[
  V_{\text{tot}} \;\approx\; -K \cos(\theta_A - \theta_B) \;-\; h \cos(2(\theta_A - \alpha)) \;-\; h \cos(2(\theta_B - \beta)).
\]
In the deep-well limit ($h \to \infty$), the pointers are pinned to $\theta_A \in \{\alpha, \alpha+\pi\}$ and $\theta_B \in \{\beta, \beta+\pi\}$. The effective energy of a configuration $(A,B)$ becomes $E_{AB} \approx -K A B \cos(\alpha - \beta)$.

If the system samples this landscape according to simple emergent thermodynamics (Boltzmann weights $P \propto e^{-E/T_{\text{eff}}}$), the joint probability is $P(A,B) \propto \exp(J A B \cos\Delta)$ with $J=K/T_{\text{eff}}$. This yields a correlation function:
\[ E(\Delta) \;=\; \sum_{A,B} P(A,B) A B \;=\; \tanh(J \cos\Delta). \]
While this satisfies no-signalling and reproduces the qualitative symmetries (zeros at $\pi/2$, extrema at $0,\pi$), it is "squarer" than the quantum cosine. For large $J$ it approaches a step function (the classical correlation limit); for small $J$ it is linear in $\cos\Delta$.

#### The Quantum Potential Fix ($Q$)

The exact cosine dependence $E(\Delta)=-\cos\Delta$ required by QM is recovered by including the **amplitude-phase backreaction** (Quantum Potential) term, which is inherent to the framework's envelope dynamics but neglected in the simplified Langevin model. As shown in the envelope limit (§2.4), the phase $\phi$ and amplitude $R$ are coupled via the continuity and modified Hamilton-Jacobi equations:
\[ \partial_t S + (\nabla S)^2 + V + Q = 0, \qquad Q = -\frac{\nabla^2 R}{R}. \]
The "Quantum Potential" $Q$ forces the phase dynamics to track the curvature of the probability density. In the measurement context, rotating the apparatus changes the boundary conditions for the mode, altering the amplitude profile $R$. The resulting force $-\nabla Q$ modifies the phase equilibrium away from the simple thermal-rotor prediction.

Clarification (scope of the claim). The existence of $Q$ is *not* a special Bell assumption: it is the generic amplitude–phase coupling that appears whenever the effective envelope dynamics are Schrödinger-like and one writes $\psi=R e^{iS}$. What is nontrivial—and therefore the real technical target—is that the **apparatus setting must actually reshape the effective phase boundary-value problem** (so that $R$ changes in the right setting-dependent way). In this framework that physical \"setting → operator\" pathway is supplied by matter–kernel coupling: the apparatus' matter configuration perturbs the local phase kernel/boundary conditions, which in turn shapes $R$ and therefore activates the $Q$ backreaction consistently.

Thus, recovering the exact cosine $E(\Delta)=-\cos\Delta$ is framed as an engineered consistency target for the full envelope + apparatus coupling, not an automatic consequence of merely writing down $Q$.

#### Cross-constraint (unification check): why “mirrors are strong” but “gravity is weak”

This same matter–kernel coupling that underlies weak-field gravity also underlies apparatus boundary conditions. Writing the kernel perturbation in linear response,
\[
  \delta K(\omega,k)\;\sim\;\chi(\omega,k)\,\rho_m(\omega,k),
\]
gravity probes the static IR limit $\chi(0,0)$, while slits/mirrors/detectors probe finite-band response $\chi(\Omega_0,k\sim 1/\lambda)$ and (crucially) dissipative/locking channels (the $\operatorname{Im}\chi$ part).

This yields a qualitative cross-constraint: if strong apparatus behavior were produced purely by a huge *static* kernel jump $\delta c_s^2\propto \chi(0,0)\rho_m$, then ordinary dense materials that act as strong boundaries would also create enormous static index shifts—i.e. they would be gravitationally bizarre. They are not. The framework is therefore pushed toward a unified resolution where $\chi(0,0)$ is small (weak gravity), but $\chi(\omega,k)$ can be large in narrow bands and in dissipative channels (strong apparatus coupling), consistent with causal dispersion constraints (Kramers–Kronig-type logic).

#### Hysteresis and Symmetric Dependence

In this framework, the visibility $V$ (and thus the CHSH value $S=2\sqrt{2}V$) is physically set by pointer overlaps and latching stability. Consequently, $S$ inherits the predicted hysteresis and symmetric environment dependence described in §10.2, providing a signature of the physical link mechanism.

### 10.4 Consolidated falsifiable predictions (no repetition)

1. Multi‑slit excess variance scaling: R\*ex(x;k)≈(σ_ε^2/k)G_k(x); 1/k fall at bright fringes with fixed flux; failure (after controls) falsifies the stochastic phase‑noise channel.

2. In‑slit detector vs physical block: transmitted pattern with m instrumented slits differs measurably from a pure (k−m)‑slit pattern (visibility/shape/peak shifts). Quantify difference vs detector admittance; pure block acts as reference.

3. Visibility–distinguishability law: V^2 + D^2 ≤ 1 with D from independent pointer tomography; near‑pure pointers saturate, mixed pointers follow fidelity bound V≤F. Fit F≈exp(−A/ħ_eff) with A∝(2C/M_p)/ħ_eff.

4. Imperfect eraser: conditional visibilities follow the eraser‑basis rotation (θ,φ); unconditional remains classical. Partial recovery matches calculated overlaps in the rotated basis.

5. CHSH with hysteresis and symmetric environment dependence: near threshold, correlation visibility shows hysteresis under slow setting sweeps; symmetric dependence on baths at A and B; no‑signalling holds for single‑wing marginals.

6. Delayed‑choice: unconditional P(x) independent of later apparatus POVM; conditional fringes recovered by post‑selection only. Time‑ordering reversals leave P(x) unchanged within statistical error.

Each item is parameter‑tied: (C,M*p,S*ξ(0),ħ*eff), ΔE_b/λ, and geometric factors (G_k, basis rotations). No new free parameters beyond calibration constants (κ*φ, overlap‑to‑fidelity map) are introduced.

### 10.5 Atomic stability and size scaling (mechanical account)

In this framework, the “tendency to decay to lower energy” acquires a mechanical explanation in terms of entanglement (phase‑link) lifetimes of a many‑body network inside the atom. Configurations are stable when their internal phase‑link landscape has deep minima and when unlatching by the internal bath is rare; they are unstable when barriers are shallow and the bath is effective.

- Kramers control of lifetimes. As in §10.2–10.3, excited configurations have lifetimes
  \[ \tau\;\sim\; \frac{2\pi\,\zeta}{\omega*0\,\omega_b}\; e^{\Delta E_b/(k_B T*{\mathrm{eff}})}. \]
  Deeper collective wells (larger effective coupling network \(\{C*{ij}\}\), smaller inertia \(M_p\)) increase \(\Delta E_b\); hotter/denser internal baths (larger modal density, stronger radiative/electronic dissipation) increase \(T*{\mathrm{eff}}\) and \(\zeta\), shortening \(\tau\).
- Size/excitation trends (qualitative scalings).
  - Increasing Z (more electrons, core excitations) adds internal modes → \(T*{\mathrm{eff}}\) and \(\zeta\) rise; spectator channels act as internal pointers, raising \(\Gamma*{\rm dephase}\approx (\Delta p*A)^2 S*\xi(0)/(2\hbar\_{\rm eff})\).
  - Open shells and configuration mixing create many near‑degenerate channels → shallower effective barriers and faster unlatching (einselection within the atom).
  - Rydberg excitation enlarges orbitals → reduced overlap with the core lowers path‑coupling \(C\) and hence \(\Delta E*b\); blackbody/near‑surface noise raises \(T*{\mathrm{eff}}\).
- Mapping to standard decay channels (mechanical view).
  - Radiative (E1/M1/E2): link rearrangement dissipates excess free energy via the phase channel’s EM‑like radiation while the internal network relatches into a deeper minimum.
  - Auger/autoionization (heavy atoms): electron–electron coupling provides the bath that unlatches a core‑excited link, transferring energy to an outgoing electron (large internal \(S\_\xi(0)\)).
  - Internal conversion/intersystem crossing: spin–orbit–mediated mixing opens additional dissipative pathways (effective \(\zeta\)↑, \(T\_{\mathrm{eff}}\)↑), accelerating unlatching.
- Ground vs excited states. Closed‑shell ground states maximize coherent overlap networks (large effective \(\sum C*{ij}\)) and minimize internal bath coupling; excited and highly mixed states have smaller \(\Delta E_b\) and larger \(T*{\mathrm{eff}}\), hence shorter \(\tau\).
- Falsifiable signatures beyond standard selection‑rule estimates.
  - Isoelectronic series: after accounting for radiative matrix‑element/relativistic trends, residual lifetime shortening with Z should correlate with proxies for internal mode density (bath power \(S\_\xi(0)\)).
  - Rydberg scaling: at fixed environment, increasing n reduces \(C\) (overlap) and increases blackbody dephasing; lifetimes drop faster in warmer/near‑surface conditions consistent with \(\tau\propto e^{\Delta E*b/(k_B T*{\mathrm{eff}})}\).
  - Symmetric bath dependence in composite links (electron–nuclear hyperfine): controlled noise injected on either constituent similarly reduces \(\tau\) (cf. §10.2), distinct from purely local models.

Thus “lower energy states” correspond, mechanically, to deeper, more robustly latched configurations of the internal phase‑link network with smaller effective bath temperature in local units, while “instability” reflects faster Kramers unlatching driven by increased bath power and reduced overlap‑controlled coupling.

Scope and symmetry independence. The stability/entanglement mechanism above relies only on the universal U(1) phase channel and the envelope/SVEA limit; no explicit SU(2) or SU(3) modeling is required for these qualitative and semi‑quantitative predictions. Higher symmetries refine selection rules and line strengths by modifying effective couplings (`C`), barriers (`ΔE_b`), damping (`ζ`), and bath power (`S_ξ(0)`), but they are not prerequisites for latching/unlatching. Hence the lifetime scalings, symmetric bath dependence, and Rydberg trends stated here hold with just the base phase rotor and effective couplings; precise splittings and branching ratios would require the full SU(2)/QED structure.

### 10.6 Joint splittings and Kramers rates under gravitational asymmetry (placeholder)

Motivation. Bundle two parameter‑tight, near‑term tests in one paper: (i) the 1/k excess‑variance scaling for variably blocked multi‑slit experiments (§§2–4), and (ii) entanglement‑link dynamics in asymmetric gravitational environments (this section). The latter targets joint (two‑end) observables that remain after conformal covariance cancels local calibration differences.

Minimal model (recap). The entanglement link is a single nonlocal coordinate θ with effective Lagrangian

\[ L*\theta = \tfrac12 M*{\rm eff} \, \dot\theta^2 - V(\theta), \qquad V(\theta) \approx -K\cos\theta + \lambda\cos^4\theta, \]

with reduced inertia \( M*{\rm eff} = (M*{p,A} M*{p,B})/(M*{p,A}+M\_{p,B}) \), coupling \(K\) from the phase kernel across the bundle, and weak stabilizer \(\lambda\). Barrier and small‑oscillation scales obey (schematic)

\[ \Delta E*b \sim C^2/\lambda, \qquad \omega_0 = \sqrt{V''(\theta*{\rm min})/M\_{\rm eff}}. \]

Tunnel splitting (double‑well) follows WKB:

\[ \Delta \;\approx\; (\hbar*{\rm eff}\, \omega_0/\pi)\, e^{-S/\hbar*{\rm eff}} , \qquad S=\int*{\theta*-}^{\theta*+}\! \sqrt{2 M*{\rm eff} (V(\theta)-E_0)}\, d\theta. \]

Asymmetric environments. Let A and B sit at potentials \(\Phi_A,\Phi_B\) (and possibly acceleration a at one arm). Conformal covariance equalizes local thresholds; joint dynamics depend on the combined environment via small, dimensionless corrections. We parameterize the effective bath and weak redshift of link dynamics as

\[ T*{\rm eff,link} \simeq T_0\,\Big[1 + \alpha_g\, \tfrac{\Phi_A+\Phi_B}{2 c_s^2}\Big] + \tfrac{\hbar}{2\pi c_s k_B}\, a, \qquad \omega*{\rm split,obs} \simeq \Delta/\hbar\_{\rm eff}\, \Big[1 + \mathcal O\big(\tfrac{\Phi_A-\Phi_B}{c_s^2}\big)\Big]. \]

Predictions (to be detailed/derived):

- Symmetric dependence: entanglement lifetime \(\tau*{\rm link} \sim \exp[\Delta E_b/(k_B T*{\rm eff,link})]\) shifts with \(\Phi_A\) and \(\Phi_B\) through the sum \(\Phi_A+\Phi_B\), reversible under arm swapping.
- Tiny frequency shift: \(\delta\omega/\omega \sim \mathcal O(\Delta\Phi/c_s^2)\) for the tunnel splitting; observable only near threshold or with long integration.
- Hysteresis near stability threshold: slow ramps of \(\Phi\) or a drive on a yield latching/unlatching with a loop characteristic of Kramers activation, distinct from smooth Markov dephasing.

Minimal protocols (no absolute standard needed):

- Swap‑asymmetry: place one arm deeper in \(\Phi\); measure \(\tau*{\rm link}\) or \(\Delta/\hbar*{\rm eff}\); swap depths and confirm reversibility of the change.
- Modulate \(\Phi\) or a: apply small, controlled \(\Delta\Phi\) (e.g., altitude change) or acceleration; look for synchronous shifts in \(\ln \tau*{\rm link}\) and \(\omega*{\rm split}\).
- Bath‑slope test: vary local bath temperature at each end; compare \(d\ln \tau\_{\rm link}/d(1/T)\); the higher‑noise end should dominate fragility.

Notes/open derivations (to be added):

- Micro→macro map for \(\alpha_g\) from the phase kernel’s \(\chi(\rho_N)\) and conformal window drift; bound its magnitude vs existing clock and gravitational redshift tests.
- Explicit \(\Delta E*b(\tau_A,\tau_B)\) and \(\Delta(\tau_A,\tau_B)\) in a simple bundle geometry; linearize small \(\Delta\Phi\) to extract \(\partial\ln\tau*{\rm link}/\partial\Phi\).
- Master‑equation treatment (Caldeira–Leggett‑type) for \(\theta\) with two baths; derive \(\Gamma*\varphi\) and \(\Gamma*{\rm esc}\) asymmetry at leading order.
- Control set for technical confounds (timing jitter, phase reference drift, polarization‑dependent loss) that could mimic tiny \(\delta\omega\) or \(\delta\ln \tau\).

Scope/observability caveat. In exact conformal covariance, dimensionless thresholds (e.g., \(\Delta E*b/(k_B T*{\rm eff})\)) are invariant; observable effects arise from small, non‑covariant residuals, asymmetries, or retarded nonequilibrium of the link. Operate near threshold to exploit Kramers exponential sensitivity.

### 10.7 Emergent Configuration Space: The Moduli of Multi-Soliton Fields

Standard quantum mechanics describes $N$ particles in a $3N$-dimensional configuration space. In this framework, the fundamental object is always the real field in 3D (or on the graph). The "configuration space" emerges as the **moduli space** of multi-soliton solutions.

- **Field vs. Wavefunction (The "Lump" Picture):** A single particle is a localized amplitude concentration ("lump") of the field, $\Psi \sim f(x-x_1)$. A multi-particle state is a multi-lump configuration, $\Psi \sim \sum_i f(x-x_i)$. The "wavefunction" $\Psi_{\text{QM}}(\mathbf{x}_1, \dots, \mathbf{x}_N)$ is not a physical wave in a high-dimensional space; it is the probability distribution over the *parameters* (moduli) of the field—specifically, the center positions $\{\mathbf{x}_i\}$—evolved by the stochastic field dynamics.
- **Entanglement in Real Space:** Entanglement corresponds to correlations in the noise driving these moduli. If the phase-link coupling (a physical force in 3D) correlates the stochastic kicks $\varepsilon_1(t)$ and $\varepsilon_2(t)$ experienced by two solitons, their trajectories $(\mathbf{x}_1(t), \mathbf{x}_2(t))$ explore the moduli space in a correlated way. This reproduces the statistics of an entangled configuration-space wavefunction without requiring a literal $3N$-dimensional space. The "potential" $V(x_1, x_2)$ in configuration space is simulated by the non-local phase forces coupling the lumps.

### 10.8 Multi-Particle Entanglement as Network Modes

Generalizing to $N > 2$ particles (e.g., GHZ or W states), entanglement is modeled not as a tensor product structure, but as collective metastable modes of a **phase-link network**.

- **Rotor Network Picture:** The phase-link bundle connecting $N$ droplets forms a fully connected rotor network (resembling a spin glass).
- **Network Modes:**
  - **GHZ States ($|000\rangle + |111\rangle$):** Corresponds to a "ferromagnetic" symmetry-broken sector. The network possesses two deep global basins (all-aligned vs. all-flipped). The system maintains coherence between these topological sectors. Measurement of one particle pins its phase, transmitting a stress wave through the rigid phase network that instantly restricts the available basins for the remaining particles (non-local constraint).
  - **W States ($|100\rangle + |010\rangle + |001\rangle$):** Corresponds to a single excitation ("spin wave") shared across the network, delocalized over all nodes. This mode is robust to the loss of a single node, matching the known robustness of W-state entanglement.
- **Scaling:** This approach scales physically with $N$. It does not require an exponentially large space, but rather a network of $N$ nodes capable of supporting exponentially many metastable latching patterns (basins), which is a generic property of frustrated networks.

## 11. Next steps

- Extend beyond small‑ε to finite‑σ_ε regime; include band‑limited noise spectra and detector integration filters.
- Model partial path correlations with a simple covariance structure for {ε_j} and derive the modified scaling.
- Explore near‑field (Fresnel) regimes and multi‑plane interferometers.
- Connect σ*ε to microscopic parameters (C, M_p, bath S*ξ(0)) to close the loop with entanglement/measurement sections in the main paper.

---

### Appendix A (sketch): Leading‑order variance derivation

Let E = (A/√k) Σ*j e^{i(φ_j + ε_j)}. For small ε, expand e^{iε_j} ≈ 1 + iε_j − ε_j^2/2. Then
I = |E|^2 = (|A|^2/k) [ |S_k|^2 + 2 Im(S_k^\* Σ_j e^{iφ_j} ε_j) + O(ε^2) ].
Taking expectation over ε_j (zero mean), E[I] ≈ (|A|^2/k) |S_k|^2 − (|A|^2/2k) Σ_j |e^{iφ_j}|^2 Var[ε_j] + …
The variance is dominated by the linear term in ε: Var[I] ≈ (|A|^4/k^2) · Var(2 Im(S_k^\* Σ_j e^{iφ_j} ε_j)) = (|A|^4/k^2) · 4 |S_k|^2 Σ_j Var[ε_j] · |e^{iφ_j}|^2 · c(x), with c(x) an O(1) geometrical phase factor. With Var[ε_j]=σ*ε^2 and |e^{iφ*j}|=1, Var[I] ≈ |A|^4 · (σ*ε^2/k) · F\*k(x). Normalizing by E[I]^2 ∝ |A|^4 · (|S_k|^4/k^2) yields the stated ∝ (σ\*ε^2/k) scaling times an O(1) geometric factor.

### Appendix B: Formal Mappings to Standard QM

#### B.1 POVM/Kraus representation of the measurement interaction

Let the initial joint state be $\rho_S \otimes |A_0\rangle\langle A_0|$ and let $U$ implement the system–apparatus interaction. Define Kraus operators by $K_m := \langle A_m| U |A_0\rangle$, where $\{|A_m\rangle\}$ is an (approximately) orthonormal pointer basis. The reduced state is

$\rho_S' = \mathrm{Tr}_A[ U (\rho_S \otimes |A_0\rangle\langle A_0|) U^\dagger ] = \sum_m K_m \rho_S K_m^\dagger$, with completeness $\sum_m K_m^\dagger K_m = I$ because $\sum_m |A_m\rangle\langle A_m|=I_A$.

The associated POVM elements on the system are $M_m := K_m^\dagger K_m$, $M_m\ge 0$, $\sum_m M_m=I$. If the apparatus states correlated with distinct system channels are non‑orthogonal, $\langle A_n|A_m\rangle = \gamma_{mn}$, then in the eigenbasis $\{|m\rangle\}$ the off‑diagonals transform as $\rho_{mn}' = \gamma_{mn} \, \rho_{mn}$ (up to normalisation), realising partial measurement with visibility set by $|\gamma_{mn}|$.

#### B.2 Lindblad master equation (pure dephasing)

For weak, Markovian coupling of a system observable $L$ (which‑way operator, e.g., $L=\sum_m m\,|m\rangle\langle m|$ or projectors onto alternatives) to a broadband phase bath, the reduced dynamics obey

$\dot\rho = -i[H_S, \rho] + \Gamma_\varphi \Big( L \, \rho \, L - \tfrac12\{L^2, \rho\} \Big)$.

In the present framework the dephasing rate is set by the macroscopic pointer current separation and bath power (cf. main text §9):

$\Gamma_\varphi \approx (\Delta p_A)^2 \, S_\xi(0) / (2 \, \hbar_{\mathrm{eff}})$, with $\Delta p_A \approx 2C/M_p$.

This reproduces overlap‑suppressed coherences, $\rho_{mn}(t) = e^{-\Gamma_\varphi t} \, \rho_{mn}(0)$, consistent with $\gamma_{mn} \approx e^{-\Gamma_\varphi T_{\mathrm{int}}}$. For mixed pointers, use fidelity $F(\rho_A^{(1)},\rho_A^{(2)})$ and Fuchs–van de Graaf bounds $1-F \le D \le \sqrt{1-F^2}$; then $V \le F$ provides an operational upper bound consistent with §4.2.

#### B.4 Two‑level (qubit) exemplar

Identify two stable phase‑current orientations with $|0\rangle$ and $|1\rangle$. A weak transverse drive induces Rabi‑like oscillations governed by an effective two‑level Hamiltonian $H = (\Omega/2)\,\sigma_x + (\Delta/2)\,\sigma_z$. Measurement couples $\sigma_z$ to the apparatus phase‑axis; partial alignment implements a weak $\sigma_z$ measurement with Kraus operators $K_\pm \propto (I \pm \eta \, \sigma_z)$, $0<\eta<1$.

#### B.5 Uncertainty relation in the envelope limit

With the Schrödinger envelope (main text §2.4), define $X$ and $P := -i \, \hbar_{\mathrm{eff}} \, \nabla$ acting on $\psi(x,t)$. The commutator $[X, P] = i \, \hbar_{\mathrm{eff}}$ implies $\Delta X \, \Delta P \ge \hbar_{\mathrm{eff}}/2$ within the SVEA regime.
