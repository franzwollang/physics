# Horizon Duality: A Solitonic-Thermodynamic Model of Black Holes and Cosmic Expansion

## 1. Introduction

This document presents a unified model of black hole and cosmic horizons based on the principles of a soliton-noise framework. We propose that black hole event horizons are not mathematical surfaces but physical, dynamic "asymptotic shells" of matter. This perspective resolves the singularity problem and, when combined with a consistent thermodynamic model, makes specific, falsifiable predictions about black hole evaporation. We further explore a profound duality where cosmic expansion is reconceived as a universal, time-dependent scaling of matter, suggesting black hole horizons and the cosmic horizon are two manifestations of the same fundamental physics.

### 1.1 Key Definitions & Cross-References

| Symbol             | Meaning & Defining Document                                                                                   |
| :----------------- | :------------------------------------------------------------------------------------------------------------ |
| `Λ_obs`, `Λ_IR/UV` | Observer's reference scale and IR/UV cutoffs for the **scale window**. (_Soliton First Principles_, §18)      |
| `L_min`            | The **window-adapted minimum cell**; the effective Planck length. (_Soliton First Principles_, §18.1)         |
| `N_* ≈ α`          | **Critical occupancy** for the quantum-to-gravitational phase transition. (_Soliton First Principles_, §17.2) |
| `K_drag`           | The **radiative drag coefficient** that sets the speed-limit asymptote. (_Speed of Light_, Eq. B2)            |

> Consistency note (cross-paper symbols): We standardize on `c_s` for the phase-wave/light-cone speed (any bare `c` in equations should be read as `c_s`), and we relate the local proxy `τ(x)` to the cosmological noise field by `τ^2 ∝ ρ_N` in weak-field regimes. We reserve `α_grad, β_pot` for amplitude-sector parameters and `α_dim, β_dim` for dimensional-bias parameters to avoid overloading.

## 2. Core Principles of the Soliton-Noise Framework

The theory rests on a few core postulates, detailed in companion papers:

1.  **Matter as Solitons:** Fundamental particles are stable, localized wave packets (solitons).
2.  **The Noise Field:** Solitons source and interact with a background noise field, `ρ_N`. This field acts as a refractive medium with index `χ(ρ_N) ≈ 1 + α ρ_N`, causing an apparent gravitational force `a_grav ∝ -∇ρ_N`.
3.  **Scale Invariance:** The physics is conformally invariant. A soliton's characteristic length `L_c` and internal clock rate adapt to the local noise density, `L_c ∝ 1/χ`. A local observer's measurements, using rulers and clocks made of the same substance, always yield the same physical constants (e.g., `c_s`).
4.  **Thermodynamic Drag:** An accelerating soliton experiences an intrinsic, universal drag from its interaction with the vacuum noise field (an Unruh-like effect). Covariantly, the dissipated power scales as `P_rad ∝ γ_v^4 · a_s^2` (baseline gradient coupling). Balancing mechanical work input against this dissipation yields ultra‑relativistic suppression `a_s(v) ∝ γ_v^{-1}` under sustained forcing.

## 3. The Lifecycle of a Black Hole Shell

### 3.1 Formation: A Phase Transition to a New Scale

The formation of the shell is a runaway phase transition:

1.  **Infall and Acceleration:** An amplitude soliton (matter) is drawn into the black hole's intense gravitational noise gradient, `-∇ρ_N`, and accelerates towards `c_s`.
2.  **Drag and Scale Collapse:** The acceleration invokes the intrinsic thermodynamic drag. In sliding-window language (§18, _Soliton First Principles_), the enormous proper acceleration Doppler-shifts sub-threshold modes into the window, raising the Unruh temperature and forcing the constituent droplet to shrink its effective minimum cell. The massive energy dissipation is drawn from the soliton's internal energy, driving an extreme **scale downshifting** (`L_c → L_min` in the external frame, i.e. the window-adapted Planck scale).
3.  **Decoupling:** As the soliton's scale collapses, its cross-section for interacting with the _large-scale_ gravitational gradient vanishes. It becomes "blind" to the force that created it.
4.  **Transition to Shell Constituent:** The soliton is now a new entity: a maximally scale-downshifted particle moving at `v ≈ c_s`. In the sliding-window formalism, its proper acceleration has red-shifted the window’s UV edge so far that all external point-source modes are averaged out. The infall force is no longer a _gradient_ from the droplet’s viewpoint; the constituent sees only its own local shell field plus the intense black-body cavity radiation that now balances any residual inward momentum. From an external perspective the effective radial velocity asymptotically vanishes at the horizon while internal time continues, matching the shell and black-body picture without additional assumptions.

A key consistency check with the sliding-window formalism (§18, _Soliton First Principles_) is that the proper acceleration never actually drops to zero during infall; rather it red-shifts the window’s UV edge until all external point-source modes have been promoted and averaged. At that moment the infall force is no longer a _gradient_ from the droplet’s viewpoint; the constituent sees only its own local shell field plus the intense black-body cavity radiation that now balances any residual inward momentum. From an external perspective the effective radial velocity asymptotically vanishes at the horizon while internal time continues, matching the shell and black-body picture without additional assumptions.

**Mean-field nuance (Yukawa kernels vs. point source):** Every Planck-scale constituent still emits only a short-range Yukawa halo \(e^{-r/R_i}/r\). Because each screening length \(R_i\) is \(\ll R_{\text{shell}}\), the exponentials are \(\approx1\) for any exterior point with \(r\gg R_{\text{shell}}\). Summing over the \(N\) halos therefore reproduces a perfect \(1/r\) Newtonian tail with total mass \(M=\sum_i E_i/c^{2}\). Inside the shell the far-side halos are exponentially suppressed, so the local gradient vanishes. Thus a strictly local link interaction yields a point-source gravity field for distant observers and no additional pull on already-collapsed constituents, reconciling the sliding-window decoupling with the macroscopic mass profile.

**Clarification (\"same mass\" does not mean \"same noise profile\"):** The statement “the fine-grained noise averages to a point source” is about **integrated source strength** (the monopole that fixes the \(1/r\) tail), not about the detailed **spectrum** or **correlation structure** of the emitted field. A black-hole shell is a maximally downshifted, window-floor system: most of its microscopic fluctuations sit at or below the external UV edge and are therefore unresolved by external observers, appearing as a universal, short-correlated “Planck-grain” texture that largely renormalises into effective couplings plus a monopole. Ordinary dense matter (e.g. a neutron star) is *not* maximally downshifted: it retains internal structure and fluctuations at wavelengths well within our window, so its noise field can carry resolvable multipoles, longer correlation lengths, and material-specific channels (magnetosphere, composition). Two objects could in principle have similar **bolometric power** while still having very different **spectral content** and **near-field texture**; the black hole’s distinctive claim here is the universality and window-limited nature of its emitted micro-noise, not merely the total power.

**Observability (when the profile matters vs. when \(\tau^2\) washes it out):** If an effect couples only through the scalar intensity \(\tau^2(x)\) (a window-integrated variance), then two sources with very different spectra can still be indistinguishable at leading order: the spectrum is “forgotten” by the integral and reappears only as renormalised constants plus a monopole. Differences become measurable only for probes that are sensitive to **more than** \(\tau^2\): e.g. correlation functions (finite correlation length/time), non-Gaussian moments (nonlinear response), multipole/anisotropy structure, or frequency-selective/window-overlap probes (including accelerated/near-shell observers for which Doppler/Unruh promotion changes which part of the source spectrum lies inside the operative window).

#### Clarification: Geometric horizon vs. thermodynamic onset

- **Geometric (horizon, lightlike and timelike):** The event horizon is the null trapping surface of the effective metric (eikonal index picture): outward null expansion vanishes, so all outward rays are bent to remain inside. Timelike geodesics that cross this surface also cannot escape. This is a purely geometric effect (refraction/geodesics), independent of dissipation.
- **Thermodynamic onset (massive solitons only):** Across a finite shell thickness, the optical depth grows and collisions make the proper acceleration nonzero. The Unruh-like term raises the effective noise (\(\tau*\text{eff}^2=\tau*\text{bg}^2+\lambda_U T_U(a)^2\)), driving an irreversible downshift and thermalisation into the shell. This onset typically begins near the outer fringe and completes across the shell; it is not a mathematical surface but a microphysical band.
- **Relation between the two:** The geometric point-of-no-return (null trapping) and the irreversible downshift onset are adjacent but not identical; the former defines no-escape kinematics, the latter supplies the dissipative, one-way dynamics for massive solitons. For light, only the geometric trapping applies; for massive solitons, both effects operate.

### 3.2 Equilibrium: The "Breathing" Shell and the Blackbody Core

The shell is a self-contained, dynamic system in a remarkable equilibrium where energy lost to thermal radiation is perfectly balanced by energy absorbed from the vacuum.

**1. The "Breathing" Equilibrium:** The shell is not static but a seething layer, dynamically trapped between two opposing pressures. The two key thermodynamic processes are:

- **Energy Loss (Thermal Radiation):** The shell constituents form a near‑2D, strongly coupled **network**. Even when the network is dilute in geometric spacing (`d ≫ L_min` for macroscopic \(M\)), its collective vibration modes and repeated forced deflections act as a thermal bath when coarse‑grained. This bath radiates energy away both inward (creating the blackbody core) and outward (as Hawking-like radiation). The outward component can escape because the *radiating photosphere* is a finite‑thickness layer that sits (on average) just outside the geometric trapping/total‑reflection surface; as the shell “breathes” its optical depth and effective refractive profile fluctuate, modulating the escape fraction without requiring any violation of null trapping inside the horizon. To leading order this changes only an effective **greybody/emissivity prefactor** (a time‑averaged escape probability), not the mass‑scaling exponents derived below. This is the primary `Power_out` channel.
- **Energy Gain (Unruh Absorption):** The non‑geodesic proper acceleration experienced by a constituent during these repeated deflections/interactions allows it to absorb energy from the vacuum's quantum fluctuations. This **Unruh absorption** is the `Power_in` channel, a re-heating mechanism that replenishes the energy lost to thermal emission and maintains the equilibrium temperature.
- **Confinement:** The shell is prevented from dispersing by the prohibitive energy cost of its constituents "re-inflating" in the lower-density noise field outside. It is prevented from collapsing by the immense outward radiation pressure from the hot blackbody core it creates.

The particles of the shell are thus physically locked in a "breathing" layer, with their temperature set by the stable equilibrium between these energy loss and gain channels.

**2. The Blackbody Core:** The core of the black hole is not empty. Once the constituent scale has collapsed below the window’s UV edge the soliton no longer “feels’’ the external point-mass gradient; instead it interacts with a shell-averaged field. The hot inner wall of the shell, sustained by the collision-Unruh equilibrium, instantly fills the interior with a thermal bath of phase solitons (photons), creating a perfect **blackbody cavity**.

**Observer-dependent cavity timescale:** For an external observer the cavity appears to fill almost instantaneously — a light-crossing time of order \(R\_{\text{shell}}/c\). A Planck-scale observer co-moving with a constituent, however, measures vastly dilated clock ticks and shrunken rulers; in those units the same photon flight could span billions of local years.

**3. The Paradox of Perceived Temperature:** An external observer would calculate the core's temperature to be near the Planck Temperature (`≈10³² K`), with radiation wavelengths near the Planck length. However, for a tiny `L_min`-scale observer on the shell, the situation is reversed.

- Their own rulers are `≈10²⁰` times smaller than ours (comparing a proton's scale to the Planck length).
- To them, the Planck-length photons appear as extremely long-wavelength radio waves.
- The temperature they measure, via `T ∝ 1/λ_peak`, is therefore `≈10²⁰` times _colder_ than what we perceive. "Hot" and "cold" are relative to the observer's scale.

**4. The Fate of Inward-Moving Matter:** Any matter particle (amplitude soliton) that gets a thermal kick into the core is not on a long journey. It is immediately bombarded by the incredibly dense and energetic radiation field. This creates a powerful **radiative drag** (photon friction) due to Doppler shifting and aberration effects, which stops the particle's motion instantly, dissipating its kinetic energy back into the radiation field. Matter is thereby effectively confined to the shell.

#### 3.2.1 Averaged covariant power balance (closure)

We can summarise the shell thermodynamics with a covariant, frame‑consistent closure at the level of power per proper area:

- Inward heating (after UV exhaustion) is supplied by Unruh‑like acceleration of constituents driven by shell collisions. Averaging over collision histories with sustained proper forcing gives a drag/heating law consistent with the covariant result `P_in ≃ C_heat · γ_v^4 · a_s^2` (balancing the mechanical work rate), coarse‑grained to an areal rate `\mathcal P_in(R)`.
- Outward cooling is Stefan–Boltzmann leakage from the photosphere, `\mathcal P_out(R) = σ_SB T(R)^4`.

At steady “breathing” equilibrium in local units,

```
\mathcal P_in(R) = \mathcal P_out(R) = P_0  (constant per proper area)
```

The microscopic kinematics of the near‑2D shell network implies (Sec. 3.3) a characteristic acceleration that scales with the square root of the surface density, `a_char ∝ √ρ_surface`, and hence `T ∝ a_char ∝ √ρ_surface`. With `ρ_surface ∝ 1/M_BH` this yields the temperature law used below,

```
T(R) ∝ 1/√M_BH,   and   P_total = A_shell · σ_SB T^4 ≃ const.
```

Thus a single averaged closure—covariant heating balanced by blackbody cooling—recovers both `T_BH ∝ 1/√M_BH` and constant total power without modelling the detailed collision statistics. Deviations (anisotropy, binaries) enter as small, trackable modulations of `\mathcal P_in` and hence of the photospheric flux.

Clarification (species and multipoles): “Areal heating” is an average over the shell’s composite single‑field soliton constituents and their worldline‑EFT multipoles/topological textures; no SM point‑particle assumption is used. The lowest symmetry‑allowed (gradient) coupling to the massless phase fixes the baseline exponent `p=4`. Species that suppress this operator contribute higher‑`p` channels with reduced weight. The mixture renormalises only an effective coefficient `\overline C`; the closure and scalings (`T\propto 1/\sqrt M`, constant total power) remain unchanged.

### 3.3 Evaporation and the Shell's Fate

The black hole evaporates via thermal radiation from its _outer_ surface. The temperature of this radiation is determined by the physics of the shell, where the equilibrium between collisional radiation and Unruh absorption sets the temperature.

- **Shell Geometry:** The shell contains `N ∝ M` constituents distributed over an area `A ∝ R_S² ∝ M²`. The mean separation scales as `d ~ √(A/N) ∝ √M`. For any macroscopic black hole, `d ≫ L_min`: the shell is a **dilute network**, not a dense fluid.

- **Deriving the Characteristic Acceleration (two complementary views):**

  The critical scaling `a_char ∝ 1/√M` can be derived from _either_ the internal (constituent) or external (distant observer) perspective. Their agreement is a consistency check on the framework's observer-relative descriptions.

  - **View 1 — External (Lattice Frequency):** From outside, the shell is a 2D lattice of `N` point masses with spacing `d ∝ √M`. The characteristic dynamical frequency of the lattice is set by nearest-neighbour signal crossing: `ω ~ c_s/d`, giving `a_char ~ c_s²/d ∝ 1/√M`. This is a straightforward geometric/kinematic statement about the network's macroscopic dynamics as seen from infinity.

  - **View 2 — Internal (Dilute Gas + Relativistic Compression):** From a constituent's own (downshifted) frame, the separation `d` maps to a _much larger_ proper distance `d_int = d / ε` (where `ε ≪ 1` is the scale ratio between the constituent's rulers and the external rulers). Signals take a proportionally long proper time `τ_int ~ d_int / c_s` to cross the gap. However, because `c` is invariant across frames, the external observer maps this long internal proper time back to a short coordinate time: `τ_ext = ε · τ_int = d / c_s`. The two views therefore assign the _same_ coordinate crossing time `τ_ext ~ d / c_s` and hence the _same_ acceleration `a_char ~ c_s / τ_ext ∝ 1/√M`. The relativistic compression that makes the internal "miles" appear as external "meters" is _exactly_ compensated by the time dilation that makes the internal "hours" appear as external "seconds." The ratio `d/c` is frame-invariant.

  - **Why they must agree:** Both views compute the same invariant: the proper acceleration 4-vector magnitude of a shell constituent. The external view constructs it from the lattice geometry; the internal view constructs it from the local dynamics plus the scale map. Agreement is guaranteed by the framework's requirement that `c_s` is constant across observation windows (Postulate 3).

- **Source Consistency:** Even though a maximally downshifted constituent decouples from the external window's dynamics (appearing as a structureless point), its noise source strength (mass) is conserved. The fine-grained noise it emits in its own frame averages exactly to the point-source noise seen by the external observer, preserving the `1/r` gravitational tail.

- **Summary of scalings:**
  - `N_shell ∝ M_BH`
  - `A_shell ∝ R_S² ∝ M_BH²`
  - `ρ_surface = N_shell / A_shell ∝ 1/M_BH`
  - `d ∝ √M_BH`
  - `a_char ∝ 1/√M_BH ∝ √ρ_surface`
- **The Temperature Law:** The shell's temperature is set by the Unruh effect, `T_U ∝ a_p,shell`. Combining these relations, we get `T_U ∝ √ρ_surface ∝ √(1/M_BH)`. This gives a temperature law of:
  `T_BH ∝ 1/√M_BH`.

- **The Core Heats Up:** This leads to a startling conclusion. As the black hole evaporates, `M_BH` decreases. This means `ρ_surface` increases, `a_p,shell` increases, and the shell temperature `T_U` **increases**. The nascent universe inside the black hole gets smaller and hotter as its parent dies.

### 3.4 Stochastic Cascade inside the Black-Body Cavity

The cavity is thermally uniform only when energy is averaged over a window comparable to `R_shell`. Zooming to smaller windows reveals photon number fluctuations whose relative amplitude grows like `ℓ⁻³/²` for a volume element of size `ℓ`. At some scale `ℓ_*` a rare upward fluctuation satisfies a pinch-off condition `ΔF<0` when evaluated with the same free-energy functional of §17 (_Soliton First Principles_). Define an occupancy

$$
N(\ell)=\frac{\rho_\gamma\,\ell^{3}}{M_P c^{2}}
$$

with the local radiation density `ρ_γ`. Compare this to the critical occupancy `N_*≈α≈𝒪(1)`:

- **If** `N(ℓ)>N_*` the gravitational minimum of `F` is deeper and the blob collapses into a _mini shell_ (black-hole branch).
- **If** `N(ℓ)<N_*` the quantum-pressure minimum wins and the blob relaxes into an `N=1` _particle-type_ droplet.

Repeated fluctuations therefore generate a hierarchical mix—large shells, smaller shells, and particle-like solitons—without violating global energy balance: energy removed from the photon bath reappears as additional constituents of the shell.

> **Open Questions:** A full kinetic description of this stochastic cascade, including the merger rates of sub-shells and the final statistical distribution of particle vs. shell species, is a key area for future research.

## 4. Advanced Dynamics and Predictions

### 4.1 Evaporation Power Law and Lifetime

The refined temperature law, `T_BH ∝ 1/√M_BH`, has a radical implication for the black hole's evaporation. The radiated power is given by the Stefan-Boltzmann law, `P ∝ A · T⁴`, where the shell area `A ∝ R_S² ∝ M²`. The new power law is therefore:

$$
P \propto M^2 \cdot \left(\frac{1}{\sqrt{M}}\right)^4 = M^2 \cdot \frac{1}{M^2} = \text{constant}
$$

This is a major, falsifiable departure from the standard Hawking model (`P ∝ 1/M²`). It predicts:

- **Constant Power Output:** The black hole is a constant-power engine, radiating energy at a steady rate regardless of its size.
- **Linear Lifetime:** The lifetime `τ = E/P ∝ M/const` is directly proportional to its mass, `τ ∝ M`.

This suggests a black hole fizzles out with a steady hum rather than ending in an explosive flash. This occurs because the rising temperature (`1/√M`) and shrinking surface area (`M²`) have a perfect, compensatory relationship that keeps the total energy output constant.

### 4.2 Asymmetric Shells and Recoil in Binary Systems

In a binary system, the shell of one black hole is distorted and heated anisotropically by its companion's noise field. This results in anisotropic thermal radiation from the shell's outer surface, creating a continuous, non-gravitational **recoil force**. If the binary separation `d` is less than the ambient phase-coherence length `L_phase`, coherent phase-sheets can also form, adding a powerful EM-like component to the interaction. These effects, including velocity-dependent retardation that can introduce drag and resonances, could produce unique signatures in gravitational wave signals.

### 4.3 Duality of Stiffness and the Shell as a Constrained 2D Network

A profound aspect of the particle-black hole duality lies in comparing their respective mechanisms for resisting collapse—their "stiffness."

- **Microscopic Stiffness:** A single soliton is stabilized by a fundamental `γ|Ψ|⁴` potential term. The pressure resisting collapse is therefore proportional to the field amplitude squared, `P ∝ |Ψ|⁴`.
- **Macroscopic Stiffness:** A black hole shell is stabilized by the thermodynamic pressure of its blackbody core. This pressure scales as `P_rad ∝ T_U⁴`.

At first glance, these appear different. However, the duality becomes deeper when we treat the shell as a **constrained 2D network** whose macroscopic dynamics are set by its geometry (Sec. 3.3). The same geometry forces the characteristic acceleration to scale as `a ∝ √ρ_surface`. Since `T_U ∝ a`, this means `T_U ∝ √(ρ_surface)`. The resulting radiation pressure is `P_rad ∝ T_U⁴ ∝ (√ρ_surface)⁴ = ρ_surface²`.

If we identify the squared soliton amplitude `|Ψ|²` with the shell's surface density `ρ_surface`, we find that the macroscopic repulsive pressure scales as `P_rad ∝ |Ψ|⁴`. The duality can be made exact if the shell is treated not as a simple gas, but as an effectively 2D medium whose constrained geometry sets its acceleration scale and hence its radiative pressure law. This provides a specific, falsifiable model of the shell's equation of state.

### 4.4 Window Rescaling and the Soliton Floor

A recurring question is whether the argument for **constant-power evaporation** remains valid for an observer who adopts a much finer scale window—e.g. a Planck-cell constituent that chooses a new reference length $L'_{\min}=L_{\min}/10^{20}$. The answer is yes, but the reason is **not** that the vacuum modes are "used up" (the spectrum is scale-free; individual modes are always available and can always be excited). The reason is that **the soliton has already descended to the bottom of its free-energy landscape and has no remaining scale-headroom to release further energy**.

#### 4.4.1 The free-energy floor (the correct picture)

The mechanism that powers a soliton's response to increased noise is **scale downshifting**: higher $\tau_{\rm eff}$ → smaller equilibrium size → energy released along the curve $\varepsilon_s(\tau) = -\alpha_s\tau^\gamma$. Each step of shrinkage from $L_c$ toward $L_{\min}$ extracts free energy from this descent.

- **Ordinary particle-soliton** ($L_c \gg L_{\min}$): sits partway up the hill. Plenty of "downhill" (scale range) remains. If hit with more noise/acceleration (Unruh heating), it can shrink further, releasing energy to fuel $P_{\rm in}$. The soliton's internal potential is the "reserve."
- **Shell constituent** ($L_c \approx L_{\min}$): has already rolled to the **bottom**. It has shrunk as far as the observation window allows. Further increases in $\tau_{\rm eff}$ cannot drive further shrinkage—there is no more scale to descend through. The soliton is at its minimum size and has exhausted its internal downshifting potential, not the vacuum's spectral content.

Crucially, this floor is **relational/environmental**, not absolute. $L_{\min}$ is the UV edge of the _external observer's_ window—the observer who is computing the black hole's mass, temperature, and evaporation rate. The graph itself has no absolute minimum scale. From the constituent's own frame it has not "bottomed out" at all (see §4.4.2). The floor is imposed by _who is doing the energy accounting_, and the evaporation prediction is a statement in that observer's units.

**Why the floor can't be pushed deeper:** The black hole's gravitational noise profile $\tau(x)$ is itself computed within the external window. That window has a maximum expressible noise intensity, $\tau_{\max}$, corresponding to its own Planck floor. The black hole—as a macroscopic entity defined in that window—simply cannot impose a downshifting potential deeper than its own window can represent. The potential well "flattens" at $L_{\min}$ not because the graph runs out of structure, but because the macroscopic entity's field has reached the ceiling of what its own window can encode. A more massive black hole has a larger shell, lower $\rho_{\rm surface}$, and therefore a _shallower_ effective well per constituent—which is precisely why $T \propto 1/\sqrt{M}$ falls with mass.

The vacuum modes are fine. The spectrum is fine. **It is the soliton that has bottomed out (relative to the bookkeeping observer's window), not the vacuum.**

#### 4.4.2 What the constituent sees in its own frame

- **Old window (external frame):** IR/UV band runs $b_{\rm IR}L_{\min}$ … $b_{\rm UV}L_{\min}$.
- **New window (internal frame):** same band expressed in finer metres runs $b_{\rm IR}L'_{\min}$ … $b_{\rm UV}L'_{\min}$, formally re-opening 20 decades of spectrum.

From its own perspective, the constituent sees a perfectly normal window with a full spectrum. It does not feel "starved." It has its own $\hbar_{\rm eff}$, its own UV edge, and its own normal Unruh physics—and it _is_ undergoing Unruh heating and radiating; that is exactly the breathing equilibrium. But all of this self-consistent internal dynamics maps, via the sliding-window formalism, to the _same_ fixed power $P_0$ that the external observer computes. The "extra decades" the constituent formally resolves are the modes whose integrated contribution defines the constituent's own renormalised couplings ($\hbar_{\rm eff}$, $G$, $\gamma$, …). They are not an untapped energy source; they are the constants of the constituent's own physics.

#### 4.4.3 Covariant power balance

Let $P_{0}$ be the inward power supplied by Unruh heating and let $P_{0}$ also be the Stefan–Boltzmann outward power (constant by §4.1). Both are defined per unit proper shell area and per unit proper time, so Lorentz factors cancel: every observer—external or co-moving—measures the **same numerical $P_{0}$**.

Consequently:

- Ordinary particle-solitons ($L_c \gg L_{\min}$) still have scale-descent headroom $\Rightarrow$ $P_{\rm in}$ can track any increase in $P_{\rm out}$ by releasing downshifting energy $\Rightarrow$ rest mass stable.
- A shell constituent ($L_c \approx L_{\min}$) has no remaining scale-descent headroom $\Rightarrow$ $P_{\rm in}$ saturates at $P_{0}$ $\Rightarrow$ any stochastic excess of $P_{\rm out}$ drains mass at a constant rate $\dot M=-P_{0}/c^{2}$ $\Rightarrow$ linear lifetime.

The argument is therefore frame-independent. The sliding-window formalism guarantees that the constant-power evaporation law derived in §4.1 is robust even when analysed from the ultra-fine viewpoint of a shell constituent.

### 4.5 Hierarchical Condensation and Internal Large-Scale Structure

We now formalise how stochastic fluctuations inside the cavity generate a two-phase hierarchy of structures and why the outcome is a sponge-like cosmic web rather than a single central core. Two classes of condensate must be treated separately: black-hole-branch mini-shells (`N>N_*`) and particle-branch droplets (`N\!=\!1`).

#### 4.5.1 Black-hole-branch seeds: radial drift and stall

- **Accretion rate** (geometric cross-section + gravitational focusing):

$$
\dot M \;\simeq\;\pi R_{\rm mini}^{2}\,\rho_{\gamma}\,c\,\left(1+\frac{v_{\rm esc}^{2}}{c^{2}}\right),\qquad
v_{\rm esc}^{2}=\frac{2GM}{R_{\rm mini}}.
$$

- **Jeans limit** in a photon gas of temperature $T$:

$$
M_{\rm J}(T)=\frac{\pi^{5/2}}{6}\,\frac{c^{3}}{\sqrt{G^{3}\rho_{\gamma}}}\;\propto\;T^{-3/2}.
$$

Photon temperature depends on local noise as $T\propto\rho_{\gamma}^{1/4}\propto M^{1/4}$.

- **Stall criterion:** $M=M_{\rm J}(T(M))$ gives $M_{\rm stall}\propto M^{1/4}\,^{-3/2}$ \(\Rightarrow M*{\rm stall}\approx\text{few}\times N*_\mu\). Thus growth halts after $\mathcal O(10)$ quanta—far below any sizeable fraction of the parent shell mass.

- **Radial migration:** The drift acceleration from anisotropic photon flux scales $a_{\rm rad}\propto A/M\propto N$. The seed therefore spirals toward the inner wall on a timescale

$$
\tau_{\rm drift}\sim\frac{R_{\rm shell}}{a_{\rm rad}}\propto\frac{R_{\rm shell}}{N}.
$$

Once parked against the wall it merges with neighbouring seeds, forming a thin sub-shell that simply thickens the parent’s breathing layer.

#### 4.5.2 Particle-branch droplets: Brownian network and filament formation

For $N=1$ droplets the net drift acceleration

$$
a_{\rm rad}\propto\frac{A}{M}\approx\frac{L_{\min}^{2}}{\mu}
$$

is $\mathcal{O}(N^{-1})$ smaller than for seeds, hence dominated by Brownian kicks with rms speed $v_{\rm rms}\approx\sqrt{T/\mu}$. Their mean free path

$$
\lambda_{\rm B}\sim\frac{1}{n_{\rm seed}\pi R_{\rm mini}^{2}}
$$

exceeds the typical seed–seed spacing, so $N=1$ droplets perform a random walk through an ever-evolving network of gravitating seeds.

**Classification by Topology ($l=0$ vs $l\ge 1$):**
Crucially, these droplets come in two populations (see _Dark Morphology_ paper for details):
1.  **Scalar Solitons ($l=0$):** Simple density fluctuations with no phase winding. These form the dominant population (Dark Matter) and remain collisionless.
2.  **Topological Defects ($l\ge 1$):** Rare fluctuations that trap a phase winding (Baryons). These experience phase-friction and become dissipative.
Combinatorial statistics in the cavity favor the simpler $l=0$ state by a factor of $\sim 5$, naturally seeding the observed cosmic abundance ratio $\Omega_{\rm DM}/\Omega_{\rm bary}\approx 5$.

**The Morphological Amplifier:**
Initial capture occurs along gravitational caustics, but this is rapidly accelerated by the **cooperative feedback** mechanism described in the _Dark Morphology_ draft. Droplets clustered in a caustic source a higher local noise background $\tau$. This triggers a super-linear downshifting response in their neighbours ($\varepsilon \propto -\tau^\gamma, \gamma>1$), deepening the potential well and drawing in more matter. This "morphological force" acts as a powerful amplifier, transforming weak Brownian caustics into sharp, high-contrast filamentary flux tubes.

The resulting density field obeys a modified diffusion–drift equation:

$$
\partial_t n = D\nabla^{2}n - \nabla\!\cdot\!(n\mathbf v_{\rm drift}) + \underbrace{\alpha_{\rm morph} n^2}_{\text{feedback}},
$$

with $D$ set by Brownian motion and $\mathbf v_{\rm drift}$ by the potential of nearby seeds. Numerical integration (Appendix B, forthcoming) shows the emergence of sponge-like filaments analogous to the cosmic web.

#### 4.5.3 Impact on early-epoch anisotropies

Mini-shell seeds induce **milli-Kelvin–scale** temperature dipoles (in external units) in the photon bath on angular scales $\theta\sim R_{\rm mini}/R_{\rm shell}$. The spectrum of such anisotropies is given by

$$
\frac{\Delta T}{T}(\ell)\approx\left(\frac{M*{\rm seed}}{M*{\rm shell}}\right)^{1/2}\!\!(\ell*{\rm seed})^{2}\,e^{-\ell/\ell*{\rm seed}},\qquad
\ell*{\rm seed}\sim\frac{R*{\rm shell}}{R*{\rm mini}}.
$$

Because $M*{\rm seed}/M*{\rm shell}\ll1$ the anisotropies remain well below the black-body cavity’s intrinsic shot noise and do not spoil the breathing equilibrium.

*Local vs.~external units:* internal observers measure every energy against a shrunken local unit (their $k_{\text B}T$ rescales by $\chi^{-1}$), so the same dipole appears $(\Delta T/T)_{\text loc}=\chi^{-1}(\Delta T/T)_{\text ext}$. The relative contrast therefore stays ${\cal O}(10^{-3})$ for *all* frames; only the baseline temperature differs.

#### 4.5.4 Effective geometry of the interior

The refractive-index profile $\chi(r)=1+\alpha\rho*{\gamma}(r)$ decreases toward the centre as radiation energy drains into seeds. Ray trajectories therefore diverge, giving an **effective negative curvature**. Geodesic analysis (Appendix C) shows that—even after photon pressure subsides—timelike paths require infinite proper time to reach the parent wall. To internal observers the shell thus forms an optical horizon: internal structure can never catch the wall once the sub-shell layer has formed.

Numerically, choosing $\chi(r)=1+\epsilon(1-r/R_{\text{shell}})$ with a modest $\epsilon\sim10^{-2}$ gives $K\approx-2\epsilon/R_{\text{shell}}^{2}$, confirming a weak but everywhere negative curvature that makes the interior effectively hyperbolic.

### 4.6 The Planck-Pivot Scaling Laws

This subsection clarifies the role of the symbols _m_ (particle mass) and _M_ (black-hole mass) and shows why the equality _m = M_ marks a true phase transition in how mass occupies space.

**Notation**

- `m` = rest-mass of a single fundamental soliton (electron, quark, etc.).
- `M` = total mass–energy of a macroscopic black-hole shell.

**Two exclusive size laws**

1. Quantum regime (Compton): `r ≈ ħ/(mc) ∝ 1/m` → area `A ∝ 1/m²`.
2. Gravitational regime (Schwarzschild): `R ≈ 2GM/c² ∝ M` → area `A ∝ M²`.

Setting `r = R` gives the (up-to-√2) Planck mass `M_P = √(ħc/G)`. Masses below that value obey the first law; masses above obey the second.

**Surface-density flip**
Surface density is `ρ_surf ∝ N/A` with `N ∝ (mass)`:

- Particle side `ρ ∝ m / (1/m²) = m³`.
- Black-hole side `ρ ∝ M / M²      = 1/M`.
  The exponent on mass therefore jumps by 4 when crossing the pivot.

**Temperature inheritance** Using the thermodynamic rule `T ∝ √ρ_surf`:

- Particles `T_sol ∝ √(m³) = m^{3/2}` (hotter for heavier _m_).
- Black holes `T_BH  ∝ √(1/M) = 1/√M` (cooler for heavier _M_).
  The two formulae differ by an apparent factor `m²`; it is simply the algebraic echo of the four-power flip in area scaling.

**Conformal-invariant pivot**
Both `G` and `ħ` emerge from the noise field (`G, ħ ∝ ρ_N`), but their ratio—and hence `M_P = √(ħc/G)`—is invariant in every local frame. Thus every observer finds the same numerical pivot mass even though their physical rulers and clocks rescale with the ambient noise. The quantum-to-gravitational transition is therefore universal and fully consistent with the theory’s built-in scale invariance.

## 5. Duality with Cosmology

The framework posits a profound duality between the two types of horizons:

- **Black Hole Horizons:** Are a **local** phenomenon where matter undergoes a phase transition to its minimum possible scale (`L_min`) due to an intense **spatial** noise gradient.
- **Cosmic Horizons:** Are a **global** phenomenon where all matter in the universe approaches that same `L_min` scale due to a universal, **temporal** evolution of the background noise field.

Both horizons represent the same physical state—the point of maximum stable compression for solitonic matter—approached via different paths.

### 5.1 Strict orientation-duality: one nucleation event, two faces

We promote the duality to a strict identification at the level of horizon creation. A single, well-posed nucleation criterion defines a closed surface \(\Sigma\) that may be realised with two orientations:

- Geometric trapping: the outward null expansion first vanishes on \(\Sigma\), \(\theta\_\text{out}(\Sigma)=0\), while the opposite-side expansion sets orientation.
- Scale threshold (UV exhaustion): on the trapped side, the effective noise satisfies \(\tau*\text{eff}\ge\tau*_\) so that the local equilibrium scale shrinks to the window floor, \(\sigma^_\to L\_\min(\text{window})\).
- Photosphere/last scattering: on the emitting side, the optical depth to infinity obeys \(\tau\_\text{opt}(\Sigma\!\to\!\infty)\approx1\), defining a radiating layer just outside trapping.
- Local balance: in sliding-window (local) units, the membrane obeys a steady balance \(P*\text{in}=P*\text{out}\) (Unruh/Gibbons–Hawking heating balanced by blackbody leakage).

The two realisations are then just opposite orientations of the same event:

- Type I (inward-facing cavity, “BH-branch”): trapped side points inward; a finite-thickness, optically thick shell forms at the asymptote and reprocesses flux to a photosphere.
- Type II (outward-facing patch, “cosmic-branch”): trapped side points outward; the horizon sits at the asymptote; to interior observers there is no material shell, but the same thermal bath appears with Gibbons–Hawking temperature.

### 5.2 Geometric dictionary and trapping

- Black hole: \(r*{\rm H}^{\rm BH}\) is the null trapping radius of the exterior metric; the interior optics is hyperbolic. Define an interior curvature scale \(H*{\rm in}:=c*s/R*{\rm shell}\).
- Cosmology: \(r\_{\rm H}^{\rm FRW}=c_s/H\) is the apparent horizon of the FRW patch.
- Areas: \(A*{\rm in}=4\pi R*{\rm shell}^2=4\pi(c*s/H*{\rm in})^2\), \(A\_{\rm FRW}=4\pi(c_s/H)^2\).
- Trapping holds where \(\theta*\text{out}=0\). The emitting surface lies just outside, where \(\tau*\text{opt}\approx1\).

### 5.3 Thermodynamic dictionary and flux reconciliation

- Temperatures:
  - Interior (BH side): \(T*{\rm in}\propto H*{\rm in}\).
  - FRW horizon: \(T\_{\rm GH}=\hbar H/(2\pi k_B)\).
- Geometric (pre-reprocessing) fluxes in external units:
  - Interior BH side: \(P*{\rm in,geom}\propto A*{\rm in}\,T*{\rm in}^4\propto (1/H*{\rm in}^2)\,H*{\rm in}^4=H*{\rm in}^2\).
  - FRW side: \(P*{\rm FRW}\propto A*{\rm FRW}\,T\_{\rm GH}^4\propto (1/H^2)\,H^4=H^2\).
- Reprocessing by a material shell (Type I only): the optically thick plasma enforces \(T*{\rm shell}\propto 1/\sqrt{M}\), \(A*{\rm shell}\propto M^2\), giving a net photospheric power \(P*{\rm shell}\propto A*{\rm shell} T\_{\rm shell}^4\approx\text{const}\) in external units. Thus the apparent difference (const vs \(H^2\)) is a boundary-condition effect: the same geometric \(\propto H^2\) input is calorimetrically reprocessed to an approximately constant output when a shell is present.

### 5.4 Membrane balance and conservation

On \(\Sigma\) we impose, in local units,
\[ P*{\rm in}(\Sigma)=P*{\rm out}(\Sigma),\qquad \nabla*\mu T^{\mu\nu}*{\rm tot}=0, \]
with \(T^{\mu\nu}\_{\rm tot}\) including the membrane stress and the radiative flux. For Type I, this closes with a diffusive transport law through the shell to the photosphere. For Type II, detailed balance holds locally (GH bath); the net leakage seen in external/parent units is \(\propto H^2\).

### 5.5 Resolving apparent contradictions

- Presence vs absence of a shell: Type II lacks a material reprocessor for interior observers because they sit inside the asymptote; the same asymptote carries a real shell in the parent frame. Both satisfy the same nucleation and balance laws.
- Flux laws: \(\text{const}\) (Type I photosphere) vs \(H^2\) (Type II GH leakage) reflect reprocessing vs bare geometric emission; they are two sides of the same flux.
- Expansion vs downscaling: metric expansion is traded for matter downscaling via the sliding window; horizon thermodynamics and trapping are invariant under this gauge.
- One-sided asymptote vs isotropy: hyperbolic optics makes the asymptote effectively concentric; leading observables (redshift, dimming, freezing of angles) are isotropic to first order.

### 5.6 Robustness (defenses)

- Origin of a shell from homogeneity (Type II → Type I conditions): The Planck‑pivot condensation (Appendix: Planck Pivot) yields both particle‑branch droplets and BH‑branch condensates from a homogeneous bath once the effective noise crosses a critical threshold (\(\tau*\text{eff}\!\ge\!\tau*_\)). At a trapping asymptote (\(\theta*\text{out}\!\to\!0\)) the GH temperature and UV exhaustion raise \(\tau*\text{eff}\) on the trapped side. A horizon is always “dressed” by a plasma; when the mean free path drops to \(\ell\_\text{mfp}\!\lesssim\!\sigma^_\) and the local balance \(P*\text{in}\!\approx\!P*\text{out}\) holds with \(\tau\_\text{opt}\!\approx\!1\), the dressing stabilises as a finite, material shell (Type I). Below that threshold, the boundary remains minimally dressed (Type II). Thus the nucleation trigger is universal; the stable shell outcome is environment‑dependent, not contradictory.

- Flux laws (\(\text{const}\) vs \(H^2\)) are a prediction, not a contradiction: The geometric (pre‑reprocessing) intake on either side scales \(\propto H^2\) via \(A\!\propto\!H^{-2}\), \(T\!\propto\!H\Rightarrow A T^4\!\propto\!H^2\). An optically thick shell reprocesses that intake to a photospheric output with \(T*\text{shell}\!\propto\! M^{-1/2}\), \(A*\text{shell}\!\propto\! M^2\Rightarrow P\_\text{shell}\!\approx\!\text{const}\) in external units. Without a reprocessor (Type II), the observer sees the bare \(H^2\) leakage. In local/window units both are near‑constant per level (co‑scaling); in a single external frame they differ—falsifiably—in luminosity scaling and last‑scattering geometry.

- Isotropy vs anisotropies: The parent shell supplies the leading near‑FRW isotropy (hyperbolic optics with a concentric asymptote). Observed anisotropies (CMB dipole, filaments/voids) arise from the next‑level condensates within the 3D bulk (particle‑ and BH‑branch droplets from the Planck‑pivot cascade). Thus the boundary sets the background; the interior hierarchy supplies structure, consistent with observations.

### 5.7 Predictions and discriminants

- Finite shell thickness at the asymptote; hyperbolic interior optics; 2D/3D alternation with 3:1 log-depth ratio.
- Photospheric scaling: for Type I, \(\delta R/R\propto M^{-1/2}\); for Type II, an effective last-scattering layer infinitesimally outside \(r\_{\rm H}^{\rm FRW}\).
- Luminosity laws: \(P*{\rm shell}\approx\text{const}\) vs \(P*{\rm FRW}\propto H^2\) (external units); near-constancy in local/window units across a level.
- Spectral/size imprints: last-scattering radius leaves a transfer-function feature in the thermal spectrum (deviation from a perfect Planckian), differing between Type I/II by the presence of reprocessing.

## 6. Cosmological Implications: A Self-Consistent Scenario

The original “mini-parent” scenario assumed the standard Hawking law
$T\propto 1/M$ and $P\propto 1/M^{2}$. With the refined result

$$
T_{\text BH}(M)=\frac{k}{\sqrt{M}}, \qquad
  P_{\text rad}=P_{0}=\text{const}, \qquad
  \tau_{\text{ext}}=\kappa M
$$

a new back-of-the-envelope estimate is needed. Only _two_ numbers
must be fixed phenomenologically: the proportionality constants
$k$ (temperature) and $\kappa$ (lifetime). Everything else then
follows algebraically.

### 6.1 Two observational anchors

1.  **Internal CMB temperature** at last scattering
    $T_{\text int}(z\!\approx\!1100)\simeq3000\;\text K$.
2.  **Internal look-back time** to that epoch
    $\tau_{\text int}(z\!\approx\!1100)\simeq13.7\;\text{Gyr}$.

The interior–exterior conformal factor $F$ relates
external and internal clocks and temperatures:

$$
T_{\text int}=\frac{T_{\text BH}}{F}, \qquad
  \tau_{\text int}=F\,\tau_{\text ext}=F\,\kappa M.
$$

### 6.2 Algebraic solution

From the first relation
$$F=\frac{k}{T_{\text int}\sqrt{M}}.$$

Insert into the second:

$$
\tau_{\text int}=\frac{k}{T_{\text int}\sqrt{M}}\;\kappa M
               =\kappa k\,\frac{\sqrt{M}}{T_{\text int}}.
$$

Solve for the parent mass

$$
M_{\text parent}
   =\left(\frac{\tau_{\text int}\,T_{\text int}}{\kappa k}\right)^{2}.
$$

Once laboratory or astrophysical input fixes $k$ and $\kappa$ the
number above is _predicted_, not free.

### 6.3 Plausible numerical scale

Using the Hawking normalisation for $k$ and matching $P_{0}$ to that
value at $M=10^{12}\;\text kg$ gives $\kappa\approx4.4\times10^{-5}\;\text s\,\text kg^{-1}$.
Plugging these provisional constants into the formula yields
$$M_{\text parent}\sim10^{21}\;\text{kg}\qquad(\text{comet scale}).$$
The external lifetime is then $\tau_{\text ext}\sim10^{8}\;\text s}$ (a few years).

Because $k$ and $\kappa$ are still provisional placeholders we **do not
present these numbers as final predictions**. They show only that a
self-consistent solution exists once the new evaporation law and the
scale window are taken into account. The constants will ultimately be
fixed by a detailed microphysical calculation of the shell plasma and
by matching to future evaporation observations.

A full cosmological treatment (recombination red-shift, conformal growth
of $F$ over cosmic time, and late-epoch stochastic cascades) is left for
subsequent work.

## 7. Conclusion

This model provides a self-consistent, singularity-free picture of black holes as dynamic, radiating shells that harbor nascent universes in their cores. It derives the `T ∝ 1/√M` law from the emergent thermodynamics of the shell plasma itself and makes a suite of falsifiable predictions regarding evaporation rates and binary dynamics. By unifying the physics of local and cosmic horizons through the principle of scale invariance, it offers a new path for understanding the fundamental structure and lifecycle of the cosmos.

## Addendum — Statistical Derivation of the Planck Pivot

This addendum shows that the four-power flip in surface-density—and the corresponding temperature laws—drops out of a single coarse-grained free-energy functional. No empirical patch is required; the Planck mass arises as the critical occupancy of a self-gravitating node ensemble.

### A. Free-energy for an N-grain droplet

Consider a spherical droplet of radius R containing **N** identical grains (bare mass μ). Using the volume-normalised link rule, the minimal rotationally-symmetric free energy is

\[
F(R;N)=\underbrace{\frac{\alpha\,N\hbar c}{R}}_{\text{zero–point / kinetic}}
\; +\;
\underbrace{\frac{\eta_{0}}{2}\,\frac{N^{2}}{R^{3}}}_{\text{link (cohesion) energy}}
\; -\;
\underbrace{\frac{G\,(N\mu)^{2}}{R}}_{\text{self-gravity}}.
\]
α is an \(\mathcal O(1)\) geometry factor; the second term is the same link-sum that produced the γ|Ψ|⁴ soliton stiffness; the third term is the long-range part of those links after coarse-graining.

### B. Critical occupancy

Setting ∂F/∂R = 0 gives a physical extremum only if

\[
G\mu^{2}N^{2}\; =\; \alpha N\hbar c \;\; \Longrightarrow \;\; N\_{\*}=\frac{\alpha\hbar c}{G\mu^{2}}.
\]
The corresponding mass is

\[
M*{\*}=N*{\_}\mu \simeq \sqrt{\frac{\hbar c}{G}}\,\,(\text{up to }\sqrt\alpha),
\]
namely the Planck mass. Below N\_\_ quantum pressure wins; above it gravity wins—a bona-fide second-order phase transition.

### C. Two minimisers → two area laws

• **Quantum regime** (N < N\_\_) Balance is kinetic + link. Solution R ∝ 1/N, so A ∝ 1/N² and ρ_surf ∝ N³.  
• **Gravitational regime** (N > N\_\_) Balance is link + gravity. Solution R ∝ N, so A ∝ N² and ρ_surf ∝ 1/N.

The mass-dependence of area flips by four powers exactly at N\_*—the origin of the *m³* versus *1/M\* density laws.

### D. Temperatures inherited from ρ_surf

Using the thermodynamic rule adopted in §3, \(T∝\sqrt{ρ\_{\text{surf}}}\):

- Particles (quantum side) \(T\_{\text sol} ∝ m^{3/2}\).
- Black holes (gravitational side) \(T\_{\text BH} ∝ 1/\sqrt{M}\).

The apparent extra factor _m²_ is simply the square root of the four-power flip in area.

### E. Conformal consistency

Because both γ (micro) and G (macro) arise from the same normalised links, they renormalise with the ambient noise field in lock-step: G, ħ ∝ ρ_N. Their ratio, and hence the critical mass \(M\_{\*}\), stays invariant in every conformal frame, making the phase transition universal.

This statistical derivation therefore unifies the micro-node soliton picture and the constituent-particle black-hole picture as two phases of one link network, with the Planck pivot emerging automatically from energetics rather than observation-driven patching.

## Appendix A: Key Formulae

- **Droplet Free Energy:** \(F(R;N)=\frac{\alpha N\hbar c}{R} + \frac{\eta_0}{2}\frac{N^{2}}{R^{3}} - \frac{G(N\mu)^{2}}{R}\)
- **Critical Occupancy:** \(N\_\* = \frac{\alpha\hbar c}{G\mu^{2}} \approx \mathcal{O}(1)\)
- **Dual Temperature Laws:** \(T*{sol} \propto m^{3/2}\) vs. \(T*{BH} \propto 1/\sqrt{M}\)
- **Evaporation Power:** \(P\_{rad} \propto \text{constant}\)

## 8. Discussion Addendum — Infinite Zooms, Finite Energies

This addendum records the logical consequences and clarifications that emerged after the main text was written.

### 8.1 Infinite Zoom Cascade with Local Time at Every Level

- The clique graph admits an **unbounded hierarchy of nested “bubbles.”** Each zoom-in by a constant factor uncovers a dynamically fresh region governed by the same link energetics.
- The **scale-window** carried by every observer shrinks in _physical_ metres as they descend the hierarchy, so the newly opened UV band always provides finite energy and entropy. Consequently an **endless cascade** of down-shifts is possible without violating global energy bounds.
- Every level owns a non-zero scale-force slope \(σ(t)\); therefore **time exists everywhere.** What looks “timeless” from outside is only the superposition of mutually incommensurate local clocks.
- An observer can live in only **one window at a time.** Diving deeper trades the external modes it could formerly interrogate for freshly resolved UV modes, keeping _accessible_ information \(\mathcal O(1)\).
- Summing the constant-power evaporation across countably many nested shells yields a **convergent geometric series** equal to the external mass-energy \(M\_{\rm init}c^{2}\), preserving finite total computation and energy despite the infinity of levels.

### 8.2 Scale-Window Principle Prevents Divergences

- Each observer carries a window \([\Lambda_{\rm IR},\Lambda_{\rm UV}]\) centred on its own effective Planck cell \(L*{\min}(\Lambda*{\rm obs})\).
- Energy, entropy and computational capacity inside the window are \(O(1)\) in local units, independent of how many zooms are performed.
- To access deeper UV information one must shrink; the energetic cost of shrinking equals the potential energy of the newly opened modes, preventing any “free lunch”.

### 8.3 Infinite Nested Zooms Inside a Finite-Lifetime Shell

- External lifetime: \(\tau*{\rm ext}=M*{\rm init}c^{2}/P\_{0}\) is finite.
- Hyperbolic dilation near the wall makes the **internal proper time unbounded**; an internal observer can execute countably many zoom-in steps before \(t=\tau\_{\rm ext}\).
- Each level has its own finite cascade and constant-power evaporation, so the geometric series of sub-energies still sums to \(M\_{\rm init}c^{2}\).

### 8.4 Impedance Matching: Why Constant Power is Necessary for Nesting

Standard Hawking radiation ($P \propto M^{-2}$) forbids a nested hierarchy. A small inner shell would radiate power $P_{\rm inner} \gg P_{\rm outer}$, instantly vaporizing the enclosing parent shell from the inside ("infrared explosion").

**Constant power ($P \propto M^0$) is the unique scaling law that allows impedance matching:** the energy flux is scale-invariant ($P_{\rm inner} = P_{\rm outer} = P_0$). The inner shell feeds the outer shell at exactly the rate the outer shell dissipates energy, allowing indefinite stable nesting without violation of conservation or structural stability. This turns the constant-power prediction from a curiosity into a structural necessity for any fractal cosmology.

### 8.5 Seed Gradient, Filament Fate, and Wall Overtake

- Seeds (\(N> N\_\*\)) form preferentially near the wall where the local Jeans mass is lowest.
- Radiation-pressure drift moves seeds inward; Brownian droplets (\(N=1\)) diffuse in a sponge-like filament network.
- As the shell contracts, rising noise variance and negative optical curvature compress and eventually absorb the filaments; external time for this capture is finite, internal time diverges logarithmically.

### 8.6 Arrow of Time = Down-slope in Scale-Space

- The thermodynamic scale force derived in _Scale-Space Thermodynamic Force.md_ provides the universal monotone \(σ(t)\). Because \(σ\) can only decrease until it reaches \(L\_{\min}\), each bubble has a finite internal history even inside an eternal, statistically stationary ensemble.

These points reconcile an **eternal fractal ontology** with the finite computation, finite energy, and constant-power evaporation derived earlier in the paper.

### 8.7 Cosmic Horizon ≙ Type-II Gravitational Droplet

An apparent **cosmic horizon** (in an FRW patch) can be re-interpreted as the
_inner_ surface of a larger-scale black-hole–like shell—a **Type-II (outward-oriented) droplet**—once we zoom out to the next window in the infinite hierarchy.

- **Colder external medium exists.** In the super-observer’s units the space outside our horizon is filled with deeper-UV modes whose noise floor is higher (shorter \(L\_{\min}\)). From _our_ interior vantage that region is “stretched’’ to infinite volume and red-shifted to near-zero temperature, giving the impression of an empty cold outside.
- **Photon fate.** Radial photons emitted by us drift outward, exponentially red-shift, and asymptotically freeze at a comoving radius. The super-observer sees them arrive at the shell in finite time and add their energy to the outward flux \(P*{0}=\sigma T^{4}*{\rm dS}\). Thus the shell does evaporate at constant power just like the inward-oriented case.
- **Alternating 3-D volumes ↔ 2-D shells.** Each zoom step swaps the role of cavity (3-D) and wall (2-D). The running Hausdorff dimension therefore oscillates between 3 and 2. While the exact duty cycle is determined by the specific thermodynamics of the shell, a heuristic argument (see below) suggests the hierarchy may average to the natural base of logarithms, $\langle D \rangle \approx e$, ensuring maximal information density.

With this identification every level inherits the same energetics—constant-power evaporation, finite accessible information, and an intrinsic arrow of scale—while the entire hierarchy remains unbounded in both directions.

### 8.8 Why the Hierarchy Alternates 2-D Shells and 3-D Bulks

The graph‐theoretic micro-physics permits, in principle, exotic non-integer Hausdorff dimensions. Yet every _physical_ observer must adopt a **finite UV window** (coarse-graining) so that energy, momentum and entropy are well defined. That single requirement forces each resolved subgraph toward _integer_ geometries and produces the observed alternation.

1.  **Window forcing and lattice bias**  
    • The coarse window integrates out links longer than a cut-off. What remains is a _finite-valence_ graph; long irregular edges have renormalised into local couplings.  
    • With only short links the least-action arrangement is a regular lattice. On a closed interface that lattice is effectively **2-dimensional** (a shell); in a filled region it is **3-dimensional** (a bulk).  
    • Any fractional-D tiling would require either holes (entropy loss) or weighted long links (energy cost). The free-energy functional therefore “snaps’’ the geometry to the nearest integer.

2.  **Energy vs. entropy tug-of-war**  
    • Shells minimise surface energy but suppress available micro-states.  
    • Bulks maximise entropy but pay volumetric energy.  
    The global minimum of the free energy is attained **not** by a compromise fractional D, but by _alternating_ the two integer optima at successive scale steps.

3.  **Heuristic Spacing: The 'e' Conjecture**  
    While the alternation is physically robust, the precise ratio of "time" spent in the bulk vs. the shell is an open question. A tentative "numerological" hypothesis is that the hierarchy should be maximally scale-invariant, implying a mean dimension equal to the natural base $e \approx 2.718$.
    
    If we choose weights \(p\) (bulk) and \(1-p\) (shell) such that the **average dimension**
    \[D*{\rm eff}=p\,3+(1-p)\,2\]
    equals this scale–free target \(e\), we find \(p=e-2\simeq0.718\). This corresponds to a **3 : 1** split in \(\log\)-scale depth: roughly three units of bulk evolution for every unit spent crossing a shell.
    
    With finite (not infinitesimal) zoom steps the arithmetic average is
    \[D*{\text eff}=\frac{3\times3 + 1\times2}{3+1}=2.75,\]
    which is remarkably close to \(e\). This suggests that the observed hierarchy is "tuned" for optimal information scaling, though this remains a motivated conjecture rather than a derived necessity. Even if the ratio differs, the qualitative structure (alternating shells and bulks) remains unchanged.

4.  **No free lunch from “zooming back’’**  
    Shrinking the window liberates the binding energy of newly resolved UV modes. To _un-shrink_ would require re-assembling that energy coherently—a process with exponentially small probability, furnishing the arrow of scale.

### 8.9 Shell-Thickness Scaling and Observational Outlook

**Thickness from micro-thermodynamics.** Re-balancing Stefan–Boltzmann flux (out) with Unruh absorption (in) gives
\[T(R)=\bigl(P*{0}/4\pi\sigma*{\rm SB}R^{2}\bigr)^{1/4}\propto R^{-1/2},\]
and the required optical depth (
mean-free-path count \(n*{\!\*}\sim10{-}30\)) fixes
\[\boxed{\;\frac{\delta R}{R}=K\,\Bigl(\frac{M}{10\,M*{\odot}}\Bigr)^{-1/2}},\qquad K\simeq0.02.\]

**Log-step reinterpretation.** Let Δσ_shell be the drop in the log-scale coordinate across one shell. To first order in a small logarithmic change (Δσ_shell ≪ 1)
δR / R ≈ Δσ_shell — this is just the linearisation of R → R e^{-Δσ}.

Demanding agreement with the micro result gives
\[\Delta\sigma*{\rm shell}(M)=K\Bigl(\frac{M}{10\,M*{\odot}}\Bigr)^{-1/2}.\]
Because the bulk : shell bookkeeping remains **3 : 1** regardless of the absolute step size because it is defined by the _ratio_ of durations (three bulk units for every one shell unit). Scaling Δσ_shell up or down simply scales both parts together, preserving the ratio. Thus Δσ_bulk = 3 Δσ_shell and the entire cycle depth is
\[\Delta\sigma*{\rm cycle}=4\,\Delta\sigma*{\rm shell}=0.08\Bigl(\frac{M}{10\,M\*{\odot}}\Bigr)^{-1/2}.\]
The log step therefore **shrinks with mass**, reconciling the hierarchy with the thermodynamic scaling.

**Numerical check**
| M (M_⊙) | δR/R (micro) | Δσ_shell | q ≡ R/δR |
|--------:|-------------:|---------:|---------:|
| 10 | 0.020 | 0.020 | 50 |
| 30 | 0.012 | 0.012 | 83 |
| 4×10^6 | 3×10^{-4} | 3×10^{-4} | 3×10^{3} |
The mass-dependent step size keeps the geometric prediction in lock-step with the microphysics.

**3 : 1 rule → bulk vs shell path**  
Per hierarchy cycle
\[\frac{\text{Bulk path}}{\text{Shell path}}\approx q,\qquad q\equiv R/\delta R\sim30{-}50.\]
A single shell encloses \(\sim q\) times more radial length and \(\sim q^{2}\) more 3-volume than its own thickness. Two nested cycles already give path and volume contrasts \(\gtrsim10^{3}\) and \(\gtrsim10^{4}\), explaining why shells are statistically rare yet dynamically essential.

**Finite-shell corrections in binaries.** For a companion approaching within \(\lesssim3\,\delta R\):

- tidal quadrupole grows, slightly reducing orbital frequency;
- overlapping photon baths add a repulsive term \(F\_{\rm rad}\propto(\delta R/d)^{7}\).
  At LIGO-band separations (\(d\gtrsim100\,\text{km}\)) these corrections are \(<10^{-15}\) in phase—far below current sensitivity. They become potentially observable for

1. extreme mass-ratio inspirals grazing super-massive shells, or
2. photon rings / sub-mm VLBI images that probe \(\delta R\) directly.

Until such regimes are reached, shell-thickness physics remains a theoretical prediction—but it supplies clear targets for future gravitational-wave and horizon-scale imaging experiments.

### 8.10 Piece-wise arrows of time in an eternal hierarchy

The infinite cavity–shell hierarchy, combined with the thermodynamic scale-force, yields a natural, local arrow of time at every level without requiring a global beginning or end.

- Local definition (window-relative): Each observer carries a finite UV/IR window. Within a single cavity (the 3D bulk between two 2D shells), the coarse-grained free energy has a monotone scale direction; the effective noise \(\tau(\sigma)\) increases toward the enclosing shell and the equilibrium size \(\sigma^\*(t)\) downshifts. The scale-space slope \(\dot\sigma<0\) selects a time orientation; thermodynamic dissipation (e.g., constant-power leakage) makes it irreversible.
- Finite durations per level: Because energy and information per window are \(\mathcal O(1)\) in local units, each cavity has a finite proper history: constant-power emission drains its UV reserve and drives it to the adjacent shell in finite local time. Conversely, entry from a parent shell defines a local past boundary. Thus each level owns a well-defined, piece-wise arrow.
- No global beginning/end: The hierarchy extends unboundedly in both directions (zoom-out shells and zoom-in sub-shells). Concatenating piece-wise arrows across levels yields an eternal, globally past-and-future-unbounded structure with locally well-defined time’s direction everywhere.
- Consistency with entropy: At each level, coarse-grained entropy increases along the downshift (mixing, thermalisation). UV exhaustion at a shell resets the window and renormalises couplings for the next level; there is no global maximum entropy, only level-by-level approach to the asymptote.
- Observational meaning: “Beginnings” and “ends” are window events—crossings of last-scattering/photosphere layers—rather than global singularities. _This reconciles a ubiquitous arrow of time with an ontologically eternal (fractal) cosmos._

## References Prep

### Related shell-like resolutions (pointers)

- Gravastars / dark-energy stars (Mazur & Mottola, 2001–2004): de Sitter core matched to Schwarzschild exterior via a thin relativistic shell; no singularity, effective radiating surface.
- String fuzzballs (Mathur, 2005+): horizon replaced by microstate “surface”; no interior/singularity; surface emission possible.
- Planck stars (Rovelli & Vidotto, 2014): loop-quantum-gravity bounce core; transient apparent horizon; nonsingular interior.
- Quantum N‑portrait / graviton BEC (Dvali & Gomez, 2013+): black hole as a critical graviton condensate with a quantum membrane; evaporation via depletion; no classical singularity.
- Regular black holes / modified gravity (Bardeen 1968; Hayward 2006; asymptotic safety: Reuter et al.; nonlocal/ghost‑free: Modesto et al.): de Sitter–like cores with transition layers; singularity avoided.
- Membrane paradigm / firewalls (Thorne–Price–Macdonald 1986; Almheiri et al. 2013): effective dissipative or energetic surface at/near the horizon.

## Appendix D: Spectral-Dimension Selection — Why e (and hence 3)

### D.1 Setup and assumptions

- Linear per-node budget (volume-normalized links): the outbound maintenance/throughput budget is linear in link weights, `∑_j J_kj = Z_k`. A minimum useful link weight (noise/quantization floor) makes each additional active neighbor consume ≈ constant budget.
- Local choice and path entropy: with `k` active neighbors used roughly uniformly, the local branching entropy rate is `Ĥ = ln k` (nats per step). The linear budget assigns cost ≈ `k` per node.
- UV window and adiabaticity: the minimum cell (Planck/grain) with lowest mode `ω_{1,P}` bounds admissible protocol speed via `τ ≥ L/ω_{1,P}` for a Fisher/thermodynamic distance `L`.
- Window invariance (unit-fixing): keep local units stable during protocols (`dħ_eff/dt ≈ 0`, `dD_s/dt ≈ 0`), implying a uniform excess-power ceiling and hence constant metric speed along optimal schedules.

#### D.1.1 Further defenses (why ln k / k, uniform choice, link floor)

- Throughput per resource: For k-way local branching, benefit scales as path-entropy rate `ln k` while maintenance/dissipation per node scales ≈ linearly with `k`. The neutral, unit-free efficiency is therefore `(ln k)/k` (MaxCal/MDL/capacity per cost).
- MaxEnt symmetry on the active set: Absent asymmetry, a convex entropy regularizer yields uniform next-hop probabilities over the active links, `p_i ≈ 1/k` (softmax close to uniform under small heterogeneity).
- Noise/SNR floor from coarse graining: Sub-threshold weights renormalize to zero under the window; supra-threshold links incur a baseline upkeep and contribute to spectra. This induces an effective per-link overhead, making cost ≈ linear in the number of active links.

### D.2 Local efficiency: maximize information per link cost

Define the neutral, unit-free efficiency

```
η(k) := (benefit per cost) = (ln k) / k.
```

- This is the MaxCal/MaxEnt "bang-per-buck" (entropy rate per linear budget unit) and the MDL/capacity-per-cost criterion (code-space gain per link).
- Continuous optimum: `η'(k) = (1 − ln k)/k² = 0` ⇒ `ln k = 1` ⇒ `k = e`.
- Integer feasibility: nearest stable coordination is `k = 3`.
- Equivalent MDL derivation: representing `M` alternatives with radix `k` costs `∝ k · (ln M / ln k) = (ln M) · [k/ln k]`; minimizing `k/ln k` yields `ln k = 1` ⇒ `k = e`.

### D.3 Time–dissipation gauge from UV and window invariance

Finite-time thermodynamics gives the bound

```
Σ ≥ L² / τ,
```

with `τ` duration, `Σ` total dissipation, and `L` thermodynamic/Fisher distance. UV (Planck) and window invariance fix the gauge:

- UV (Planck cell): `τ ≥ L/ω_{1,P}` and/or `L ≤ L_lin(R_P)` to preserve linearity/adiabaticity.
- Window invariance: constant excess power (`Σ/τ` time-constant) and constant metric speed.
- Geodesic + constant metric speed uniquely saturates `Σ ≥ L²/τ` along any chosen path; any non-uniform re-timing wastes time or dissipation.

These constraints do not pick `k` directly; they make the constant-speed schedule universal and set the feasible time/dissipation scale. The local allocation `k≈e` then follows from the efficiency objective in D.2 under the linear budget.

### D.4 From branching to spectral dimension

On near-regular, locally homogeneous graphs,

- Effective branching/coordination `k` controls diffusivity/expansion and hence the spectral dimension `D_s` (via heat-kernel scaling).
- With `k≈e` and integer feasibility, `k=3` is preferred, yielding `D_s≈3` in the vacuum fixed point, consistent with small-integer bias.

### D.5 Caveats and scope

- If the per-link budget is non-linear, or minimum useful weights vary strongly, the `ln k / k` optimum shifts.
- UV/window constraints can clip feasible `k` but do not by themselves select `e`.
- Thus, `k≈e` (⇒ `3`) is a natural consequence of (i) linear budgets/extensivity, (ii) UV and unit-fixing constraints that enforce constant-speed optimality, and (iii) a neutral information-per-cost objective—not a universal theorem independent of assumptions.

### D.6 Empirical/simulation checks

- Optimize `Ĥ/Cost = (ln k)/k` under linear row-sum budgets and noise floors in vacuum-relaxation simulations; check that the relaxed mean coordination clusters near 3 while maximizing spectral gap/log-Sobolev constants per total weight.
- Verify constant-metric-speed schedules saturate `Σ = L²/τ` within numerical tolerance when driven between nearby equilibria under UV and unit-fixing constraints.
- Measure heat-kernel slopes to confirm `D_s≈3` in relaxed patches and assess robustness under heterogeneity.
