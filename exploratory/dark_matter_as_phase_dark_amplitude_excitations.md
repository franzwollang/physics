# Dark matter as phase-dark amplitude excitations (stability + non-conversion)

This note sketches a conservative “next branch” after the MOND/kinematic-sampling dead end: **keep the Newtonian/Poisson window as real**, and explain galactic missing mass via an additional population of **gravitating but electromagnetically quiet excitations** of the *same substrate*.

The goal is to stay aligned with what the drafts already commit to:

- Gravity paper: a single order parameter \(\Psi=w e^{i\phi}\) with free-energy \(\mathcal F[w,\phi]\), and two generic low-energy sectors: **gapped amplitude** \(w\) and **massless phase** \(\phi\).
- ICG SI: amplitude is gapped (finite core energy, screening length \(\ell\)); phase is Goldstone (massless) with a shift symmetry; defects require \(w\to 0\) locally.
- Gravity SI: a U(1)-gauge-like structure is associated with **phase defects/sheets**, and “matter–kernel coupling” gives Poisson equivalence \(\nabla^2\delta\tau^2=-C\rho_m\).

---

## 1) What “phase-only”, “amplitude-only”, and “mixed” mean in *your* equations

### 1.1 Starting point: the coupled Euler–Lagrange equations

From the Gravity paper appendix (free-energy solutions), the stationary equations are
\[
-\alpha_{\rm grad}\,\nabla^2 w + V'(w) + \frac{\kappa}{2}\,|\nabla\phi|^2\,w = 0,\qquad
\nabla\cdot(\kappa w^2\nabla\phi)=0,
\]
with a Mexican-hat potential
\[
V(w)=-\beta_{\rm pot}w^2+\gamma w^4,\qquad w_*^2=\beta_{\rm pot}/(2\gamma).
\]

So the two fields are coupled by the **phase-gradient pressure term** \((\kappa/2)|\nabla\phi|^2 w\) and by the \(w^2\) “conductivity” in the phase equation.

### 1.2 Phase-only excitations

In the phase-only regime you clamp \(w\approx w_*\) (valid on scales \(\gg\ell\)), leaving
\[
\mathcal F_\phi\approx \frac{\kappa w_*^2}{2}|\nabla\phi|^2
\]
and, dynamically, \(\phi\) obeys the massless wave equation with speed \(c_s^2=\kappa w_*^2\).

**Two subtypes exist:**

- **Radiation (smooth phase waves):** localized wave packets, not static lumps.
- **Topological phase defects:** phase windings/sheets, which *force* \(w\to 0\) somewhere (a core) so the phase becomes undefined.

That second subtype is crucial: “phase-only” in the strict sense cannot carry nontrivial winding without an amplitude core.

### 1.3 Amplitude-only excitations

Amplitude-only means the phase is (approximately) constant in space on the lump scale:
\[
\nabla\phi\approx 0\quad\Rightarrow\quad \nabla\cdot(\kappa w^2\nabla\phi)=0\ \text{trivially}.
\]
Then the amplitude equation reduces to a gapped scalar equation (plus whatever nonlocal cohesion you include in the scale-space reduction).

In your framework, amplitude-lump stability is not being attributed to a pure local \(|\nabla w|^2+V(w)\) model alone; it is attributed to the **boundary-dominated cohesion / shell interaction** in the scale-space derivation:
\[
E(\sigma,x)=A_\gamma\sigma^{-2}-\widetilde A_\eta f(\tau(x))\sigma^{-1},
\]
which yields a true size minimum \(\sigma^*(x)\) and a stable localized object.

So an “amplitude-only soliton” is simply a soliton core supported by the amplitude stiffness + short-range cohesion balance, with no significant phase texture.

#### 1.3.1 Downscaling universality as a configuration-space attractor statement (DM vs baryons)

If the **shell interaction picture is formulated in configuration space** (i.e. it is a coarse description of how localized excitations couple to their environment through a boundary layer / window), then **the object does not need to be literally spherical in real space**. What matters is that, after coarse-graining over internal microstates, the excitation has an effectively **isotropic attractor basin** in the reduced variables that control its environmental coupling.

Operationally, “DM downscales the same way as baryons” should mean:

- **Universal piece (same functional response):** the environment enters only through the same scalar \(f(\tau)\), so \(\sigma^*(x)\propto 1/f(\tau(x))\) for *any* massive soliton family (baryons, DM candidates).
- **Non-universal piece (different basins):** the profile-dependent constants \(A_\gamma,\widetilde A_\eta\) can differ between families (different internal textures/topologies), changing the absolute size/mass scale but not the **form** of the downscaling law.

This makes “DM downscales” a **robust consequence** of being an amplitude-core excitation, while allowing DM to be structurally distinct (e.g. knotted tubes, neutral composites): those details largely renormalize prefactors rather than changing the \(\tau\)-dependence.

**When can this fail?** If a candidate is strongly non-ergodic internally (very rigid topological structure, slow internal relaxation), then a single-\(\sigma\) adiabatic reduction can break locally and produce small sector-dependent corrections (hysteresis, extra slow modes). But the leading expectation remains: any amplitude-core object experiences downscaling in higher \(\tau\).

### 1.4 Mixed amplitude+phase excitations

Mixed excitations are amplitude solitons **with nontrivial phase structure** (currents, winding, sheets terminating on cores, etc.). In your papers, baryonic matter is naturally discussed as amplitude solitons coupled to a universal phase sector; it is therefore natural to view ordinary visible matter as “mixed” in this taxonomy.

---

## 2) What keeps each class stable?

### 2.1 Radiation (smooth phase waves)

- Protected by the Goldstone nature of \(\phi\): no mass term in the homogeneous vacuum.
- Not a static bound object; stability is kinematic (propagation at \(c_s\)), not a local minimum of a static energy.

### 2.2 Phase defects (winding/sheets)

Gravity SI makes the key point: changing winding requires \(w\to 0\) somewhere.

- **Topological protection:** \(\pi_1(U(1))=\mathbb Z\). Winding cannot be removed by continuous deformations of \(\phi\).
- **Energetic pinning:** the amplitude mode is gapped; forcing \(w\to 0\) over a region \(\sim\ell\) costs a finite core energy.

So defects can be long-lived because “unwinding” requires a nonperturbative core event.

### 2.3 Amplitude-only solitons

In your scale-space picture, stability is a literal energy minimum in \(\sigma\):

- \(\sigma^{-2}\) repulsion from gradient stiffness
- \(\sigma^{-1}\) cohesion from shell coupling (volume normalization + boundary interaction)

This produces \(\sigma^*(x)\) and therefore a stable lump.

The amplitude gap also implies:

- isolated amplitude disturbances have Yukawa tails with screening length \(\ell\), providing a natural core thickness and finite self-energy.

### 2.4 Mixed excitations

Mixed excitations inherit both stabilizers:

- the amplitude-lump \(\sigma\) minimum,
- and any topological constraints associated with their phase texture.

So “mixed” can be more stable than “amplitude-only” for given conserved/topological charges.

---

## 3) Why don’t they rapidly transition into each other?

This is the critical question for using “amplitude-only” as dark matter.

### 3.1 Phase-only \(\to\) amplitude soliton is suppressed by the amplitude gap

To convert phase-wave energy into a localized amplitude core you must excite \(w\) significantly away from \(w_*\). Because the amplitude mode is gapped (mass \(m_\xi^2=2\beta_{\rm pot}\)), this requires **local energy density above a threshold**.

In weak-field astrophysical environments, generic phase fluctuations are not expected to self-focus above the gap scale unless a specific nonlinear instability exists. Absent such an instability, the conversion rate is small.

### 3.2 Amplitude-only \(\to\) phase-defect requires creating phase texture *and* cores/sheets

To create a phase defect you need either:

- a winding configuration of \(\phi\), which costs gradient energy \(\propto w^2|\nabla\phi|^2\), and
- (if the winding is nontrivial) an amplitude core where \(w\to 0\) to regularize it.

This is not forbidden, but it is not entropically automatic either: it is a specific, structured excitation.

### 3.3 Mixed \(\leftrightarrow\) amplitude-only is blocked by topological barriers when winding is present

If the mixed excitation carries nonzero winding (or sheet endpoints), removing the phase structure requires a core-mediated event. That is an energetic barrier.

If the mixed excitation carries only smooth currents (no net winding), then dissipation/relaxation can in principle shed the phase texture as radiation; but in that case the mixed state is not protected and may relax to amplitude-only over long times.

This suggests a natural separation:

- **Visible-sector matter:** mixed states maintained by chemistry/phase-locked structures and interactions.
- **Dark-sector candidates:** amplitude-only (or minimally phase-textured) states that do not efficiently couple into the phase-defect network.

---

## 4) Why amplitude-only could be “dark” but still gravitate

Your Gravity SI’s U(1) gauge-like discussion ties gauge structure to **phase defects/sheets** (distributional \(\nabla\phi\) configurations). That suggests a plausible route:

- EM-like long-range interactions are associated with phase-sector defect/sheet structure.
- Objects with negligible phase texture (no winding, no sheet endpoints, minimal currents) couple weakly to that gauge channel.

But gravitational coupling in your framework is tied to:

- stress-energy (both \(T_w\) and \(T_\phi\) source curvature), and
- the matter–kernel coupling \(\delta K\propto\rho_m\) leading to Poisson equivalence \(\nabla^2\delta\tau^2=-C\rho_m\).

An amplitude-only lump has nontrivial \(w\) gradients and potential energy, so it has \(T_w\neq 0\) and thus acts as a gravitational source.

So the conceptual target is: **phase-dark amplitude lumps** that contribute to \(\rho_m\) (and thus \(\Phi,\delta\tau^2\)) but do not participate strongly in the phase-defect/gauge network.

---

## 5) What would we need to show (next derivations)

To make this a credible “dark matter from the substrate” mechanism, the next concrete derivations are:

1. **Existence of at least two metastable lump families** at fixed total norm/energy:
   - mixed (phase-textured / interacting)
   - phase-dark amplitude-only (weakly interacting except gravitationally)

2. **A barrier argument** that suppresses interconversion in the late universe:
   - amplitude gap + topological constraints + energy/momentum conservation.

3. **A production story**:
   - early-universe formation of amplitude-only lumps (e.g., during relaxation after symmetry breaking / defect annihilation cascades), leaving a relic abundance.

4. **Phenomenology**:
   - gravitational clustering (halos)
   - lensing consistency (same \(\Phi\) sources rays via \(\chi(\rho_N)\))
   - minimal non-gravitational coupling (no strong EM scattering)

---

## 6) Clarifications: decay channels and why “phase-dark” is not the same as “cannot radiate phase waves”

A useful distinction:

- **"Has no net long-range phase/gauge coupling"** (neutral at infinity)
- vs.
- **"Cannot exchange energy with the massless phase sector"** (radiatively decoupled)

Those are not the same.

### 6.1 Why an amplitude-dominated object can still radiate the phase sector

Even if a lump has \(\nabla\phi\approx 0\) in its *static* configuration, the full coupled equations show that time-dependent amplitude dynamics can source phase excitations.

From the Gravity paper appendix, the coupled stationary equations include the term
\(
(\kappa/2)|\nabla\phi|^2 w
\)
that couples phase gradients back into the amplitude equation, and the phase equation has an effective conductivity \(\propto w^2\):
\[
\nabla\cdot(\kappa w^2\nabla\phi)=0\quad(\text{static}).
\]
In dynamical settings (time-dependent \(w\)), the phase sector is not generically inert: changes in \(w\) modulate the phase stiffness and allow parametric emission of phase waves unless a symmetry/selection rule forbids it.

Moreover, Gravity SI explicitly uses a **gradient coupling** between a soliton internal coordinate \(q\) and the phase field:
\(
H_{\rm int}=-\lambda\,q\,e_i\,\partial_i\phi
\)
when deriving dissipative radiation/drag. That is a concrete statement that amplitude-sector degrees of freedom can radiate into the massless phase bath.

So: an “amplitude-only” lump is not automatically “dark.” It is only dark if it is **effectively decoupled** from the phase/gauge channel relevant for EM-like interactions.

### 6.2 What “dark” should mean in this framework

A conservative definition consistent with your Gravity SI gauge-like construction is:

- “bright” matter participates in the **phase-defect / sheet network** that yields EM-like long-range structure;
- “dark” matter does not carry the relevant long-range defect charge/current (or carries it only in neutralized combinations), so it has no strong long-range coupling in that channel.

This still allows **phase radiation** in principle; what is suppressed is strong *charged* interaction and efficient thermal coupling to the EM-like sector.

---

## 7) Making “mixed but phase-neutral” non-handwavy

You are right that “nontrivial internal phase but no net coupling” sounds like an enormous space of solutions. The way to make it solvable is to focus on a **small set of neutral composite topologies** whose existence is generic, and whose far field cancels by construction.

The correct analogy is not “some arbitrary neutral object” but **a neutral bound state built from charged constituents**, like atoms.

### 7.1 Minimal neutral composites (the ones worth focusing on first)

1. **Defect–antidefect bound pairs (dipoles):**
   - sheet endpoints carry opposite orientation/charge;
   - a bound pair has no net Gauss-law flux at infinity;
   - far field is dipolar (shorter range / weaker coupling).

2. **Closed sheets / closed loops (no endpoints):**
   - a sheet that closes on itself has no endpoints by definition;
   - no net enclosed charge for large surfaces;
   - still stores internal phase structure (tension/energy) and can be massive.

3. **Vortex–antivortex / winding-neutral composites:**
   - net winding at infinity is zero;
   - unwinding requires core events, but far-field cancels.

These are not arbitrary; they are the natural neutral objects implied by the phase-defect picture.

### 7.2 Why they can be *more* stable than pure amplitude lumps

Neutral composites can carry quasi-conserved labels that block decay pathways:

- topological linking/knotting of closed structures,
- separation barriers for defect–antidefect annihilation (core energy + geometry),
- angular momentum stored in internal phase currents.

This is the sense in which “mixed but neutral” can be a better DM candidate than “pure amplitude”:

- it can be *dark at long range* (neutral),
- yet metastable due to topological/energetic barriers.

### 7.3 Why this is still tractable

Even though the full solution space is large, you can make progress by:

- defining the **simplest neutral composite ansatz** (e.g., a defect–antidefect dipole with separation \(d\), or a closed loop of radius \(R\)),
- estimating its energy as a function of \(d\) or \(R\) (core energies + sheet tension + gradient energy),
- checking whether there is a metastable minimum.

This is exactly how one makes “many solutions” into a solvable program: start with the smallest neutral object and see if it has a robust metastable basin.

---

## 8) Taxonomy of Stability: Baryons vs. Dark Matter

This framework supports a clean topological distinction between ordinary baryonic matter and dark matter, based on the **nature of the closure** and the **mechanism of decay**. Both are "neutral composites" (zero net topological charge at infinity), but their stability arises from different sources.

### 8.1 The Baryon (Neutron): A "Composite of Partial Twists"

- **Structure:** A closure of three **fractional ribbons** (partial phase defects, e.g., \(2\pi/3\)). Each ribbon cannot exist alone (infinite sheet energy), so they *must* join to form a closed \(2\pi\) geometry.
- **Neutrality:** The three fractional windings sum to integer (0 mod 1), so there is no mismatched sheet extending to infinity.
- **Decay Mechanism (Fragile):** Because the object is built of **segments**, the vacuum allows the nucleation of a "spare part" (a phase-anchor dipole, \(\pm 1\)). A stray \(+1\) anchor can locally "swap in" to repair or re-bind one of the fractional ribbons.
  - *Mechanism:* Local Nucleation + Re-binding \(\to\) Identity change (Beta decay).
  - *Context:* Neutrons are unstable in vacuum because this "swap" is energetically accessible. They are stabilized inside nuclei because the environment changes the energy budget of the swap channel, not because the particle changes structure.

### 8.2 Dark Matter: A "Topological Knot"

- **Structure:** A **single continuous closed loop** (or knotted sheet) that has no "segments." It is a phase-defect loop (integer winding) that is twisted or knotted on itself (e.g., a Trefoil knot of the phase core).
- **Neutrality:** It has **no endpoints**. A closed loop has zero net flux through a large sphere. It is electromagnetically "dark" because it carries no net topological charge to infinity.
- **Stability Mechanism (Robust):** It cannot be destroyed by swapping a part because it has no parts. To destroy a knot, you must **cut the rope**.
  - *Mechanism:* "Cutting" the phase field requires forcing the amplitude \(w \to 0\) at the crossing point (reconnection).
  - *Protection:* In the cold, late-time universe, the energy density required to create a \(w \to 0\) core event is exponentially suppressed by the **Amplitude Gap** (\(V_{\rm pot}\)). The knot is "frozen" because the universe is too cold to melt the field lines.

### 8.3 Summary Table

| Feature | **Neutron (Baryon)** | **Dark Matter (Knot)** |
| :--- | :--- | :--- |
| **Topology** | **Composite** of 3 partial ribbons | **Single** knotted/linked loop |
| **Constraint** | "Must close to avoid infinite sheet" | "Cannot untie without cutting" |
| **Decay Mode** | **Local Swap:** Anchor nucleation changes parts | **Global Cut:** Reconnection requires core creation |
| **Stability** | **Metastable:** Vulnerable to "spare parts" | **Gap-Protected:** Frozen by \(V_{\rm pot}\) |

This categorization leverages the framework's **Geometric Grades** (ribbons vs. loops) and **Amplitude Gap** to naturally explain why one family is chemically active and metastable, while the other is inert and essentially stable on cosmic timescales.

---

## 9) Uniformity \u21d4 high-energy (precise): geometrogenesis, effective temperature, and why \u201czooming in\u201d does not make the vacuum \u201cyoung\u201d

The phrase \u201c**cold, late-time vacuum**\u201d should not be read as \u201cold at large scales only.\u201d In this framework, \u201cold vs. young\u201d is not primarily a statement about *length scale* but about **which term dominates the local coarse-grained free energy** and therefore what kinds of fluctuations are typical.

### 9.1 Two notions of \u201cuniformity\u201d exist, and only one is \u201chigh-energy\u201d

ICG derives vacuum organisation from the slow link functional (paper, Sec. \u201cEmergence of Geometric Structure and Locality\u201d). After marginalising fast phase variables, the vacuum relaxation functional for links is
\[
F_{\rm vac}[J] \;=\; \frac12\sum_{k<j}J_{kj}(\phi_k^\star-\phi_j^\star)^2 \;-\; H[J] \;+\; \beta(\langle d_{\rm int}\rangle[J]-d^\star)^2,
\]
with the convex entropy regulariser \(H[J]=\sum (J_{kj}/Z_k)\ln(J_{kj}/Z_k)\) and \(\phi^\star\) solving the graph Laplace equation. Stationarity yields a Boltzmann-like link form
\[
J_{kj}\propto \exp\!\Big[-\tfrac12\,\Theta_k(\phi_k-\phi_j)^2\Big],
\]
where \(\Theta_k\) plays the role of an effective temperature/constraint scale for the local outgoing weight budget.

From this, there are two distinct \u201cuniform\u201d regimes:

- **Hot-uniform (entropy-dominated / pre-geometric):** a maximally mixed, highly fluctuating regime where the system has not yet self-organised into a near-regular sparse backbone. In the Quantum-Graphity analogy explicitly cited in ICG, this corresponds to the \u201chigh-energy, maximally connected, non-local\u201d phase before geometrogenesis. In this regime, strong fluctuations can routinely access \(w\approx 0\) patches, and topology-changing reconnections are common.
- **Cold-uniform (ground-state / near-regular vacuum):** the late-time near-regular lattice-like state selected by vacuum relaxation: local isotropy holds (SI Laplacian convergence), \(w\approx w_\star\) almost everywhere, and the amplitude mode has a finite gap \(m_\xi^2=2\beta_{\rm pot}\) with screening length \(\ell\). This is \u201cuniform\u201d because it is close to a fixed point of \(F_{\rm vac}\), not because it is hot.

So: **uniformity \(\not\Rightarrow\) high energy in general**. \u201cHigh energy\u201d corresponds specifically to the **hot-uniform** phase where entropy/temperature-like control parameters make large fluctuations typical.

### 9.2 Fractal filamentarity as \u201cenergy condensation\u201d on top of a cold-uniform vacuum

The ICG picture already separates three layers:

1. **Vacuum backbone (near-regular, low valence):** selected by minimising \(F_{\rm vac}[J]\); this produces operational locality and a smooth Laplace\u2013Beltrami phase kernel on patches.
2. **Excitations/droplets (localized departures):** matter are localized amplitude departures \(w\neq w_\star\) with Yukawa tails; they are \u201cenergy condensed\u201d into small regions embedded in the vacuum.
3. **Large-scale clustering:** many droplets cluster under the long-range window (\(R_{\rm cl}\ll r\ll\ell\)), producing effective \(1/r\) far fields and the usual hierarchical structure.

The resulting universe is therefore naturally **filamentary/fractal at the level of where excitations sit**, while the intervening voids are \u201ccold-uniform vacuum\u201d in the sense above.

### 9.3 Why \u201czooming in\u201d does not find a \u201cyoung bath\u201d in a void

Zooming changes the *analysis window* and what you call \u201clocal,\u201d but it does not by itself move you from cold-uniform to hot-uniform. In a void, at any resolution above the microscopic cutoff, you keep finding:

- \(w\approx w_\star\) (vacuum), and
- suppressed amplitude excursions because the amplitude mode is gapped.

The place you *do* find something \u201cyoung-like\u201d (in the sense of \(w\to 0\) access) is inside **cores/defects**, i.e. inside matter/defect structures where \(w\) is already displaced. But those are precisely the localized, condensed excitations, not generic vacuum patches.

### 9.4 Reconnection/freeze-out statement in ICG variables (what \u201ccold\u201d means operationally)

The \u201cfrozen knot\u201d statement should be read as: in the cold-uniform vacuum, the probability that a patch of the substrate accidentally enters the \(w\approx 0\) sector required for reconnection is strongly suppressed. Operationally, for a topology-changing event you need a local amplitude-core excursion; its rate is of Kramers/Boltzmann form
\[
\Gamma_{\rm rec}\;\sim\;\Gamma_0\,\exp\!\Big(-\Delta E_{\rm core}/\Theta_{\rm loc}\Big),
\]
where \(\Delta E_{\rm core}\) is controlled by the amplitude potential stiffness (ultimately by \(\beta_{\rm pot},\gamma\)) and \(\Theta_{\rm loc}\) is the local effective noise/temperature scale for that degree of freedom (set by the bath/noise environment and window). In the hot-uniform regime, \(\Theta_{\rm loc}\) is large enough that \(\Gamma_{\rm rec}\) is appreciable; after geometrogenesis and subsequent cooling/relaxation, \(\Theta_{\rm loc}\ll \Delta E_{\rm core}\) in typical vacuum patches, and knots become long-lived \u201cfossils.\u201d

This gives the desired embedding: **\u201cyoung/high-energy\u201d is the regime where the entropy/temperature-like controls in the coarse vacuum functional make topology-changing excursions typical; \u201cold/low-energy\u201d is the regime where vacuum relaxation has produced a near-regular backbone and the amplitude gap makes such excursions exponentially rare in generic (void) regions.**

---

## 10) Minimal energy model for a stable knotted loop: tension + bending/twist \u21d2 \(R_\* > 0\)

This section sketches a minimal, solvable energy model showing how a **closed, knotted defect tube** can have a **finite metastable size** in this framework.

### 10.1 Effective description: a defect tube as an elastic string

The Forces/Particles picture treats gauge-like structure as carried by phase defects/sheets, while the ICG/Gravity free energy gives a gapped amplitude sector with a finite core thickness \(\ell\). A knotted DM candidate can therefore be idealized as a **thin tube** (core radius \(a\sim \ell\)) in which \(w\) is suppressed and the phase texture is concentrated.

At scales \(R \gg a\), its long-wavelength energetics can be parameterized as an elastic string/rod with:

- **Tension** \(T\) (energy per unit length): the cost to maintain a tube of suppressed \(w\) and phase gradients across its cross-section.
- **Bending modulus** \(B\): additional cost to curve the tube (curvature forces extra gradients/strain in the core).
- **Twist modulus** \(C\): cost to maintain internal phase twist/current along the tube (when a conserved/topological twist number is present).

This is standard \u201cthin-tube\u201d reasoning: you integrate the 3D field energy density over the cross-section and keep the leading invariants of the tube geometry.

### 10.2 The energy functional

Let the tube be a closed loop of length \(L\), with centerline curvature \(\kappa(s)\) and (optional) twist density \(\omega(s)\). A minimal free energy is
\[
E \;\approx\; E_{\rm core} \;+\; T\,L \;+\; \frac{B}{2}\int \kappa(s)^2\,ds \;+\; \frac{C}{2}\int \omega(s)^2\,ds.
\]

- \(E_{\rm core}\) collects non-extensive core formation energy (its detailed split does not affect the radius minimization below).
- The \(T\,L\) term drives shrinkage.
- The bending/twist terms grow as curvature/twist grows; for a shrinking loop, they can oppose collapse.

**Microscopic origin (scaling).** These coefficients ultimately descend from the same stiffnesses in the continuum free energy,
\[
\mathcal F[w,\phi] \sim \alpha_{\rm grad}|\nabla w|^2 + V(w) + \frac{\kappa}{2}w^2|\nabla\phi|^2,
\]
after integrating over a cross-section of radius \(a\sim \ell\). Dimensional estimates give \(T\sim \sigma_{\rm sheet}\) and \(B,C\sim T\,a^2\) (up to \(O(1)\) profile factors).

### 10.3 Circular loop estimate and metastable radius

To show existence of a finite optimum, it suffices to evaluate a circular loop of radius \(R\) (length \(L=2\pi R\)) with constant curvature \(\kappa=1/R\).

If the loop carries a conserved integer twist \(n\) (e.g., a trapped phase rotation along the tube consistent with a fixed linking number), then a simplest uniform-twist ansatz gives \(\omega \sim 2\pi n/L = n/R\).

Plugging these into the energy yields the canonical \u201c\(R + 1/R\)\u201d structure:
\[
E(R) \;\approx\; E_{\rm core} \;+\; 2\pi T R \;+\; \frac{B}{2}\,(2\pi R)\frac{1}{R^2} \;+\; \frac{C}{2}\,(2\pi R)\frac{n^2}{R^2}
\;=\; E_{\rm core} + 2\pi T R + \frac{\pi(B + C n^2)}{R}.
\]
Define \(A := \pi(B + C n^2)\). Then
\[
E(R) = E_{\rm core} + 2\pi T R + \frac{A}{R}.
\]
Minimizing,
\[
\frac{dE}{dR} = 2\pi T - \frac{A}{R^2} = 0
\quad\Rightarrow\quad
R_\* = \sqrt{\frac{A}{2\pi T}} = \sqrt{\frac{B + C n^2}{2T}}.
\]
Since \(T>0\) and \(B,C\ge 0\), any nonzero \(B\) (bending stiffness) or nonzero twist number \(n\) yields a finite \(R_\*>0\) and \(d^2E/dR^2 = 2A/R^3>0\), i.e. a **local minimum**.

**Interpretation.**

- Without bending/twist stiffness, a pure-tension string collapses.
- A defect tube in this framework is not an ideal Nambu string: it has a finite thickness \(a\sim \ell\) and therefore generically inherits bending/twist rigidity from the underlying \(|\nabla w|^2\) and \(w^2|\nabla\phi|^2\) costs.
- If the loop is also knotted/linked, topology further suppresses elimination because \u201cunknotting\u201d requires reconnection (a \(w\to 0\) event at a crossing), but that is a **separate longevity mechanism** from the existence of \(R_\*\).

### 10.4 What this model does and does not prove

This proves only a minimal point: **the effective field theory does not force closed tubes to shrink to zero size**; a finite-radius minimum is generic once you include the leading geometric stiffness terms that thin tubes inherit.

To make this quantitative (and to connect to a DM mass scale), the next steps would be to:

- derive \(T,B,C\) (even just scaling laws) from the underlying free-energy parameters \((\alpha_{\rm grad},\beta_{\rm pot},\gamma,\kappa)\) and the core profile \(w(r)\),
- estimate whether \(R_\*/a\) is large enough to justify the thin-tube approximation,
- and map the resulting mass \(M\sim E(R_\*)\) and interactions to halo phenomenology.

---

## 11) Analytical Sanity Checks & Simulation Strategy

Before committing to full nonlinear simulations, two "back-of-the-envelope" checks confirm that this candidate fits the coarse observational constraints (mass dominance and self-interaction limits).

### 11.1 The "Geometric Probability" Argument (Relic Ratios)

Why should Dark Matter (knots/loops) be 5$\times$ more abundant than Baryons (junctions)?

- **DM (Knotted Loops):** Requires a single flux tube to loop back on itself. In a "spaghetti" vacuum of random phase filaments formed during a phase transition, closed loops are the generic "garbage" left over when the network relaxes.
- **Baryons (3-Ribbon Junctions):** Requires the coordination of *three* specific fractional defects meeting at a point (or a closure of complex topology).
- **Prediction:** Purely combinatorial arguments (analogous to polymer physics or cosmic string networks) suggest that **simple loops should vastly outnumber complex junctions**.
- **Result:** This naturally fits the observed mass dominance of Dark Matter without fine-tuning. If the ratio were reversed (Baryons > DM), the theory would struggle to explain why the "complex" object is more common than the "simple" one.

### 11.2 The "Puffiness" Check (SIDM Constraints)

Dark Matter cannot be arbitrarily "fluffy" or it would violate self-interaction constraints (e.g., Bullet Cluster, halo sphericity). The astrophysical limit for Self-Interacting Dark Matter (SIDM) is roughly $\sigma_{\rm int}/m < 1 \text{ cm}^2/\text{g}$.

If we assume a DM knot has roughly nuclear mass ($M \sim 1$ GeV $\sim 10^{-24}$ g) and nuclear size ($R \sim 1$ fm $\sim 10^{-13}$ cm), its geometric cross-section is $\sigma_{\rm int} \approx \pi R^2 \sim 10^{-26}$ cm$^2$.

\[
\frac{\sigma_{\rm int}}{m} \sim \frac{10^{-26} \text{ cm}^2}{10^{-24} \text{ g}} \sim 0.01 \text{ cm}^2/\text{g}
\]

- **Result:** This lands safely within the allowed window ($\lesssim 1$ cm$^2$/g).
- **Implication:** If the solitons were much larger (e.g., atomic size) for the same mass, they would be ruled out. The fact that the amplitude gap scale $\ell$ (which sets the core size) is likely the nuclear/QCD scale naturally places these candidates in the correct "compactness" regime.

### 11.3 Simulation Roadmap: Related Fields & Procedures

To move beyond these estimates, numerical simulations should leverage methods from three established fields:

#### A. Lattice Field Theory / Scalar Electrodynamics

*For simulating the phase transition, core formation, and defect statistics.*

- **Model:** Discretize the complex scalar field $\Psi = w e^{i\phi}$ on a 3D lattice.
- **Procedure:** Initialize with random phases (high $T_{\rm eff}$), then cool/relax the system (quench) to observe defect formation.
- **Goal:** Measure the density of closed loops vs. junctions vs. monopoles.
- **Tools:** Standard algorithms for the XY model or Abelian Higgs model (e.g., Metropolis-Hastings, Hybrid Monte Carlo).

#### B. Cosmic String Network Evolution

*For understanding the long-term stability and scaling of the loop population.*

- **Context:** Cosmic strings are topological defects in gauge theories. While your defects are "fat" tubes, the scaling laws for loop chopping and network evolution are relevant.
- **Procedure:** Simulate the evolution of a network of strings in an expanding background (or effective cooling bath).
- **Goal:** Determine if the loop distribution reaches a "scaling solution" or if small loops rapidly decay. (Your "frozen knot" argument suggests they don't decay, unlike standard cosmic strings which radiate GWs, but the network dynamics are similar).

#### C. Polymer Physics / Knot Theory

*For analyzing the stability of specific knotted configurations.*

- **Context:** Polymers are modeled as elastic chains with excluded volume and topological constraints.
- **Procedure:** Model the defect tube as a discrete chain of beads with bending/twist stiffness (as in Section 10).
- **Goal:** Compute the energy landscape for different knot types (Trefoil, Figure-8) to verify the existence of the metastable radius $R_*$ and the height of the energy barrier for unknotting.
