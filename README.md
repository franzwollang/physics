# README

## Overview

This repository collects my speculative physics papers along with supporting files. These are **drafts**: incomplete, likely incorrect but earnest attempts to wrestle with long-standing questions I have carried for years. The drafts are also posted to Zenodo and clearly labeled as drafts there.

The goal is not polished correctness but translation: turning a messy, intuitive web of ideas into a form that can at least be communicated, critiqued, or dismissed on their own terms.

## Why This Exists

I studied physics at the university level for about five years without earning a degree. Many of my friends continued on to obtain Masters and PhDs, which means I have lived in a community of scientific discussion while remaining on the margins of formal academia. Over the years, I have tried and failed to fully communicate my own half-formed ideas even to people with such backgrounds.

This repository is an attempt to break that logjam. With LaTeX, LLM assistance, and stubborn persistence, I’m writing these drafts into existence. I’m not claiming to have “the answers”. I’m claiming that these questions keep resurfacing for me — and I want to give them and some kind of roughly coherent response form, whether it ends up being completely accurate or not.

## Goals

1. Provide a coherent record of a speculative framework that can be evaluated on its content.
2. Connect the papers through the personal/philosophical motivations that drive them.
3. Document the evolution of my thought, mistakes included.
4. Make critique possible by putting everything in the open.

## What Counts as Success (Rigor vs Coherence)

My target is not full formal rigor everywhere. I aim for a mostly consistent, coherent narrative where the pieces work together to form a compelling gestalt.

I will keep improving accuracy, cross‑checks, and consistency, and I do intend to confront data as time allows. Publication is not the goal although it would be a welcome cherry on top. Without an institutional affiliation, I can take more risks and I intend to use that freedom to the fullest as it's one of the few advantages I have here.

## How to Read These Drafts

- Treat them as **drafts** and thought experiments.
- Do not assume correctness but do assume sincerity.
- Look for predictions or anomalies as handles for critique.
- Appreciate the elegance of some of the concepts and approaches in this framework, even if you take them with a huge grain of salt.

## Repository Contents

- `/papers/` — LaTeX sources for individual papers.
- `/figures/` — Images, diagrams, and supporting visuals.
- `/notes/` — Background notes, scratchwork, and intermediate ideas.

## Current Drafts

- Spacetime from First Principles: Free-Energy Foundations on the Infinite-Clique Graph (Draft) — DOI [10.5281/zenodo.15843979](https://zenodo.org/records/17000342)
- A Sea of Noise: Relativity from a Thermodynamic Force in Scale-Space (Draft) — DOI [10.5281/zenodo.17000665](https://zenodo.org/records/17000666)

## Transparency and Limitations

- Drafts are LLM-assisted, with a high likelihood of mathematical errors.
- The point is not perfection but communication.
- I acknowledge that this project is in many ways indistinguishable from crank/crackpot physics. My way of countering that is total transparency: show the drafts, show the reasoning and motivations, include as much maths as I can, include falsifiable predictions when possible, leave the work open for critique, and make meaningful attempts to address flaws that are pointed out.
- My background is informal but informed: years of study without a degree, embedded in a community of physicists.

## How I Use AI (Workflow and Models)

- I lean into LLMs' strengths: fast iteration, broad fluency across a huge swath of fields, and the ability to connect proposals to known patterns or suggest new ones for me to interrogate. Rapid cycles help discard clearly bad ideas and converge on not‑immediately‑terrible ones worth deeper work.
- I use multiple state‑of‑the‑art models. This project only became viable in the last ~12–18 months as models crossed qualitative thresholds: more abstract reasoning and better coherence across more complex contexts. They still require constant guidance and a tight leash; I remain the editor and final arbiter.
- Even with the heavy use of LLMs, these drafts represent months of work. That is already only a small fraction of the time it would have taken otherwise.

---

## Motivations (Problems I Keep Returning To)

- **Gravitational Potential.** Where is gravitational potential “stored”? In GR it isn’t stored—gravity is geometry—but in Newtonian reductions we talk about a scalar Φ(x) and global “potential energy” bookkeeping. For electromagnetism and other short-range forces, potential energy makes sense as a bookkeeping mechanism of local underlying dynamics. But for gravity it feels reified into an ontic entity. Edge cases make the discomfort concrete: naive pairwise sums over infinite or unbounded matter are ill‑defined and need background subtraction (the “Jeans swindle”), and for large‑N systems the total can hinge on boundary prescriptions. So even though the two‑body potential is perfectly finite for r>0, the global bookkeeping—what the energy is, where it “lives,” and how it “moves”—remains philosophically unsatisfying.

- **Indistinguishability in Quantum Mechanics.** Why should particles be considered strictly indistinguishable, aside from quantum numbers? This postulate/derivation (depending on framework) predisposes interpretation toward probability and indeterminacy. If the apparent randomness of QM is to be reconciled with determinism (which I see as the only consistent philosophical stance), then randomness must emerge from chaos in hidden internal states—while acknowledging Bell/Kochen–Specker require any such completion to be nonlocal/contextual. Particles should have some kind of extent and internal, time-dependent configuration — not be idealized pointlike excitations with only a probability amplitude over configurations.

- **Before the Big Bang.** In standard ΛCDM/GR, extrapolating the classical solution places the origin of our cosmic time at the Big Bang; “before” lies outside the model’s domain. But there being a "start" to time creates a huge philosophical problem. And even if expansion can be pictured conformally (a unit sphere expanding into itself), the philosophical problems remain unless something like conformal self-similarity applies to everything. Hours-long Mandelbrot set zoom videos on Youtube and multiple LSD trips gave me the hunch that both views — eternal universe and Big-Bang-like phenomenology — might be reconciled. A single fractal substance, endlessly self-dividing, gives us an eternal universe with no dualism problem, and still allows local patches to behave like our observed big bang. But how?

- **Bell’s Trilemma.** Realism, locality, and determinism are too simple, practical, and pure to give up. Yet Bell’s theorems show that, under standard assumptions (locality, measurement independence, outcome definiteness), you cannot keep **all three simultaneously**. In years of reading philosophy and arguing determinism online, Bell’s results were always the final refuge of inherent randomness (Copenhagen) or free-will proponents, the final wall that I was never able to fully undermine since it relied on an open problem in physics. I’ve always believed there _must_ be a way through that keeps all three by rethinking assumptions.

- **Dimensionality of Reality.** Why is our universe three-dimensional? Number bases are arbitrary but 3D seems to be uniquely stable in the physics we perceive. This does not appear to be a biological artifact. It seems to be embedded in the fundamental fabric of our world. What makes 3D special?

- **Evolution as Algorithm.** Evolution is substrate- and implementation-neutral. It suggests—at least as a motivating analogy—that similar selection dynamics could apply in physics itself. Freak waves form from random fluctuations in fluids — 200-foot walls of water from chance. Why couldn’t unstable or semi-stable structures form in a pre-law physics substrate, then alter their environment and recursively stabilize, just as cyanobacteria oxygenated Earth and reshaped the fitness landscape? If laws (or stable environments that give rise to the perception of the laws) evolve, what would their selection rules look like?

---

## Inspirations

These are not proofs but **hunches and narrative outlines** of how solutions to the motivating problems might hang together in a coherent picture.

- **Particles as Solitons.** The idea is that indistinguishability can be reframed if all particles are some kind of standing waveform, effectively solitons stabilized by a quartic potential. This picture gives particles extent and internal structure. But if you go this route, you then have to deal with a huge cascade of questions: how do interference, non-commutativity, and operator algebra still emerge if particles are solitons? Can the known behavior of quantum fields really be rederived from this basis? This is a monumental risk but it feels like the right direction.

- **Gravitational Potential Recast.** Instead of asking where gravitational potential is stored, imagine that the kinetic energy was always there within infalling matter just “latent”, and only becomes expressed depending on curvature. If the standing waveform has a stable attractor in phase-space and curvature shifts where that attractor lies, moving through curvature would change how the total energy is expressed. This resonates with what we already accept for photons: redshift and blueshift as curvature changes their expression of energy. The leap is to extend that to all particles. In this picture, contraction of a waveform means some hidden portion of its energy comes out as kinetic. But keeping this consistent with known physics pushes strongly toward a conformal framework.

- **Fractal Universe.** If the gravitational problem demands conformal covariance and the big bang/eternal universe reconciliation problem requires a fractal-type evolutionary symmetry, then we are naturally led toward a fractal universe picture: zoom into any point, or out from any point, and the structure persists. Like an endless Mandelbrot zoom but with physical content. This might reconcile an eternal fractal substrate with the local phenomenology of a Big Bang. The details are tricky — how to make such a framework yield predictions rather than just metaphors — but the intuition is strong.

- **Thermodynamic Free-Energy Functional.** A natural starting point for building an evolutionary picture of physics is to posit a free-energy functional over some graph-like primitive. Fixed points of this functional correspond to interesting or observed outcomes. Such a picture should be single-substance (to avoid dualism issues), compatible with infinite self-subdivision and recombination in a hyperbolic way, and capable of producing phenomenology from an undifferentiated substrate once dynamics are applied. Energy minimization and entropy maximization are so core to physics that using them as the basis of the extremized operator feels like the right minimal foundation.

- **Bell’s Theorems and Locality.** Perhaps realism, determinism, and locality can all survive if locality itself is redefined. Locality may be emergent, not fundamental — a property that arises from a higher-dimensional substrate (this idea is not novel). The challenge then is twofold: recover the 3D lattice-like behavior we observe and simultaneously accommodate higher-dimensional behavior that looks non-local but remains compatible with relativity.

- **3D as Fixed Point.** If three dimensions are some kind of attractor or fixed point, what is the primordial substrate? Both the dimensionality question and the Bell trilemma motivate seeking a minimal starting point that can generate the needed degrees of freedom. One candidate is an infinite-clique graph over an uncountable set, which could in principle reduce or build up to the stable 3D structure we observe.

---

## Resolutions

The following drafts turn these hunches into more concrete frameworks, addressing key aspects:

### Spacetime from First Principles: Free-Energy Foundations on the Infinite-Clique Graph

- **Single substance** → Concrete: An infinite clique of real‑valued nodes with a volume‑normalized free‑energy functional avoids dualism while generating rich emergent structure.
- **3D as fixed point** → Concrete: Spectral dimension emerges as a statistical fixed point at d\*=3 from minimizing an information‑theoretic cost functional, not assumed a priori.
- **Bell's theorems and locality** → Concrete: Entanglement as phase‑link hysteresis breaks locality in a higher‑dimensional sense (phase bundle) while preserving relativistic causality and no‑signalling, allowing realism + determinism.
- **Evolutionary Picture of Physics** → Concrete: The extremization principle over graph links produces vacuum self‑organization and emergent locality, and measurable correlations with predictable lifetimes.

### A Sea of Noise: Relativity from a Thermodynamic Force in Scale-Space

- **Particles as solitons** → Concrete: Matter emerges as localized amplitude excitations of a single complex field; radiation as massless phase waves. Internal structure without losing universal kinematics.
- **Gravitational potential recast** → Concrete: Scale‑space thermodynamics makes solitons shrink in higher noise regions, creating an attractive force toward \(\nabla\tau\). This replaces global pairwise bookkeeping with a local free‑energy density: define an effective potential Φ\*eff(x) ∝ E_eq(τ(x)), so F = −∇Φ_eff(x).
- **Conformal framework** → Concrete: Local units co‑scale with the field's equilibrium parameters, making dimensionless physics approximately invariant under observational window changes.

### Scale-Downshifting as a Unified Origin for Dark Matter and Dark Energy

- **Fractal universe** → Concrete: A static, scale‑free vacuum spectrum \(S(k)\) accessed through a sliding observational window reconciles an eternal background with local Big Bang phenomenology. Dark matter and dark energy emerge as spatial vs temporal responses to the same noise field, avoiding Occam's razor violations.

### Horizon Duality (Black Holes)

- **Before the Big Bang** → Concrete: Black hole interiors as nascent universes with their own CMB, structure formation, and thermal history. Horizon duality makes cosmic expansion equivalent to universal matter downscaling, reconciling eternal fractal substrate with local Big Bang phenomenology.
- **Fractal universe (refined)** → Concrete: Infinite zoom hierarchy alternating 2D shells and 3D bulks with average dimension \(D\_{\rm eff} = e\), but local patches always integer‑dimensional. Each level has finite energy/computation via the scale window principle.
- **The Start of Time Paradox** → Concrete: Thermodynamic scale‑space force \(\dot{\sigma} < 0\) provides a universal arrow of time at every hierarchical level without there ever being any beginning or end. Gravity becomes the manifestation of this scale‑space reduction, making time's direction emergent from infinite fractal conformal downshifting.

---

## Open Questions / Unresolved Points

Here are ongoing puzzles that drafts have not yet satisfactorily (to my loose standards) touched upon.

- How to make soliton-based particles consistent with known quantum interference and algebraic structures?
- How to recover higher symmetries/gauge forces in a clean, organic way from within the framework?
- Why does the gauge rotor ladder stop at SU(3)? A fixed-point or thresholding argument similar that used for explaining "why 3D" would be ideal.
- How to derive the fractal vacuum's S(k) power spectrum from first principles?
- How to connect the alternating 2D/3D hierarchy to cosmological observables?

---

## Contributing (Issues and Moderation)

- If you find the direction interesting and want to contribute, please open an issue. Constructive criticism is welcome; I’ll try to make amendments or clarifications in good faith.
- Non‑constructive criticism or schizoposting/numerology/sacred‑geometry type takes will be closed immediately; repeat offenders may be blocked/blacklisted.
- This project is squarely within the realm of fringe theories, so I do not wish to gatekeep or fall prey to No True Scotsman fallacies when it comes to unorthodox ideas. However, please display some minimum understanding of standard modern theories and some minimum of respect for them; the goal of science is to _build_ on what has come before, not completely overturn it or find "gotcha"s. Almost all scientific paradigm shifts throughout history have been a result of insightful reinterpretations and theory-crafting on known foundations and data.

---

## License

Released under a permissive license (e.g. CC-BY or MIT, TBD). You may reuse, critique, and build upon it with attribution.
