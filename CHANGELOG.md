# Changelog

## 2026-07-20  -  Phase 8C pure-gauge code gate

- Added the neutral-link pure-gauge `3 x 3` toric-code preflight. The imposed
  reference passes the exact stabilizer, ground-space, distance, syndrome-gap
  and complete below-distance local-indistinguishability gates.
- Recorded the explicit boundary: this is an Abelian reference code with
  inserted neutral links/checks, not an ANTLER derivation, twist, fusion space
  or non-Abelian braid.
- Added T3a: a global `e<->m` duality calibration that maps all stars to
  plaquettes and rejects static-sign and bounded-basis-change false walls.
  The required finite local dislocation with endpoints remains explicitly open.
- Added T3b-reference: a finite non-CSS `[[7,1,3]]` triangle patch with
  X/Z-to-mixed string deformation through a central junction. It is a
  twist-to-boundary reference, not a separated-defect fusion or braid result.
- Added T4: odd-distance rotated planar code calibrations at `d=3,5` with four
  boundary twists. Both retain a GSD-2 locally protected code space with
  distance equal to the corner separation plus one; this is explicitly a
  boundary-twist code-space calibration, not physical fusion, mobile defects
  or a braid.
- Added T5a: a graph-based external reference with four degree-three interior
  twists. The `L=4,6` static parents have `k=3`, GSD `8`, a syndrome gap and
  distances `3,4`; a neighbouring graph changes five checks and has a
  noncommuting replacement, explicitly leaving code deformation and holonomy
  as the next gate.
- Added T5b: an exact dense `3x3` negative control for the direct linear
  check interpolation. Its eight-state band stays gapped, but is locally
  readable (`0.2937` at the midpoint), so the path is rejected as protected
  transport rather than being misreported as a braid precursor.
- Added T5c: a four-measurement stabilizer deformation that maps the two
  `3x3` interior-twist graph codes while retaining rank `6`, GSD `8` and no
  one-qubit logical action at any recorded stage. It remains a reference
  protocol pending outcome frames, logical transport and hardware derivation.
- Added T5d: complete logical-Pauli transport and all sixteen outcome-frame
  words for one T5c reference move. A second adjacent deformation and any
  commutator/braid conclusion remain explicitly open.
- Repaired the stabilizer-span equality key used by T5c--T5f: raw GF(2)
  row-echelon vectors were insertion-order dependent. All affected scripts
  were replayed; the target-check-only closing leg remains rejected and all
  earlier positive endpoint matches remain valid.
- Added T5g: two closed, outcome-conditioned external measurement loops now
  induce noncommuting GF(2) logical symplectic maps. The result needs three
  declared single-vertex measurements and is documented only as an audit-stack
  calibration, never as a physical non-Abelian braid.
- Added T5h: a bounded exhaustive anti-false-positive control. No protected
  return loop appears through six arbitrary one-vertex Pauli measurements, so
  the high-weight graph checks in T5f/T5g cannot be replaced by this short
  local-only grammar.
- Added T5i: every required high-weight graph check now has an exact one-ancilla
  parity-measurement circuit (4--6 CNOTs, residual below `1.45e-15`). This is
  an explicitly inserted external gate set, not an ANTLER or fault-tolerant
  measurement realization.

## 2026-07-20  -  Literature decision through July 2026

- Added a primary-literature map for non-Abelian twists, number-conserving
  Floquet ladders, adaptive gauging and recent D4/Fibonacci/S3 hardware
  demonstrations. It confirms the current `Z2` twist route as the nearest
  noncommutative target, while recording `D(S3)` plus fusion as a separate,
  long-term universal architecture requiring six-state non-Abelian links.
- Confirmed that the 2025 particle-conserving Floquet ladder is already the
  external Phase-8 benchmark; its paper does not reopen the registered ANTLER
  direct-realization obstruction without a new exact microscopic symmetry.

## 2026-07-20  -  External `Z2` proposal triage and Phase 8C-T0

- Audited the external `Z2`/Majorana/Higgs/branch-cut proposals. The neutral
  edge link is retained as a declared new reference resource; the claims that
  a three-arm star or a fixed `pi` sign alone provides non-Abelian fusion are
  rejected. The only non-Abelian research target retained is an `e <-> m`
  twist-defect route after a 2D `Z2` code has passed independent gates.
- Implemented and executed the separated 128D neutral-link star audit. At
  fixed `N=2` it has a 6D physical sector, exact Gauss and U(1) commutators,
  and zero projected bare hopping. A projected local density is non-scalar,
  preventing any false code/fusion/braid promotion. The next preflight is one
  square plaquette, not a T-junction dynamics run.
- Completed that one-plaquette T1 preflight: exact Gauss/U(1) remains, the
  magnetic `B_p` is gauge invariant and has a separated static flux sector,
  while a local density is still readable. The branch advances only to a 2D
  code-patch audit; no flux result is reclassified as an anyon or a tresse.

## 2026-07-20  -  Phase 8C handoff

- Sealed a next-session contract for a separated matter--link `Z2` reference
  model. The neutral link is a declared new resource and must not reuse the
  ANTLER charge-two mediator.
- Corrected the Gauss-law guardrail: for
  `c_v^dag tau^z_vw c_w+h.c.`, the generator is
  `G_v=(-1)^n_v product_(e incident v) tau^x_e`; a pure product of link flips
  is not a symmetry of the dressed hopping.
- Pre-registered a 4-vertex/3-edge, fixed-`N=2` 48D ED T0 audit before any
  fusion, T-junction or braid work. The current shared-matter pulse grammar
  remains closed; no RL/deep-RL search is justified inside it.

## Unreleased  -  Phase 8 canonical Floquet reproduction and finite-pulse bridge

- Corrected the Phase-8B shared-matter pulse interpretation. The registered
  integer virtual-Rabi closure has a small absolute error only because its SW
  signal is small: at `g/Delta=0.025` the relative target error is `1.0063`
  for all local `X/Y/Z` axes and the physical polar operation is nearly
  scalar. The same correction rejects the eight-segment C3 parity group as a
  gate (`relative error=1.00012` at the deepest point). Static/downfolded
  link and group-algebra results remain, while all pulse-gate wording is
  withdrawn. Future schedules must pass a relative-signal audit.
- Added a bounded exact duration search as a non-heuristic fallback before
  proposing a new pulse ansatz. Across `t_A,t_B in [0,8]` and
  `g/Delta=0.10,0.05,0.025`, no nonzero `XX` rotation meets the registered
  low-leakage/fidelity screen; the best rows are effectively identity. This
  rejects the direct two-segment pulse box, not off-resonant dressed or
  composite schedules in general.
- Closed a second, independent abrupt-control box before considering any
  learning-based pulse search. Exact repeated `A B` Floquet switching over
  20 rows (`g/Delta=0.05,0.025`; total duration `20,100`; segment duration
  `0.025..0.4`) yields no passing row: its best signal-relative error is
  `14.2956`, versus the pre-registered `0.1` threshold, and its leakage is
  never below `1e-4`. This excludes the stated abrupt window, not smooth,
  dressed or composite controls that have not yet been physically derived.
- Identified and measured the more fundamental shared-matter Floquet
  obstruction: averaging after Schrieffer--Wolff is not the same operation as
  downfolding the rapid microscopic average. For the registered sign echo,
  the static average carries `XX~(g/Delta)^2`, while the rapid microscopic
  average has exactly zero `XX` and a residual `IZ` that is `49.9542` times
  the static signal. Exact fast products converge to this microscopic average.
  A 30-row smooth `sin^2` echo preflight and a three-seed two-duration
  classical differential-evolution baseline also contain no local gate. This
  makes an RL search premature until a new, microscopically valid control
  algebra is derived.
- Tested the minimal balanced four-sign repair of that obstruction. Its static
  SW group is exactly selective and its microscopic first average is entirely
  scalar, removing the two-word `IZ` spectator. Nonetheless all 30 exact
  smooth-pulse rows return a nearly scalar logical operation (best
  signal-relative error `0.999239`). Sign balancing alone is therefore closed;
  a future route must derive a nonzero conditional term at the microscopic or
  explicitly controlled Magnus level before any optimizer is meaningful.
- Stress-tested a tempting conditional-JW bridge before promoting it. The
  local `XZ` gate is excellent only if a charge-two mediator is assigned an
  un-derived odd string weight. Under the frozen adjacent-rung convention
  (weight zero) and physical U(1)-charge counting (weight two), `c_XZ` is
  exactly zero. The result is therefore demoted to an algebraic control and
  identifies an explicit neutral/odd `Z2` gauge link as the missing resource;
  it is not a native ANTLER gate.

- Tested the first non-factorizing microscopic bridge after the disjoint
  reservoir obstruction. In a fixed-charge shared-matter block, a correlated
  Peierls-sign / pair-channel-phase echo cancels the isolated rail and walker
  flips in the Schrieffer--Wolff Hamiltonian while retaining a selective
  `X_rail tensor X_walker` coefficient. The deep fit is order `2.0000`, with
  low-frame capture increasing to `0.99914` at `g/Delta=0.025` and no
  serialized non-scalar companion. This is a static/downfolded local compiler
  result only: finite-pulse closure, full walker compilation, code integration
  and all twist/fusion/braid claims remain open.
- Ran a local virtual-Rabi population-return control for that bridge. One
  period per echoed segment at `g/Delta=0.025` gives leakage `5.04e-7` and
  polar distance to the piecewise SW target `3.95e-5`; single-cycle leakage
  scales as `(g/Delta)^4.894`. The registered one-block timing target passes
  over `[-0.2%, +0.3%]` and fails at `-0.3%` / `+0.4%`. A physical rung
  Peierls phase continuously rotates the derived conditional axis from `X` to
  `Y`; the independently pulsed `pi/2` control reproduces the same closure.
  The later relative-signal correction above reclassifies these as population
  closure controls, not gates; they are local deterministic data only.
- Completed the static one-link operator vocabulary rather than assuming the missing
  axis: a signed physical rail-potential bias, echoed with the pair-channel
  phase, derives `Z_rail tensor X_walker` with the same order-two deep scaling
  and a population-return pulse control (`5.41e-7` leakage, `7.91e-5` absolute
  polar distance).
  Together with real and `pi/2`-phase rung hopping, the registered block now
  compiles selective `X`, `Y`, and `Z` rail actions on its two-state walker.
  This remains a one-link library; no common multi-state walker or code is
  claimed.
- Ran the first explicit multi-link composition gate rather than extrapolating
  the one-link library. The 1488-state shared-walker C3 ring does generate
  `XXX`, but it is rejected because `XIX` is larger throughout the deep-SW
  regime: the fitted target/parasite powers are `7.05` and `4.56`, and a
  detuning scan `0.5..10` has no selective point (best parasite/target
  `3.93`). This closes the direct C3 grammar and prevents a false C4/code
  promotion; a derived cancellation or a different walker encoding is now
  required.

- Added a canonical TeNPy rung-MPO implementation with total `U(1)` and
  branch-parity `Z2` as exact tensor charges. It prevents the previously
  detected sector drift in unconstrained DMRG and reproduces the `L=8` ED
  sector energies below `2.7e-14`.
- Extended the registered fixed-density external point `U0=-2, alpha=0.5`
  through `L=16`. Its parity split falls from `9.10e-3` (`L=8`) to
  `2.01e-4` (`L=16`), and is stable from `chi=256` to `384`. The neutral gap
  falls with size, as expected for the reference's gapless total-charge mode;
  it is no longer misused as a topological rejection test.
- Added the periodic charged-sector curvature diagnostic. At the published
  `U0=-1.5, alpha=0.5, nu=1/3` point it is `0.381` (`L=6`) and `0.359`
  (`L=9`). This is small-system external evidence only, not an infinite-MPS
  topological-gap certificate.
- Added a fixed-density MPS scaling at that published point through `L=18`.
  Its long-size `chi=512` convergence run exceeded the execution limit before
  it wrote a result; it is explicitly retained as incomplete.
- Added an exact finite-pulse stroboscopic audit of the actual external
  Floquet sequence. At `eta=pi/2`, branch parity is conserved to numerical
  precision and the one-period effective-Hamiltonian error scales as
  `T^1.976`; a `0.1`-rad pulse-angle offset visibly breaks parity. This is a
  pulse-protocol benchmark, not an ANTLER-native derivation or a braid claim.
- Derived a local native ANTLER bridge to the Floquet starting Hamiltonian.
  In the rung-major convention, `theta=pi`, uniform leg hopping and two
  **separate** charge-two mediators per link produce the required attractive
  leg-density interaction at SW order without same-order cross-rail pair
  hopping. The local exact spectrum converges as `(g/Delta)^1.966`.
- Composed that bridge in an explicit `L=3,N=2` mediator block. An analytically
  registered virtual-Rabi closure reaches leakage `2.16e-6` and logical
  distance `7.33e-5`, with a symmetric strict timing window through `+/-0.3%`.
  Its coherent multi-cycle leakage grows approximately quadratically; Strang
  symmetrization improves the coherent unitary error but not the leakage. The
  next item is therefore composite/dressed leakage cancellation, not a phase
  or hardware promotion.
- Replaced the intermediate ideal rail rotation in the second Floquet segment
  by a direct coherent rotated-pair-mediator synthesis. Its local interaction
  factorization residual is `3.51e-16`, but the fixed-ratio leakage remains:
  this closes the hypothesis that the abrupt frame change was the primary
  source of the registered accumulation.
- Added exact direct-segment size controls through `L=6,N=3` and a
  fixed-logical-duration deep-SW convergence audit. At comparable duration,
  leakage and logical distance scale as `(g/Delta)^1.970` and
  `(g/Delta)^1.989` respectively over `0.05 -> 0.0125`; the accompanying
  `Delta=600 -> 9600`, `g=30 -> 120` demand is retained as an unqualified
  dynamic-hardware requirement.
- Added a direct-channel finite-ramp control at the deep-SW point. The imposed
  orthonormal channel path is numerically converged under ramp segmentation;
  both linear and `sin^2` transitions pass the registered local target through
  1% of a period and fail at 3%. A separate endpoint-angle audit brackets the
  coherent channel calibration requirement to `+/-0.003` rad under the same
  ideal finite-block conditions.
- Added microscopic nearest-link direct-channel crosstalk controls. The
  deterministic bracket passes through `0.3%` of `g` and fails at `1%` on the
  tested blocks. A seeded 100-realisation, fully correlated complex-error
  ensemble resolves the strict pass curve as `100,100,89,72,39` percent for
  RMS levels `0.1,0.2,0.3,0.5,1.0%`; this is retained as a coherent-control
  specification, not a material-noise qualification.
- Added an exhaustive microscopic-parity audit of the reused direct-channel
  mediator slots. `H0` and `H1` require incompatible fixed mediator charges,
  so the registered realization has no common microscopic rail `Z2`. This is
  an architecture-specific obstruction, not a no-go for separate species.
- Corrected the proposed segment-boundary error law. On the explicit
  Rabi-closed `L=3,N=2` cycle, the projected parity-flip metric scales as
  `(g/Delta)^5.931`, not `(g/Delta)^2`; a four-species fixed-charge control
  conserves parity at numerical precision. The Phase-8 T-junction theory
  document is revised accordingly and remains conditional/non-derived.
- Added the exact local Gauss-law algebra preflight for the next resource
  decision. A neutral gauge ancilla can forbid bare rail tunnelling exactly,
  but its single-vertex physical doublet remains readable by a local
  gauge-invariant dressed operator; it is retained as a gate, not a code.
- Independently stress-tested the displayed Phase 8B coarse-grained
  single-mediator Gauss gadget. The low-energy sector contains only even
  mediator-conversion orders, so its claimed order-three `X_L P_B X_R` term
  is exactly absent; this closes T3 as written while retaining T1/T2 and the
  possibility of a differently declared microscopic primitive.
- Audited the repaired repaired Lambda-walker loop. Its abstract internal layer
  produces only the intended coarse Gauss word at fourth order (deep fitted
  power `3.973`) with machine-zero non-scalar companions. The required
  density-conditioned neutral-walker hop is explicitly retained as a new,
  non-derived primitive; no physical Phase-8B code claim is promoted.
- Embedded that inserted primitive in a 96D fixed-charge two-rung Fock block.
  The block commutes exactly with its coarse Gauss generator; bare rail
  tunnelling has zero in-sector projection and the Gauss gap retains the
  fourth-order scaling. Inter-block pair transport remains the next gate.
- Completed that two-block transport gate in a 3584D fixed-charge block.
  Pair-only boundary transport commutes exactly with both coarse Gauss
  generators and is physical-sector preserving; the single-particle boundary
  hop is exactly sector-changing. This is not yet a paired-phase result.
- Propagated one pair through that two-block system exactly. Target population
  reaches `0.99755` at `g/Delta=0.025` with Gauss-sector leakage below
  `3.5e-13`; this validates an ideal single-pair shuttle only.
- Ran the necessary many-pair preflight before promoting that shuttle to a
  phase. In the projected periodic hard-core-pair chain, the mobile fluid has
  neutral and pair-addition gaps closing with fitted powers `-0.916` and
  `-1.043`. A strong-repulsion control has an isolated finite-size doublet,
  but a local density remains non-scalar on it (`0.4683` at `B=12`): this is a
  CDW control, not a protected code. The projection and interaction remain
  explicitly inserted rather than native.
- Replaced the inserted two-block pair-transfer term by an explicit
  positive-detuned charge-two mediator shared by the two local pair channels.
  The 29D exact bridge has SW coefficient error scaling as `(g/Delta)^1.982`.
  Its 3712D fixed-charge integration with both Lambda walkers transfers one
  pair with `0.99750` population at `g/Delta=0.025`, preserves both Gauss
  generators exactly and keeps leakage below `4.4e-11`. The four-leg pair
  channel and conditioned walkers remain declared, non-frozen resources.
- Added a neutral-walker route to a two-dimensional stabilizer parent. A
  1024D joint star/plaquette overlap yields only the commuting algebra
  `{I,A_s,B_p,A_sB_p}`: star and plaquette coefficients scale as fourth order,
  their product as eighth order. Under an explicitly conditional independent-
  gadget tiling, those coefficients retain the fourfold `3 x 3` toric ground
  space and a nonzero effective syndrome gap. This is an abstract Abelian
  stabilizer construction, not a native ANTLER phase or a non-Abelian braid.
- Closed the tempting coefficient-only braid loophole: every registered
  walker-derived parent is a polynomial in the same commuting stabilizers and
  acts as a scalar on its fourfold code. The next non-Abelian gate is now
  pre-registered as a microscopically derived e/m defect-deformation primitive,
  not an inserted braid or Hadamard matrix.
- Started that defect-resource gate with an exact five-step mixed-Pauli walker.
  It isolates the non-CSS pentagon `XZXZX` at fifth order; a joint 2560D
  pentagon/plaquette overlap retains only its commuting generated algebra.
  These are local building blocks motivated by twist-defect surface codes, not
  a complete branch cut, fusion space or braid.
- Formalized the programmable neutral-walker loop lemma: for a declared
  length-`ell` Pauli-conditioned cycle, the first non-scalar low-energy word
  is the full cycle at order `ell`; shorter closed paths are scalar. The
  measured four- and five-step controls validate the construction locally and
  turn the twist route into a specified geometry/crosstalk problem.
- Extended that local compiler to mixed checks without inserting a `Y` term:
  a `pi/2` phase on an X-conditioned walker hop is exactly `Y=iXZ`. The 64D
  single-loop audit isolates `YXZX` at fourth order, and the 1024D joint
  `YXZXII`/`XZIIXZ` audit closes on its commuting local algebra. The coherent
  conditional phase and full dislocation geometry remain declared resources;
  no fusion or braid claim is promoted.
- Completed the first finite mixed-check code gate. A six-check, seven-link
  algebraic completion seeded by the audited pair has rank six, one encoded
  qubit, distance three and exact syndrome gap `2J`. All 15 check pairs then
  passed explicit 2048D two-walker Schur audits. The C4 closed-walk selection
  rule further confines the formal simultaneous-walker SW series to the
  stabilizer group; its leading fourth/eighth-order syndrome-gap bound is
  positive. This remains an ideal walker grammar, not a native ANTLER phase,
  dislocation geometry, fusion space or braid.
- Closed the minimal native-bridge loophole before treating the new walker as
  derived. A fixed-charge rail-qubit plus disjoint pair-reservoir/charge-two
  mediator block factorizes exactly, including arbitrary tested complex-phase
  Floquet segments. It cannot create `W tensor X/Y/Z`; a non-separable
  shared-matter coupling is now explicitly required rather than assumed.

## Unreleased  -  Phase 7 audit harness

- Added an independently replayed Phase 7C quartet correction. Exact monomer
  projection identifies the crossed primitive as the two-monomial all-flip
  operator `F_p`, not full `XXXX`; it reproduces the archived
  `c_XXXX=-2.0304689728e-7`, establishes a squared alignment ceiling `1/8`,
  and separates the closure of the audited static grammar from any universal
  claim about all charge-two Hamiltonians.

- Added a theory-to-audit contract for parent-Hamiltonian candidates, including
  local projector algebra, exact symmetries, charge sectors, local probes,
  edge metrics, and boundary-anchor separation.
- Added a deliberately unprotected commuting-projector self-test. It has an
  exactly conserved edge operator but is rejected by local distinguishability,
  preventing a false positive from edge conservation alone.
- Independently transcribed and audited the supplied Phase 7 OBC fixed-point
  parent. Its doublet and flat gaps reproduce through `L=8`, but the local
  operator `X_j` acts logically at both edge and bulk. The candidate is
  therefore retained as a symmetry-restricted Ising/cat benchmark, not
  promoted as a localized-edge or topological code. The cell constraint is
  commuting PSD but not an idempotent projector.
- Calibrated the revised Ising strong-zero-mode recurrence against an equal
  bulk construction. The required alternating `(-r)^n` sign is recorded; the
  left residual follows the analytic support law while the bulk residual stays
  near `0.2`. This is a diagnostic control only.
- Added exhaustive weighted-charge-conserving local-algebra generation for
  Phase 7. It tests every local occupation matrix unit and Hermitian
  transition, including rail transfer, and rejects the Ising benchmark on a
  single rung as required.
- Added a separate Phase 7B exact 2D stabilizer reference control. The
  charge-frozen `3 x 3` torus has commuting star/plaquette projectors, GSD 4,
  distance 3 and no non-scalar Pauli action below that distance. It is
  expressly an imposed Abelian reference: no native ANTLER low-body mediator
  derivation, braid, or non-Abelian claim is promoted from this control.
- Added Phase 7C's constrained local microscopic environment for external
  solver proposals. Its fixed-charge four-edge block downfolds exactly to the
  monomer manifold and scores complete Pauli operator content, symmetries,
  capture and gap. The registered disconnected two-pair seed is a passing
  structural but failing target control: it produces no `XXXX` or `ZZZZ`
  stabilizer coefficient, so a local gap cannot be rewarded as a false hit.
- Added the Phase 7C replayable monitor (console, JSONL, CSV and Matplotlib)
  plus a complete optimizer-ready `S_0` payload for the registered `XXXX`
  negative control. The dashboard exposes loss/reward, monomer gap/capture,
  operator selectivity and parameter trajectories without creating a physical
  candidate claim.
- Hardened Phase 7C candidate-payload validation: a declared mediator-parity
  signature is now checked against the pair terms rather than ignored. The
  first mixed-species ring is archived as a rejection; correcting its
  false parity metadata leaves a connected but nonselective two-body model
  with a four-body target coefficient near `1e-10`.
- Added a classical discrete topology catalog (368 connected charge-two
  channel graphs) and exact coupling-order measurements. No catalog member
  clears the registered four-body signal threshold at fixed scale; two
  representative motifs scale as `XXXX~g^6` while their unwanted content
  scales as `g^2`. This corrects the external claim of a universally eighth
  order ring while preserving the tested-grammar selectivity rejection.
- Added a separate audit of a crossed charge-two mediator with plaquette
  support. It realizes a local `XXXX~g^4` component and calibrated static
  Z/ZZ cancellation, but the surviving four-body operator is not pure `XXXX`
  and remains far below the target amplitude. It is recorded as a qualified
  local direction, not a tiled ANTLER parent or topological-code result.
- Added an exact six-rung overlap stress-test for two crossed motifs. Their
  dominant all-flip interaction survives repeated local calibration, but the
  overlapping effective operators have spectral commutator norm `1`; this
  geometry is rejected as a commuting stabilizer parent.
- Closed the corresponding algebraic branch: exact logical-space tests show
  that identical all-flip terms commute only on disjoint or identical supports,
  and have spectral commutator norm `1` for every nontrivial overlap.
- Added a primary-literature repositioning before further search. The static
  crossed-mediator grammar is frozen; Floquet pair-hopping and a genuinely
  enlarged 2D code-gadget architecture are recorded as separate, unvalidated
  Phase 7D hypotheses with explicit preconditions before any RL run.
- Added the Phase 7D Floquet two-rung compiler control. With an explicitly
  declared sign-modulated interaction resource, it compiles stroboscopic
  number-conserving pair hopping exactly; this is not yet a derivation from
  frozen ANTLER hardware or a topological-phase result.
- Added a concrete positive-detuning charge-two bridge for that signed
  interaction. Complete Rabi excursions through same-rail and opposite-rail
  mediator channels provide inverse logical `ZZ` phases, and rail rotations
  compile their pair gate with residual `3.19e-16`. The bare rail parities are
  restored stroboscopically on the logical block; this is a time-multiplexed
  dynamic extension, not a new static symmetry.
- Added the four-rung even/odd multi-link compiler audit. Its noncommuting
  layers have a measured second-order Trotter residual while logical rail
  parities and unitarity remain at numerical zero. Crosstalk, pulse leakage
  and a many-body phase are deliberately not inferred from this control.
- Added the corresponding 472-state, four-rung explicit-mediator pulse audit.
  Closed pulses factorize to numerical precision at zero inactive coupling;
  under the registered coherent 1%-of-`g` inactive-channel error, the final
  monomer leakage is `5.67e-8`. This qualifies only that error model for the
  time-multiplexed control primitive, not a hardware-noise budget or phase.
- The same full-space audit rejects leaving intraleg hopping unrefocused during
  those pulses: it disrupts the opposite-rail Rabi closure, producing leakage
  `2.02e-4` already at `t_leg=0.1` and `1.98e-2` at `t_leg=1`. The next
  physical requirement is therefore a derivation of hopping freeze/refocusing
  or a calibrated full-Hamiltonian pulse, not a phase or RL scan.
- Added and exactly tested a staggered rail-plus-mediator phase echo. Its
  toggle identity leaves every nearest-link pair conversion invariant while
  changing the sign of intraleg hopping. The symmetric sequence gives
  controlled `n^-2.028` logical and `n^-4.058` leakage convergence. At
  `t_leg=1`, 1%-of-`g` crosstalk and 32 subcycles, leakage is `5.74e-8` and
  logical parity residual `4.29e-6`. This is an ideal finite-block control
  result pending an independent hardware derivation of the phase kicks.
- Added the conservative finite-duration potential-kick realization audit.
  Leaving the pulse Hamiltonian active during 16-subcycle kicks fails for
  `kappa=100..800`, with leakage remaining at `1.66e-2..4.59e-2`. The ideal
  echo is therefore not promoted to hardware: a switchable/fast or composite
  phase-kick mechanism remains a prerequisite.
- Added a continuous pair-preserving coherent-destruction-of-tunnelling
  alternative. A link-weighted mediator potential commutes exactly with every
  pair conversion while driving leg hopping to the first Bessel zero. The
  exact 472-state control reduces worst-case leakage from `1.98e-2` to
  `4.35e-4` at four cycles, but does not yet qualify a gate and requires a
  bandwidth/finite-frequency optimization before further promotion.
- Closed the corresponding narrow finite-frequency amplitude scan. The first
  Bessel zero remains the leakage optimum at four cycles across
  `xi=2.10..2.75`; nearby amplitude changes do not yield a hidden high-fidelity
  pulse and are not used to motivate a larger blind search.
- Added a direct Peierls-sign implementation of the ideal leg-hopping echo.
  In the explicit-mediator block it is algebraically identical to the ideal
  staggered echo and reaches leakage `5.74e-8` and parity residual `4.29e-6`
  at 32 subcycles with 1%-of-`g` crosstalk. It is retained as a conditional
  bridge because it requires 256 synchronous `0/pi` leg-phase switches per
  full schedule; no switching hardware derivation is claimed.
- Registered the first smooth finite-Peierls-ramp stress test as incomplete
  after it exceeded the local execution limit before final serialization.
- Completed the finite-Peierls-ramp block audit on the RTX 4070 Ti using
  exact `complex128` spectral propagation. With pair conversion active,
  ramp fractions `0.005` and `0.02` give leakage `1.06e-6` and `1.62e-5`.
  A 2/4/8-segment refinement now gives observed leakage convergence order
  `2.17` at both fractions, with final 8-segment leakage `8.82e-7` and
  `1.33e-5`. This promotes only the numerical ramp representation; a physical
  synchronization, bandwidth and error model remains required.
- Added a deterministic Peierls switching-error audit in the same exact 472D
  block. For registered `1e-4` local-control targets, a common `pi`-plateau
  phase offset passes through `5 deg` and first fails at `7.5 deg`; the
  signed-time imbalance passes through `5%` and first fails at `10%`. These
  are scan brackets only, not material-noise or hardware specifications.
- Completed 50-realisation-per-level quasi-static global Peierls-control noise
  audit. The strict local target passes `50/50` at `1%`, `35/50` at `3%`,
  `20/50` at `5%`, `8/50` at `10%` and `2/50` at `20%`; this remains a local
  control statistic rather than a hardware or code-noise claim.
- Added Phase 8 constant-density and sparse-Lanczos Floquet replication gates.
  The registered external candidate `U0=-2, alpha=0.5` is validated through
  `L=8,N=4`: `split/gap` falls from `0.1049` at `L=4,N=2` to `0.0301`, while
  the neutral gap remains `0.302`. Both registered `alpha=1` and `U0=0`
  controls are gapless. This is explicitly finite-size external-model evidence,
  not an ANTLER-native or topological-computing claim.
- Added a small-ED benchmark of the complete high-frequency
  number-conserving Floquet ladder of Defossez *et al.* rather than claiming
  that pair-hopping alone is sufficient. All registered `L=6` prolongations
  preserve the stroboscopic parity but fail the independently required local
  edge-localization gate; the sampled benchmark window is closed without
  refuting the published larger-size analysis.

## Unreleased  -  Phase 6L-6M protection controls

- Added Phase 6L-6M controls before any larger support run. The projected
  edge-commutator harness passes an exact Kitaev positive control
  (`epsilon_edge=0`) and rejects a bulk control (`2.0`). The external Iemini
  neutral gap falls to `0.2719927` at `L=10`; this does not prove an asymptotic
  closure, but provides no finite-size gap saturation.
- Corrected the Phase 6K interpretation: maximal support changed with size.
  At fixed `j=3`, the normalized residual is effectively flat (`6.3184` at
  `L=8`, `6.3057` at `L=10`), so the external finite-support generator remains
  unqualified as a protected edge mode or physical braid. The proposed `L=12`
  extension is not run as evidence for protection.
- Added the external-calibration status document and updated the status matrix
  and Phase 6 protection boundaries.

## Unreleased  -  Phase 5 pair-hopping preflight and composition branch audit

- Added an explicitly separate, number-conserving two-wire pair-hopping
  reference family and a finite-size local-indistinguishability preflight.
  The scanned minimal bridge is rejected: its low-splitting parity doublet is
  not edge-localized, so no braid calculation is promoted from it.
- Reproduced the published Iemini et al. number-conserving parent Hamiltonian
  on its exact λ=1 line through L=8. The edge-transfer locality trend is now a
  verified external benchmark, explicitly not an ANTLER derivation.
- Added finite-support Iemini braid audits through L=10, including a
  rung-major fermionic-order convention check, a nonzero-commutator gate and
  raw Yang--Baxter residuals. Finite-support leakage remains explicit, so the
  result is a calibration trend rather than an exact finite-size braid claim.
- Started Phase 6 native ANTLER discovery with a scalar three-leg family and a
  distinct chiral plaquette ring-exchange family. Both fail the mandatory
  local-indistinguishability preflight and are archived as no-go results.
- Added a Phase 6C local flux-interferometer derivation. At pi flux it cancels
  one-particle cotunnelling while a mediator interaction produces an isolated
  fourth-order pair transfer; it is recorded as a microscopic ingredient only,
  not as a tiled ANTLER code, topological phase, or braid.
- Added a charge-two mediator route with exact branch parities and a derived
  second-order pair transfer. Its explicit first ladder is also archived as a
  no-go scan: the best L=4 near-doublet fails the spectral threshold and does
  not improve at L=5.
- Added a protection-first Phase 6 protocol. The rank-three factorization of
  the external Iemini bond interaction now yields an explicit multiplet-
  mediator bridge whose large-detuning limit converges to that external parent.
  Finite-support edge-operator and charge-sector audits remain explicit: the
  former is not yet quasi-conserved at accessible sizes, and the latter is
  globally compressible without charge-sector control.
- Added a first-order Schrieffer--Wolff dressing audit for the mediator-bridge
  edge operator. It cancels the divergent virtual-mediator contribution but
  leaves a finite projected commutator, so it is retained as a no-go/diagnostic
  rather than promoted to a braid.
- Added a compiled matrix-free parent action and extended the external
  finite-support edge audit through L=10. The residual decreases at maximal
  available support, but has not reached a controlled asymptotic protection
  bound and remains explicitly non-promoted.
- Corrected the Phase 4.7 `n=8` composition analysis. The principal
  square-root used in the odd-in-theta extraction aliases at large accumulated
  phase; the summary now reconstructs the physical branch from cached raw
  `+/-theta` differentials and records the branch guard explicitly.

## Unreleased  -  Phase 4.7 closure and Phase 5 audit

- Added Phase 4.7 campaign scripts for deep-limit / timestep separation, path deformations and multi-cycle composition; their status is explicitly incomplete pending final controls.
- Added the digital-transfer lemma and a dedicated Phase 4.7 campaign-status document.
- Added N=3 pinned-mediator spectral, holonomy and consolidated no-go audits.
- Added a commutator gate to the N=3 braid diagnostic: a Yang--Baxter residual is not promoted when the two controls commute below threshold.
- Recorded the L=14 synthetic spinor code as locally distinguishable and therefore not a passive topological memory.
- Reclassified synthetic SU(2) link results as faithful transport of imposed Wilson links, not emergent non-Abelianity.
- Added effective-Majorana and phase-biased Kitaev T-junction benchmark audits with explicit non-ANTLER claim boundaries.
- Added Phase 5 literature positioning before any number-conserving microscopic extension.
- Refreshed the README read-first list, result-status matrix, archive manifest and checksums.

## v0.7  -  2026-07-18

- Added Gaussian path-deformation audit and explicit finite-depth path-dependence no-go.
- Added compact digital shuttle and sequential rung-transfer correction.
- Added distance-invariance tests for `R=3,4,5`.
- Added localization-depth series, power-law extrapolation and localization-weight correlation.
- Added exact Fock-path string-counting enumerator.
- Added static-disorder tests through `sigma/J2=0.20` with two seeds per level.
- Added analytic zero odd-dynamical-response argument at `theta=0`.
- Added exact isolated two-level handoff lemma with multiple ramp families.
- Added French executive summary, result-status matrix, latest scientific synthesis and incomplete-run register.
- Preserved partial `D=8`, multi-cycle and high-noise logs as incomplete rather than promoting them.
- Added digital `dt=0.125`, `theta=0.6`, and `theta=0.9` completed audits.
- Added explicit Abelian braid-algebra no-go and non-Abelian target comparison.
- Added Phase 4.6 accelerated local-Trotter prototype as an exploratory branch.

## v0.6  -  2026-07-18

- Added complete logical-subspace Strang propagation.
- Added exact orthonormal dressed logical frame.
- Added singular-value leakage audit and polar decomposition.
- Added odd-matrix extraction for the cleaned logical phase gate.
- Added theta=0.3, 0.6, 0.9 linearity runs.
- Added timestep and shuttle-length checks.
- Added recapture harness showing the earlier large leakage was primarily numerical.
- Consolidated all phases and original uploaded sources into one reproducible archive.

## v0.5

- Validated scalar dynamic shuttle exchange and differential odd phase.
- Added convergence and path-length studies.

## v0.1–v0.4

- Static ED model, edge-code construction, symmetry no-go, filament diagnosis, proximity lemma, cat encoding, and manual holonomy studies.
