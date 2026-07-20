# Phase 5 to 8 research archive

This archive records tested routes from the validated Abelian ANTLER phase
primitive toward a protected non-Abelian operation. Negative results are kept
because they determine the next architecture more reliably than an unqualified
positive scan.

## Status by phase

| Phase | Question | Public result | Boundary |
|---|---|---|---|
| 5 | Can the pinned three-particle mediator route generate a usable noncommuting exchange? | No. The tested operations commute within the registered numerical tolerance and the small internal split prevents an adiabatic qualification. | This rejects the tested static mediator route, not every three-particle construction. |
| 6 | Can the native ladder reproduce a protected number-conserving reference phase? | No native candidate in the registered scalar, density and three-leg scans passes the joint spectral and local-indistinguishability tests. | External reference models are calibration controls, not ANTLER realizations. |
| 7 | Can charge-two mediators produce a selective two-dimensional stabilizer parent? | No for the registered static link-local grammar. The desired four-body signal occurs at higher perturbative order than dominant two-body terms. | This is a grammar-specific no-go, not a no-go for every future ANTLER extension. |
| 8 | Can the existing microscopic resources derive the neutral link required for a protected non-Abelian route? | No. The direct Jordan-Wigner bridge gives zero physical conditional neutral-link coupling for the registered string weights. | A neutral $Z_2$ link is a missing resource, not an adjustable parameter. |

## Phase 5

The static $N=3$ mediator route was tested both algebraically and dynamically.
Its exchange candidates have a commutator below the pre-registered threshold,
and the braid relation is not interpreted when the commutator is null. The
near-crossing in the same route also prevents a credible adiabatic braid
claim. The relevant scripts and JSON files are under `experiments/phase5/`
and `results/phase5/`.

## Phase 6

External number-conserving reference models were reproduced only as audit
calibrations. The native ANTLER scans then tested local distinguishability,
charged and neutral gaps, and edge versus bulk response. The registered native
candidates do not meet all criteria simultaneously. The central documents are
`PHASE6_NATIVE_STATUS.md`, `PHASE6_EXTERNAL_CALIBRATION_STATUS.md` and the
matching scripts and results under `experiments/phase6/` and `results/phase6/`.

## Phase 7

The static charge-two route was tested by exact local downfolding, complete
Pauli decomposition, perturbative-order fits and overlapping-plaquette
audits. It fails for three independent reasons: insufficient four-body
selectivity, incomplete all-flip operator content, and noncommuting overlap
terms. A commuting fixed-point parent is retained only as an Ising/cat
diagnostic benchmark because local operators remain logically readable.

The public Phase 7 archive includes the exact scripts, JSON outputs and the
following key documents:

- `PHASE7_ADVERSARIAL_AUDIT.md`
- `PHASE7_FIXED_POINT_PARENT_AUDIT.md`
- `PHASE7C_QUARTET_CORRECTION_ADDENDUM.md`
- `PHASE7C_PERTURBATIVE_ORDER_AUDIT.md`

## Phase 8

The archive separates external reference controls from native derivations.
External gauge-code and measurement references verify that the audit stack can
recognize code-space, commutator and fault-channel conditions. They are not
presented as ANTLER hardware. The native microscopic bridge is closed by a
scoped no-go: the required odd neutral conditional link is absent from the
frozen charge-two and Jordan-Wigner resources.

The new architecture contract is explicit in `ANTLER_Z2_EXTENSION_CONTRACT.md`:
a neutral two-level link, a derived dressed hop, a two-dimensional protected
code and physical twist defects must be established in that order before any
non-Abelian claim.

## Publication boundary

This archive supports a second paper on constraints and resource requirements
for extending a number-conserving correlated-hopping ladder beyond its Abelian
phase primitive. It does not support a paper claiming a non-Abelian braid,
a protected fusion space, universal topological computation or an experimental
implementation.
