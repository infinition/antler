# ANTLER Digital Braiding Lemma

## Scope

This statement applies to the **sequential digital shuttle** in the strict
localized and adiabatic limit. It does not assert that the finite-depth ladder
is already a topological quantum computer.

## Local handoff

Consider an isolated two-state transfer between localized configurations
`|a>` and `|b>`:

\[
H(s)=\begin{pmatrix}
-D[1-q(s)] & -J e^{+i n\theta}\\
-J e^{-i n\theta} & -D q(s)
\end{pmatrix},\qquad s\in[0,1].
\]

Here `q(0)=0`, `q(1)=1`, and `n` is the Jordan-Wigner string occupation seen
by the directed hop. Under adiabatic parallel transport, the state initially
localized on `a` reaches `b` with the oriented link phase

\[
|a\rangle\longrightarrow e^{-i n\theta}|b\rangle
\]

up to a common dynamical/global phase. The result is independent of the ramp
shape `q(s)`. `run_phase4_5_two_level_lemma.py` verifies this for linear,
`sin^2`, and smoothstep ramps with a maximum phase error below `3e-14`.

## Path composition

For a sequence of isolated handoffs, phases multiply, so the odd statistical
phase is

\[
\phi_{\rm odd}=-\theta\sum_j \eta_j n_j,
\]

where `eta_j=+1` for a forward oriented crossing and `-1` for its reverse.

- Trivial round trip: every activated string crossing has its inverse, hence
  `sum eta_j n_j = 0`.
- Sequential exchange: one uncancelled activated crossing remains, hence
  `sum eta_j n_j = 1` and `phi_odd=-theta`.

## Finite-depth corrections

At finite trap depth the instantaneous state has tails on additional
configurations. The isolated two-state reduction is no longer exact, and the
measured phase becomes

\[
\phi_{\rm odd}=-\theta+\delta\phi(D,J,J_\perp,L),
\]

with corrections controlled by leakage from the localized transfer manifold.
The numerical data show:

- path-length spread for `R=3,4,5`: `3.54e-4` in slope;
- monotonic approach toward `-1` as the trap depth increases;
- exact dynamic odd response at `theta=0` equal to zero because
  `dH/dtheta|_0` is purely imaginary antisymmetric while the instantaneous
  eigenstates can be chosen real.

Therefore the remaining correction is geometric and caused by finite
localization, not by an uncancelled odd dynamical phase.

## Candidate asymptotic statement

The data and local lemma support, but do not yet constitute a full many-body
proof of,

\[
U_{\rm odd}\xrightarrow[D/J\to\infty,\;T\to\infty]{}
\exp(-i\theta Z/2)
\]

up to a global phase, with leakage tending to zero when the adiabatic time is
scaled with the trap depth.
