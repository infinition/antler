"""Constrained local environment for Phase 7C microscopic discovery.

The environment works on a four-edge patch, the support of one target toric
star or plaquette.  It never substitutes a global 3x3 microscopic ED: that
fixed-charge Hilbert space is already billions of states before mediators.
Instead it exactly diagonalizes a small charge-two-mediator block, downfolds
its isolated monomer manifold, and measures its complete four-qubit Pauli
content against a requested four-body target.

The grammar is deliberately restrictive.  Every generated term conserves
weighted U(1); branch-parity labels are either the bare rail parities or their
mediator-dressed generalization.  Passing this local environment is a search
filter, not a proof of a global ANTLER phase.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import comb
from typing import Iterable

import numpy as np


EDGE_COUNT = 4
LOW_MODE_COUNT = 2 * EDGE_COUNT
LOW_MANIFOLD_DIMENSION = 1 << EDGE_COUNT
PAULI_MATRICES = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    "Y": np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
    "Z": np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
}


def rail_mode(edge: int, rail: int) -> int:
    if not 0 <= edge < EDGE_COUNT or rail not in (0, 1):
        raise ValueError("edge must be 0..3 and rail must be 0 (a) or 1 (b)")
    return 2 * edge + rail


def _pair_signature(first: int, second: int) -> tuple[int, int]:
    """Parity carried by a charge-two mediator converting this low pair."""
    return (
        int(((first & 1) == 0) + ((second & 1) == 0)) % 2,
        int(((first & 1) == 1) + ((second & 1) == 1)) % 2,
    )


def _mode_label(mode: int) -> str:
    return f"{'a' if mode % 2 == 0 else 'b'}{mode // 2}"


@dataclass(frozen=True)
class PairConversionChannel:
    """One hard-core charge-two mediator and its allowed pair conversions."""

    name: str
    detuning: float
    coupling: float
    phase: float
    pair_terms: tuple[tuple[int, int, complex], ...]

    def parity_signature(self) -> tuple[int, int]:
        if not self.pair_terms:
            raise ValueError(f"channel {self.name} has no pair conversion")
        signatures = {_pair_signature(first, second) for first, second, _ in self.pair_terms}
        if len(signatures) != 1:
            raise ValueError(
                f"channel {self.name} mixes pair branch parities; it would not carry a definite mediator parity"
            )
        return next(iter(signatures))

    def validate(self, maximum_sw_ratio: float) -> None:
        if not self.name:
            raise ValueError("mediator channel needs a name")
        if self.detuning <= 0.0:
            raise ValueError(f"channel {self.name} must have positive detuning")
        if abs(self.coupling) / self.detuning > maximum_sw_ratio:
            raise ValueError(f"channel {self.name} leaves the registered Schrieffer-Wolff regime")
        seen: set[tuple[int, int]] = set()
        for first, second, coefficient in self.pair_terms:
            if not (0 <= first < LOW_MODE_COUNT and 0 <= second < LOW_MODE_COUNT and first != second):
                raise ValueError(f"channel {self.name} references an invalid low-mode pair")
            pair = tuple(sorted((first, second)))
            if pair in seen:
                raise ValueError(f"channel {self.name} repeats pair {pair}")
            seen.add(pair)
            if abs(complex(coefficient)) <= 1e-15:
                raise ValueError(f"channel {self.name} contains a zero conversion coefficient")
        self.parity_signature()


@dataclass(frozen=True)
class RailHop:
    """A same-rail number-conserving microscopic hopping on the local patch."""

    left_edge: int
    right_edge: int
    rail: int
    amplitude: float
    phase: float = 0.0

    def validate(self, mott_u: float, maximum_sw_ratio: float) -> None:
        rail_mode(self.left_edge, self.rail)
        rail_mode(self.right_edge, self.rail)
        if self.left_edge == self.right_edge:
            raise ValueError("a rail hop must connect distinct edges")
        if abs(self.amplitude) / mott_u > maximum_sw_ratio:
            raise ValueError("rail hop leaves the registered Mott perturbative regime")


@dataclass(frozen=True)
class ZZCoupling:
    """A local density interaction, written as Z_i Z_j in the monomer sector."""

    left_edge: int
    right_edge: int
    strength: float

    def validate(self, mott_u: float, maximum_sw_ratio: float) -> None:
        rail_mode(self.left_edge, 0)
        rail_mode(self.right_edge, 0)
        if self.left_edge == self.right_edge:
            raise ValueError("a ZZ coupling must connect distinct edges")
        if abs(self.strength) / mott_u > maximum_sw_ratio:
            raise ValueError("ZZ coupling leaves the registered Mott perturbative regime")


@dataclass(frozen=True)
class MicroscopicCandidate2D:
    """A four-edge ANTLER-compatible local candidate in a constrained grammar.

    The grammar deliberately excludes direct four-rung stabilizers.  They may
    only appear in the low-energy effective Hamiltonian through virtual
    charge-two conversion and/or Mott-violating hopping processes.
    """

    mott_u: float
    channels: tuple[PairConversionChannel, ...]
    rail_hops: tuple[RailHop, ...] = ()
    zz_couplings: tuple[ZZCoupling, ...] = ()
    rail_biases: tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)
    maximum_sw_ratio: float = 0.15

    def validate(self) -> None:
        if self.mott_u <= 0.0:
            raise ValueError("Mott penalty must be positive")
        if not 0.0 < self.maximum_sw_ratio <= 0.25:
            raise ValueError("registered Schrieffer-Wolff ratio must lie in (0, 0.25]")
        if not self.channels:
            raise ValueError("candidate needs at least one explicit charge-two mediator")
        if len({channel.name for channel in self.channels}) != len(self.channels):
            raise ValueError("mediator channel names must be unique")
        for channel in self.channels:
            channel.validate(self.maximum_sw_ratio)
        for hop in self.rail_hops:
            hop.validate(self.mott_u, self.maximum_sw_ratio)
        for coupling in self.zz_couplings:
            coupling.validate(self.mott_u, self.maximum_sw_ratio)
        if len(self.rail_biases) != EDGE_COUNT:
            raise ValueError("one rail bias is required per local edge")
        if any(abs(bias) / self.mott_u > self.maximum_sw_ratio for bias in self.rail_biases):
            raise ValueError("rail bias leaves the registered Mott perturbative regime")

    def parity_type(self) -> str:
        return "bare_rail_parities" if all(channel.parity_signature() == (0, 0) for channel in self.channels) else "mediator_dressed_parities"


def seeded_perturbative_candidate() -> MicroscopicCandidate2D:
    """Return the nonzero but intentionally insufficient Phase 7C seed.

    Its two disjoint same-branch conversion channels induce pair processes at
    second order.  It cannot manufacture a connected four-edge stabilizer by
    itself, making it a useful negative reward control rather than a guessed
    solution.
    """
    return MicroscopicCandidate2D(
        mott_u=20.0,
        maximum_sw_ratio=0.15,
        channels=(
            PairConversionChannel(
                name="even_pair_01", detuning=10.0, coupling=0.5, phase=0.0,
                pair_terms=((rail_mode(0, 0), rail_mode(1, 0), 1.0), (rail_mode(0, 1), rail_mode(1, 1), 1.0)),
            ),
            PairConversionChannel(
                name="even_pair_23", detuning=10.0, coupling=0.5, phase=0.0,
                pair_terms=((rail_mode(2, 0), rail_mode(3, 0), 1.0), (rail_mode(2, 1), rail_mode(3, 1), 1.0)),
            ),
        ),
    )


def _weighted_basis(channel_count: int, total_charge: int = EDGE_COUNT) -> tuple[np.ndarray, dict[int, int]]:
    mode_count = LOW_MODE_COUNT + channel_count
    states = np.asarray([
        state for state in range(1 << mode_count)
        if (state & ((1 << LOW_MODE_COUNT) - 1)).bit_count()
        + 2 * (state >> LOW_MODE_COUNT).bit_count() == total_charge
    ], dtype=np.int64)
    return states, {int(state): position for position, state in enumerate(states)}


def _fermionic_sign(state: int, mode: int) -> float:
    return -1.0 if (state & ((1 << mode) - 1)).bit_count() & 1 else 1.0


def _annihilate_low(state: int, mode: int) -> tuple[int, float] | None:
    if not ((state >> mode) & 1):
        return None
    return state ^ (1 << mode), _fermionic_sign(state, mode)


def _create_low(state: int, mode: int) -> tuple[int, float] | None:
    if (state >> mode) & 1:
        return None
    return state | (1 << mode), _fermionic_sign(state, mode)


def _annihilate_pair(state: int, first: int, second: int) -> tuple[int, float] | None:
    first, second = sorted((first, second))
    item = _annihilate_low(state, first)
    if item is None:
        return None
    interim, amplitude = item
    item = _annihilate_low(interim, second)
    if item is None:
        return None
    final, sign = item
    return final, amplitude * sign


def _apply_hop(state: int, source: int, destination: int) -> tuple[int, float] | None:
    item = _annihilate_low(state, source)
    if item is None:
        return None
    interim, amplitude = item
    item = _create_low(interim, destination)
    if item is None:
        return None
    final, sign = item
    return final, amplitude * sign


def monomer_embedding(states: np.ndarray, index: dict[int, int]) -> np.ndarray:
    """Columns are the ordered one-particle-per-edge rail computational basis."""
    frame = np.zeros((len(states), LOW_MANIFOLD_DIMENSION), dtype=complex)
    for column, rails in enumerate(product((0, 1), repeat=EDGE_COUNT)):
        state = sum(1 << rail_mode(edge, rail) for edge, rail in enumerate(rails))
        frame[index[state], column] = 1.0
    return frame


def build_local_candidate_hamiltonian(candidate: MicroscopicCandidate2D) -> tuple[np.ndarray, np.ndarray, dict[int, int], np.ndarray, np.ndarray]:
    """Build the full fixed-charge local block and its monomer embedding."""
    candidate.validate()
    states, index = _weighted_basis(len(candidate.channels))
    H = np.zeros((len(states), len(states)), dtype=complex)
    mott_cost = np.zeros(len(states), dtype=float)
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        for edge in range(EDGE_COUNT):
            n_a, n_b = (state >> rail_mode(edge, 0)) & 1, (state >> rail_mode(edge, 1)) & 1
            cell_cost = float((n_a + n_b - 1) ** 2)
            mott_cost[column] += cell_cost
            H[column, column] += candidate.mott_u * cell_cost
            H[column, column] += candidate.rail_biases[edge] * (n_a - n_b)
        for channel_index, channel in enumerate(candidate.channels):
            mediator = LOW_MODE_COUNT + channel_index
            H[column, column] += channel.detuning * ((state >> mediator) & 1)
        for coupling in candidate.zz_couplings:
            z_left = ((state >> rail_mode(coupling.left_edge, 0)) & 1) - ((state >> rail_mode(coupling.left_edge, 1)) & 1)
            z_right = ((state >> rail_mode(coupling.right_edge, 0)) & 1) - ((state >> rail_mode(coupling.right_edge, 1)) & 1)
            H[column, column] += coupling.strength * z_left * z_right
        for hop in candidate.rail_hops:
            source, destination = rail_mode(hop.left_edge, hop.rail), rail_mode(hop.right_edge, hop.rail)
            item = _apply_hop(state, source, destination)
            if item is not None:
                new_state, sign = item
                amplitude = -hop.amplitude * np.exp(1j * hop.phase) * sign
                H[index[new_state], column] += amplitude
                H[column, index[new_state]] += amplitude.conjugate()
        for channel_index, channel in enumerate(candidate.channels):
            mediator = LOW_MODE_COUNT + channel_index
            if (state >> mediator) & 1:
                continue
            for first, second, coefficient in channel.pair_terms:
                item = _annihilate_pair(state, first, second)
                if item is None:
                    continue
                low_state, sign = item
                new_state = low_state | (1 << mediator)
                amplitude = -channel.coupling * np.exp(1j * channel.phase) * complex(coefficient) * sign
                H[index[new_state], column] += amplitude
                H[column, index[new_state]] += amplitude.conjugate()
    if not np.allclose(H, H.conj().T, atol=1e-11):
        raise RuntimeError("local candidate compiler produced a non-Hermitian Hamiltonian")
    return H, states, index, monomer_embedding(states, index), mott_cost


def _parity_labels(candidate: MicroscopicCandidate2D, states: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    first, second = [], []
    signatures = [channel.parity_signature() for channel in candidate.channels]
    for raw_state in states:
        state = int(raw_state)
        n_a = sum((state >> rail_mode(edge, 0)) & 1 for edge in range(EDGE_COUNT))
        n_b = sum((state >> rail_mode(edge, 1)) & 1 for edge in range(EDGE_COUNT))
        for channel_index, signature in enumerate(signatures):
            occupied = (state >> (LOW_MODE_COUNT + channel_index)) & 1
            n_a += signature[0] * occupied
            n_b += signature[1] * occupied
        first.append(-1.0 if n_a & 1 else 1.0)
        second.append(-1.0 if n_b & 1 else 1.0)
    return np.asarray(first), np.asarray(second)


def pauli_word(label: str) -> np.ndarray:
    if len(label) != EDGE_COUNT or any(letter not in PAULI_MATRICES for letter in label):
        raise ValueError("target Pauli label must contain exactly four I/X/Y/Z characters")
    operator = np.asarray([[1.0]], dtype=complex)
    for letter in label:
        operator = np.kron(operator, PAULI_MATRICES[letter])
    return operator


def _pauli_coefficients(matrix: np.ndarray) -> dict[str, float]:
    dimension = matrix.shape[0]
    rows: dict[str, float] = {}
    for letters in product("IXYZ", repeat=EDGE_COUNT):
        label = "".join(letters)
        coefficient = np.trace(pauli_word(label) @ matrix) / dimension
        if abs(coefficient.imag) > 1e-9:
            raise RuntimeError(f"non-Hermitian Pauli coefficient leaked into {label}")
        rows[label] = float(coefficient.real)
    return rows


def candidate_payload(candidate: MicroscopicCandidate2D) -> dict:
    """Compact, JSON-safe candidate description for an external optimizer."""
    candidate.validate()
    return {
        "local_patch": "four charge-frozen rung qubits supporting one candidate star or plaquette",
        "mott_u": candidate.mott_u,
        "maximum_sw_ratio": candidate.maximum_sw_ratio,
        "parity_type": candidate.parity_type(),
        "channels": [
            {
                "name": channel.name,
                "detuning": channel.detuning,
                "coupling": channel.coupling,
                "coupling_over_detuning": channel.coupling / channel.detuning,
                "phase": channel.phase,
                "mediator_parity_signature": list(channel.parity_signature()),
                "pair_terms": [
                    [_mode_label(first), _mode_label(second), [complex(coefficient).real, complex(coefficient).imag]]
                    for first, second, coefficient in channel.pair_terms
                ],
            }
            for channel in candidate.channels
        ],
        "rail_hops": [hop.__dict__ for hop in candidate.rail_hops],
        "zz_couplings": [coupling.__dict__ for coupling in candidate.zz_couplings],
        "rail_biases": list(candidate.rail_biases),
    }


def _parse_mode(value: int | str) -> int:
    if isinstance(value, int):
        if not 0 <= value < LOW_MODE_COUNT:
            raise ValueError("numeric low-mode identifier outside the four-edge patch")
        return value
    if isinstance(value, str) and len(value) == 2 and value[0] in {"a", "b"} and value[1].isdigit():
        return rail_mode(int(value[1]), 0 if value[0] == "a" else 1)
    raise ValueError("mode must be an integer 0..7 or a label a0..a3/b0..b3")


def _parse_complex(value: object) -> complex:
    if isinstance(value, (int, float)):
        return complex(value)
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return complex(float(value[0]), float(value[1]))
    raise ValueError("complex coefficient must be a real number or [real, imaginary]")


def candidate_from_payload(payload: dict) -> MicroscopicCandidate2D:
    """Parse the canonical JSON-like candidate format used by an external optimizer.

    This is intentionally a parser, not an unconstrained code-execution hook:
    the resulting dataclass must still pass all U(1), parity-label, locality,
    and perturbative-window validators before an audit can run.
    """
    try:
        parsed_channels = []
        for item in payload["channels"]:
            channel = PairConversionChannel(
                name=str(item["name"]),
                detuning=float(item["detuning"]),
                coupling=float(item["coupling"]),
                phase=float(item.get("phase", 0.0)),
                pair_terms=tuple(
                    (_parse_mode(term[0]), _parse_mode(term[1]), _parse_complex(term[2]))
                    for term in item["pair_terms"]
                ),
            )
            # A solver may echo the signature supplied in its prompt.  When it
            # does, make it an auditable claim rather than ignored decoration.
            if "mediator_parity_signature" in item:
                declared = tuple(int(value) for value in item["mediator_parity_signature"])
                if declared != channel.parity_signature():
                    raise ValueError(
                        f"channel {channel.name} declares mediator parity {declared}, "
                        f"but its pair terms require {channel.parity_signature()}"
                    )
            parsed_channels.append(channel)
        channels = tuple(parsed_channels)
        hops = tuple(RailHop(
            left_edge=int(item["left_edge"]), right_edge=int(item["right_edge"]), rail=int(item["rail"]),
            amplitude=float(item["amplitude"]), phase=float(item.get("phase", 0.0)),
        ) for item in payload.get("rail_hops", ()))
        zz = tuple(ZZCoupling(
            left_edge=int(item["left_edge"]), right_edge=int(item["right_edge"]), strength=float(item["strength"]),
        ) for item in payload.get("zz_couplings", ()))
        candidate = MicroscopicCandidate2D(
            mott_u=float(payload["mott_u"]), channels=channels, rail_hops=hops, zz_couplings=zz,
            rail_biases=tuple(float(value) for value in payload.get("rail_biases", (0.0,) * EDGE_COUNT)),
            maximum_sw_ratio=float(payload.get("maximum_sw_ratio", 0.15)),
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"invalid candidate payload: {error}") from error
    candidate.validate()
    return candidate


def interaction_connectivity(candidate: MicroscopicCandidate2D) -> dict:
    """Return the edge-graph connectivity needed by the four-body curriculum."""
    candidate.validate()
    adjacency = {edge: set() for edge in range(EDGE_COUNT)}
    def join(left: int, right: int) -> None:
        if left != right:
            adjacency[left].add(right)
            adjacency[right].add(left)
    for channel in candidate.channels:
        channel_edges: set[int] = set()
        for first, second, _ in channel.pair_terms:
            join(first // 2, second // 2)
            channel_edges.update((first // 2, second // 2))
        # A mediator shared by otherwise disjoint pair conversions is a local
        # plaquette hyperedge, not two independent links.  Record that physical
        # connection rather than misclassifying the crossed primitive as
        # disconnected.
        for left in channel_edges:
            for right in channel_edges:
                join(left, right)
    for hop in candidate.rail_hops:
        join(hop.left_edge, hop.right_edge)
    for coupling in candidate.zz_couplings:
        join(coupling.left_edge, coupling.right_edge)
    unseen = set(range(EDGE_COUNT))
    components = []
    while unseen:
        root = unseen.pop()
        component, frontier = {root}, [root]
        while frontier:
            edge = frontier.pop()
            new = adjacency[edge] & unseen
            unseen -= new
            component |= new
            frontier.extend(new)
        components.append(sorted(component))
    return {"components": sorted(components), "component_count": len(components), "patch_connected": len(components) == 1}


STATE_VECTOR_FIELDS = (
    "minimum_monomer_overlap_singular_value",
    "low_to_high_gap",
    "gap_over_mott_u",
    "maximum_low_state_mott_violation",
    "parity_commutator_frobenius_max",
    "observed_target_coefficient",
    "target_coefficient_ratio",
    "target_alignment",
    "fixed_scale_spectral_algebraic_residual",
    "best_scale_operator_residual",
    "unwanted_pauli_norm_over_target",
)


def compact_optimizer_observation(audit: dict, maximum_unwanted_terms: int = 4) -> dict:
    """Compress an audit to a stable state vector for a context-limited solver."""
    state = audit["state_vector"]
    return {
        "schema": "antler.phase7.optimizer-observation.v1",
        "target": audit["target"],
        "reward": audit["reward"],
        "state_fields": list(STATE_VECTOR_FIELDS),
        "state_vector": [state[field] for field in STATE_VECTOR_FIELDS],
        "pauli_weight_square_histogram": state["pauli_weight_square_histogram"],
        "interaction_connectivity": interaction_connectivity(candidate_from_payload(audit["candidate"])),
        "top_unwanted_paulis": audit["top_unwanted_paulis"][:maximum_unwanted_terms],
        "hard_failures": audit["hard_failures"],
        "next_action": (
            "repair hard failure before changing target overlap"
            if audit["hard_failures"] else
            "increase connected target coefficient while reducing the listed unwanted Pauli terms"
        ),
    }


def registered_search_space() -> dict:
    """Pre-registered curriculum bounds, expressed without a solver dependency."""
    return {
        "invariants_enforced_by_grammar": [
            "weighted U(1) charge conservation",
            "finite four-edge local support",
            "hard-core charge-two mediators only",
            "definite mediator branch-parity signature",
            "no direct four-edge stabilizer insertion",
        ],
        "box_constraints": {
            "mott_u": [16.0, 64.0],
            "channel_detuning": [8.0, 48.0],
            "absolute_coupling_over_detuning": [0.005, 0.15],
            "absolute_hop_over_mott_u": [0.0, 0.15],
            "absolute_zz_over_mott_u": [0.0, 0.15],
            "absolute_bias_over_mott_u": [0.0, 0.05],
            "phase": [-3.141592653589793, 3.141592653589793],
            "maximum_channel_count": 4,
        },
        "curriculum": [
            {
                "stage": 0,
                "name": "structural and downfolding gate",
                "allowed": "one or two even-parity mediator channels; no rail hopping",
                "advance_when": "parity residual <1e-10, overlap singular value >=0.98, positive low/high gap",
            },
            {
                "stage": 1,
                "name": "connected four-edge virtual hypergraph",
                "allowed": (
                    "up to four channels plus same-rail hops; a mediator may couple multiple pair conversions within "
                    "one four-edge plaquette, but no direct four-body stabilizer is allowed"
                ),
                "advance_when": "no hard failure and a sign-correct target Pauli coefficient exceeds numerical noise",
            },
            {
                "stage": 2,
                "name": "operator selectivity",
                "allowed": "phase and density tuning within the same perturbative box",
                "advance_when": "target alignment >=0.80 and unwanted Pauli norm/target <=0.25 on both XXXX and ZZZZ primitives",
            },
            {
                "stage": 3,
                "name": "analytic and tiling gate",
                "allowed": "only candidates with a repeated local motif",
                "advance_when": "independent SW derivation, linked-cluster cross-check, and a scalable 2D tensor/matrix-free plan",
            },
        ],
    }


def evaluate_local_candidate(
    candidate: MicroscopicCandidate2D,
    target_label: str = "XXXX",
    target_strength: float = 1.0,
    include_pauli_coefficients: bool = False,
) -> dict:
    """Return a dense, compact reward audit against one four-body stabilizer.

    ``target_strength`` is the stabilizer coupling ``J``.  Thus the target
    traceless local contribution from ``(I-S)/2`` has coefficient ``-J/2``.
    The effective Hamiltonian is obtained by exact diagonalization and polar
    alignment of its isolated 16-dimensional monomer manifold.  This is an
    exact finite-block downfolding diagnostic, not a fourth-order analytic SW
    proof.
    """
    if target_strength <= 0.0:
        raise ValueError("target strength must be positive")
    candidate.validate()
    H, states, _, monomers, mott_cost = build_local_candidate_hamiltonian(candidate)
    values, vectors = np.linalg.eigh(H)
    low_values, low_vectors = values[:LOW_MANIFOLD_DIMENSION], vectors[:, :LOW_MANIFOLD_DIMENSION]
    overlap = monomers.conj().T @ low_vectors
    left, singular_values, right_dagger = np.linalg.svd(overlap)
    aligned_frame = left @ right_dagger
    effective = aligned_frame @ np.diag(low_values) @ aligned_frame.conj().T
    effective = 0.5 * (effective + effective.conj().T)
    identity_shift = float(np.trace(effective).real / LOW_MANIFOLD_DIMENSION)
    traceless = effective - identity_shift * np.eye(LOW_MANIFOLD_DIMENSION)
    coefficients = _pauli_coefficients(traceless)
    target_operator = pauli_word(target_label)
    target_coefficient = -0.5 * target_strength
    observed_target = coefficients[target_label]
    target_hamiltonian = target_coefficient * target_operator
    target_norm = float(np.linalg.norm(target_hamiltonian))
    residual_fixed_scale = float(np.linalg.norm(traceless - target_hamiltonian) / target_norm)
    amplitude = float(np.trace(target_operator @ traceless).real / np.trace(target_operator @ target_operator).real)
    best_scale_residual = float(
        np.linalg.norm(traceless - amplitude * target_operator) / max(np.linalg.norm(traceless), 1e-15)
    )
    nonidentity = {label: value for label, value in coefficients.items() if label != "IIII"}
    unwanted_square = sum(value * value for label, value in nonidentity.items() if label != target_label)
    total_square = sum(value * value for value in nonidentity.values())
    target_alignment = float(observed_target * observed_target / total_square) if total_square > 1e-24 else 0.0
    weights = {str(weight): 0.0 for weight in range(1, EDGE_COUNT + 1)}
    for label, value in nonidentity.items():
        weights[str(sum(letter != "I" for letter in label))] += value * value
    parity_a, parity_b = _parity_labels(candidate, states)
    symmetry_residual = max(
        float(np.linalg.norm(H * (parity_a[None, :] - parity_a[:, None]))),
        float(np.linalg.norm(H * (parity_b[None, :] - parity_b[:, None]))),
    )
    low_mott = np.real(np.sum(np.abs(low_vectors) ** 2 * mott_cost[:, None], axis=0))
    low_high_gap = float(values[LOW_MANIFOLD_DIMENSION] - values[LOW_MANIFOLD_DIMENSION - 1])
    capture_defect = max(0.0, 1.0 - float(np.min(singular_values)))
    signed_progress = max(0.0, observed_target / target_coefficient) if target_coefficient < 0.0 else max(0.0, observed_target / target_coefficient)
    unwanted_relative = float(np.sqrt(unwanted_square) / abs(target_coefficient))
    # The reward deliberately favours the right *operator*, not merely an
    # isospectral low block.  Every component remains bounded and interpretable.
    reward = (
        3.0 * target_alignment
        + np.tanh(signed_progress)
        - 1.5 * min(unwanted_relative, 4.0)
        - 2.0 * min(capture_defect / 0.02, 2.0)
        - 1.0 * max(0.0, 0.1 - low_high_gap / candidate.mott_u) / 0.1
    )
    hard_failures = []
    if symmetry_residual > 1e-10:
        hard_failures.append("exact parity symmetry residual")
    if low_high_gap <= 1e-9:
        hard_failures.append("unisolated 16-state monomer manifold")
    if np.min(singular_values) < 0.90:
        hard_failures.append("monomer capture below 0.90")
    if hard_failures:
        reward = min(float(reward), -10.0)
    top_unwanted = sorted(
        ({"label": label, "coefficient": value} for label, value in nonidentity.items() if label != target_label),
        key=lambda row: abs(row["coefficient"]), reverse=True,
    )[:8]
    out = {
        "schema": "antler.phase7.local-microscopic-reward.v1",
        "candidate": candidate_payload(candidate),
        "target": {
            "label": target_label,
            "target_strength": target_strength,
            "required_traceless_pauli_coefficient": target_coefficient,
        },
        "state_vector": {
            "block_dimension": int(len(states)),
            "low_manifold_dimension": LOW_MANIFOLD_DIMENSION,
            "minimum_monomer_overlap_singular_value": float(np.min(singular_values)),
            "mean_monomer_overlap_singular_value": float(np.mean(singular_values)),
            "low_to_high_gap": low_high_gap,
            "gap_over_mott_u": low_high_gap / candidate.mott_u,
            "maximum_low_state_mott_violation": float(np.max(low_mott)),
            "parity_commutator_frobenius_max": symmetry_residual,
            "observed_target_coefficient": observed_target,
            "target_coefficient_ratio": observed_target / target_coefficient,
            "target_alignment": target_alignment,
            "fixed_scale_spectral_algebraic_residual": residual_fixed_scale,
            "best_scale_operator_residual": best_scale_residual,
            "unwanted_pauli_norm_over_target": unwanted_relative,
            "pauli_weight_square_histogram": weights,
        },
        "top_unwanted_paulis": top_unwanted,
        "reward": float(reward),
        "hard_failures": hard_failures,
        "decision": (
            "local candidate survives the structural gate; it may enter the next local-order and tiling audit"
            if not hard_failures else
            "reject before any global tiling: " + "; ".join(hard_failures)
        ),
        "claim_boundary": (
            "This reward evaluates one four-edge local block. It does not prove a controlled Schrieffer-Wolff expansion, "
            "a tiled two-dimensional parent, global topological order, anyon braiding, non-Abelian statistics, or fault tolerance."
        ),
    }
    if include_pauli_coefficients:
        out["full_traceless_pauli_coefficients"] = coefficients
    return out


def computational_budget() -> dict:
    """Exact fixed-charge dimensions that set the optimizer's hierarchy."""
    local_dimensions = {
        str(channel_count): sum(
            comb(channel_count, mediator_count) * comb(LOW_MODE_COUNT, EDGE_COUNT - 2 * mediator_count)
            for mediator_count in range(min(channel_count, EDGE_COUNT // 2) + 1)
        )
        for channel_count in range(1, 5)
    }
    return {
        "local_four_edge_fixed_charge_dimensions": local_dimensions,
        "global_3x3_without_mediators_fixed_charge_dimension": comb(36, 18),
        "global_3x3_with_4_charge_two_mediators_fixed_charge_dimension": sum(
            comb(4, mediator_count) * comb(36, 18 - 2 * mediator_count)
            for mediator_count in range(5)
        ),
        "global_3x3_with_18_charge_two_mediators_fixed_charge_dimension": sum(
            comb(18, mediator_count) * comb(36, 18 - 2 * mediator_count)
            for mediator_count in range(10)
        ),
    }
