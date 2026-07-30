#!/usr/bin/env python3
"""Fine-lowest-K-type theta paths for :math:`Mp(2n, R)`.

The existing :mod:`theta2_pbp` calculator ends at ``O(n,n+1)``.  This
module uses the same VALUE-mode lowest-harmonic transitions but stops at a
metaplectic stage.  For an all-even dual orbit with ``r`` rows it retains
exactly the ``2**r`` concrete twist histories, labels them by subsets of the
rows (numbered from bottom to top), and attaches the type-M painted
bipartition selected by Ma's exact lift ancestry.

The displayed tableau is the actual element of ``PBP_M(O^vee, wp)``.  Its
``special_reference_drc`` is kept separately because the recursive associated
cycle is evaluated on that reference together with ``wp``.
"""

from __future__ import annotations

import contextlib
import io
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from typing import Any, Hashable, Iterable

from combunipotent.LS import char_twist_B, lift_B_M, lift_M_B
from combunipotent.drc import gp_form_M, reg_drc, verify_drc
from combunipotent.drclift import gp_form_B_ext, lift_DRCLS
from standalone import (
    build_pbp_bijection,
    compute_AC,
    dpart2Wrepns_with_wp,
    dpart_to_bipartition,
    drc_shape,
)
from theta2_pbp import (
    CHAR_TWISTS,
    DET_TWISTS,
    Drc,
    MaAncestryEdge,
    MaEndpointKey,
    OrthogonalTransition,
    PBPLadder,
    RowMultiplicityProfile,
    SubsetBits,
    TwistHistory,
    ValueList,
    determinant_twist_wedge,
    format_mp_weight,
    mp_to_o_values,
    o_character_to_mp,
    o_values_to_mp,
    pairs_with_sum,
    row_multiplicity_profile,
    trace_ma_ancestry,
    unitary_characters_o,
    validate_dual_orbit,
)


MpWeight = tuple[Fraction, ...]
AssociatedCycle = tuple[tuple[int, tuple[tuple[int, int], ...]], ...]


@dataclass(frozen=True)
class MetaplecticFineChain:
    """One grouped VALUE chain whose last stage is metaplectic."""

    transitions: tuple[OrthogonalTransition, ...]
    final_mp_total: int
    final_mp_weight: MpWeight

    def __post_init__(self) -> None:
        if not self.transitions:
            raise ValueError("a metaplectic chain must contain an orthogonal lift")
        final_transition = self.transitions[-1]
        if (
            final_transition.mp_total != self.final_mp_total
            or final_transition.mp_weight != self.final_mp_weight
        ):
            raise ValueError("the final transition does not match the final Mp data")
        if len(self.final_mp_weight) != self.final_mp_total // 2:
            raise ValueError("the final Mp weight has the wrong rank")
        if any(abs(value) != Fraction(1, 2) for value in self.final_mp_weight):
            raise ValueError("a final fine Mp weight must have only +/-1/2 entries")

    @property
    def fine_degree(self) -> int:
        """Number of ``-1/2`` entries in the final fine Mp weight."""

        return sum(value == Fraction(-1, 2) for value in self.final_mp_weight)

    def display_steps(
        self,
        history: TwistHistory | None = None,
    ) -> tuple[str, ...]:
        """Return the group-qualified lowest-harmonic chain."""

        if history is not None and len(history) != len(self.transitions):
            raise ValueError("history must select one twist at every O stage")
        steps: list[str] = []
        for index, transition in enumerate(self.transitions):
            selected = None if history is None else history[index]
            steps.extend(transition.display_steps(selected))
        return tuple(steps)

    def history_text(self, history: TwistHistory) -> str:
        """Format the selected twists with their orthogonal groups."""

        if len(history) != len(self.transitions):
            raise ValueError("history must select one twist at every O stage")
        return " -> ".join(
            f"O({transition.p},{transition.q}) {twist_label}"
            for transition, twist_label in zip(self.transitions, history)
        )


@dataclass(frozen=True)
class MetaplecticPaintedBipartition:
    """A type-M extended painted bipartition ``(tau_wp, wp)``."""

    raw_drc: Drc
    primitive_pair_indices: tuple[int, ...]
    special_reference_drc: Drc

    def __post_init__(self) -> None:
        indices = self.primitive_pair_indices
        if any(
            not isinstance(index, int)
            or isinstance(index, bool)
            or index < 0
            for index in indices
        ):
            raise ValueError("primitive-pair indices must be nonnegative integers")
        if tuple(sorted(set(indices))) != indices:
            raise ValueError("primitive-pair indices must be sorted and unique")

    @property
    def p_columns(self) -> tuple[str, ...]:
        return self.raw_drc[0]

    @property
    def q_columns(self) -> tuple[str, ...]:
        return self.raw_drc[1]

    @property
    def primitive_pairs(self) -> tuple[tuple[int, int], ...]:
        """Return type-M primitive pairs in one-based row numbering."""

        return tuple(
            (2 * index + 1, 2 * index + 2)
            for index in self.primitive_pair_indices
        )

    @property
    def parameter_key(self) -> tuple[Drc, tuple[int, ...]]:
        return self.raw_drc, self.primitive_pair_indices

    def associated_cycle(
        self,
        cache: dict[Any, Any] | None = None,
    ) -> AssociatedCycle:
        """Compute the associated cycle from the special reference tableau."""

        terms = compute_AC(
            self.special_reference_drc,
            frozenset(self.primitive_pair_indices),
            "M",
            cache=cache,
        )
        return tuple(
            (int(coefficient), tuple(ils))
            for coefficient, ils in terms
        )

    def one_line(self) -> str:
        primitive_pairs = (
            "empty"
            if not self.primitive_pairs
            else "{" + ", ".join(
                f"({left},{right})" for left, right in self.primitive_pairs
            ) + "}"
        )
        return (
            f"P={self.p_columns!r}, Q={self.q_columns!r}, "
            f"primitive_pairs={primitive_pairs}"
        )


@dataclass(frozen=True)
class MetaplecticConcretePath:
    """One selected twist history and its exact type-M endpoint."""

    number: int
    chain: MetaplecticFineChain
    history: TwistHistory
    subset_bits: SubsetBits
    selected_half_row_sum: int
    painted_bipartition: MetaplecticPaintedBipartition

    @property
    def subset_indices(self) -> tuple[int, ...]:
        return tuple(
            index
            for index, bit in enumerate(self.subset_bits, start=1)
            if bit
        )

    @property
    def subset_label(self) -> str:
        return "{" + ", ".join(map(str, self.subset_indices)) + "}"

    @property
    def name(self) -> str:
        return f"Path {self.subset_label}"

    @property
    def fine_degree(self) -> int:
        return self.chain.fine_degree

    @property
    def final_weight(self) -> MpWeight:
        return self.chain.final_mp_weight

    @property
    def final_weight_text(self) -> str:
        return format_mp_weight(self.final_weight)

    @property
    def display_steps(self) -> tuple[str, ...]:
        return self.chain.display_steps(self.history)

    @property
    def history_text(self) -> str:
        return self.chain.history_text(self.history)


@dataclass(frozen=True)
class MetaplecticCalculation:
    """Complete path calculation for one type-M dual orbit."""

    partition: tuple[int, ...]
    totals: tuple[int, ...]
    tower: tuple[str, ...]
    paths: tuple[MetaplecticConcretePath, ...]

    @property
    def rank(self) -> int:
        return sum(self.partition) // 2


@dataclass(frozen=True)
class _PathSeed:
    chain: MetaplecticFineChain
    history: TwistHistory
    subset_bits: SubsetBits
    selected_half_row_sum: int
    endpoint_key: MaEndpointKey


def derive_metaplectic_chain_totals(
    partition: Iterable[int],
) -> tuple[int, ...]:
    """Return the alternating ``O, Mp, ..., Mp`` group totals.

    An odd number of orbit rows gets one terminal zero before the rows are
    read from bottom to top.  Thus the printed tower always has even length
    and ends at ``Mp(sum(partition))``.
    """

    part = validate_dual_orbit(partition)
    induction_rows = part + ((0,) if len(part) % 2 else ())
    suffix_size = 0
    totals: list[int] = []
    for stage_number, row in enumerate(reversed(induction_rows), start=1):
        suffix_size += row
        totals.append(suffix_size + 1 if stage_number % 2 else suffix_size)
    if len(totals) % 2 or totals[-1] != sum(part):
        raise AssertionError("the derived tower does not end at Mp(|O^vee|)")
    return tuple(totals)


def enumerate_metaplectic_chains(
    chain_totals: tuple[int, ...],
) -> tuple[MetaplecticFineChain, ...]:
    """Enumerate grouped VALUE chains ending in a fine Mp weight."""

    if not chain_totals or len(chain_totals) % 2:
        raise ValueError("an M-ending chain must have positive even length")
    if any(total % 2 == 0 for total in chain_totals[0::2]):
        raise ValueError("every orthogonal total must be odd")
    if any(total % 2 for total in chain_totals[1::2]):
        raise ValueError("every Mp total must be even")

    o_totals = chain_totals[0::2]
    mp_totals = chain_totals[1::2]
    output: list[MetaplecticFineChain] = []

    def continue_from_mp(
        mp_index: int,
        mp_weight: MpWeight,
        transitions: tuple[OrthogonalTransition, ...],
    ) -> None:
        if mp_index == len(mp_totals) - 1:
            if all(abs(value) == Fraction(1, 2) for value in mp_weight):
                output.append(
                    MetaplecticFineChain(
                        transitions=transitions,
                        final_mp_total=mp_totals[-1],
                        final_mp_weight=mp_weight,
                    )
                )
            return

        next_o_total = o_totals[mp_index + 1]
        next_mp_total = mp_totals[mp_index + 1]
        next_mp_rank = next_mp_total // 2
        for p, q in pairs_with_sum(next_o_total):
            try:
                base_left, base_right = mp_to_o_values(p, q, mp_weight)
            except ValueError:
                continue

            grouped_twists: dict[
                tuple[ValueList, ValueList, MpWeight],
                list[str],
            ] = {}
            for det_left, det_right, twist_label in DET_TWISTS:
                try:
                    left = (
                        determinant_twist_wedge(base_left, p)
                        if det_left
                        else base_left
                    )
                    right = (
                        determinant_twist_wedge(base_right, q)
                        if det_right
                        else base_right
                    )
                    next_weight = o_values_to_mp(
                        next_mp_rank,
                        p,
                        q,
                        left,
                        right,
                    )
                except ValueError:
                    continue
                grouped_twists.setdefault(
                    (left, right, next_weight),
                    [],
                ).append(twist_label)

            for (left, right, next_weight), twist_labels in grouped_twists.items():
                transition = OrthogonalTransition(
                    p=p,
                    q=q,
                    base_left=base_left,
                    base_right=base_right,
                    twist_labels=tuple(sorted(twist_labels)),
                    twisted_left=left,
                    twisted_right=right,
                    mp_total=next_mp_total,
                    mp_weight=next_weight,
                )
                continue_from_mp(
                    mp_index + 1,
                    next_weight,
                    transitions + (transition,),
                )

    first_o_total = o_totals[0]
    first_mp_total = mp_totals[0]
    first_mp_rank = first_mp_total // 2
    for p, q in pairs_with_sum(first_o_total):
        for det_left, det_right, twist_label in unitary_characters_o(p, q):
            try:
                first_weight = o_character_to_mp(
                    first_mp_rank,
                    p,
                    q,
                    det_left,
                    det_right,
                )
            except ValueError:
                continue
            transition = OrthogonalTransition(
                p=p,
                q=q,
                base_left=(),
                base_right=(),
                twist_labels=(twist_label,),
                twisted_left=(1,) * (p if det_left else 0),
                twisted_right=(1,) * (q if det_right else 0),
                mp_total=first_mp_total,
                mp_weight=first_weight,
                initial_character=True,
            )
            continue_from_mp(0, first_weight, (transition,))

    return tuple(dict.fromkeys(output))


def metaplectic_histories(
    chain: MetaplecticFineChain,
) -> tuple[TwistHistory, ...]:
    """Select one concrete twist from every grouped O-stage option."""

    return tuple(
        tuple(history)
        for history in product(
            *(transition.twist_labels for transition in chain.transitions)
        )
    )


def _decode_metaplectic_raw_subset_bits(raw_bits: SubsetBits) -> SubsetBits:
    """Decode the rows contributing ``-1/2`` with type-M boundary ``epsilon=0``."""

    row_count = len(raw_bits)
    if not row_count or any(bit not in (0, 1) for bit in raw_bits):
        raise ValueError("raw subset bits must be a nonempty binary tuple")

    raw = (0, *raw_bits)
    subset = [0] * (row_count + 1)
    epsilon = 0
    lower_bound = 0 if row_count % 2 == 0 else 1
    for index in range(row_count - 1, lower_bound, -2):
        subset[index] = raw[index] ^ epsilon
        subset[index + 1] = raw[index + 1] ^ epsilon
        epsilon ^= raw[index] ^ raw[index + 1]
    if row_count % 2:
        subset[1] = raw[1] ^ epsilon
    # VALUE's raw convention records the rows contributing +1/2.  The
    # metaplectic fine degree is normalized here by the exterior-degree
    # convention: k counts the -1/2 entries, so the displayed subset is the
    # complementary set of rows.
    return tuple(1 - bit for bit in subset[1:])


def _raw_pair_from_value_stage(
    transition: OrthogonalTransition,
    twist_label: str,
) -> tuple[int, int] | None:
    if transition.p < transition.q:
        return {"tt": (0, 0), "dt": (1, 0)}.get(twist_label)
    if transition.p > transition.q:
        return {"td": (0, 1), "tt": (1, 1)}.get(twist_label)
    return None


def metaplectic_subset_bits(
    partition: Iterable[int],
    chain: MetaplecticFineChain,
    history: TwistHistory,
) -> SubsetBits:
    """Return bottom-up row bits whose half-row lengths count ``-1/2`` entries."""

    part = validate_dual_orbit(partition)
    row_count = len(part)
    if len(history) != len(chain.transitions):
        raise ValueError("history must select one twist at every O stage")
    expected_stage_count = (row_count + 1) // 2
    if len(chain.transitions) != expected_stage_count:
        raise ValueError("the chain length is incompatible with the orbit rows")

    raw = [0] * (row_count + 1)
    transition_index = 0
    if row_count % 2:
        first = chain.transitions[0]
        if (first.p, first.q) == (0, 1):
            raw[1] = 0
        elif (first.p, first.q) == (1, 0):
            raw[1] = 1
        else:
            raise ValueError("an odd-row M tower must begin with O(0,1) or O(1,0)")
        transition_index = 1
        paired_indices = range(2, row_count, 2)
    else:
        paired_indices = range(1, row_count, 2)

    paired_transitions = chain.transitions[transition_index:]
    if len(paired_transitions) != len(tuple(paired_indices)):
        raise ValueError("the paired theta stages do not match the orbit rows")
    paired_indices = (
        range(2, row_count, 2)
        if row_count % 2
        else range(1, row_count, 2)
    )
    for row_index, transition in zip(paired_indices, paired_transitions):
        twist_label = history[transition_index]
        if twist_label not in transition.twist_labels:
            raise ValueError("history uses a twist unavailable in its VALUE chain")
        pair = _raw_pair_from_value_stage(transition, twist_label)
        if pair is None or abs(transition.p - transition.q) != 1:
            raise ValueError("a paired theta stage has an incompatible branch")
        raw[row_index], raw[row_index + 1] = pair
        transition_index += 1

    subset_bits = _decode_metaplectic_raw_subset_bits(tuple(raw[1:]))
    selected_sum = sum(
        (row // 2) * bit
        for row, bit in zip(reversed(part), subset_bits)
    )
    if selected_sum != chain.fine_degree:
        raise RuntimeError(
            "the M subset decoder disagrees with the final number of -1/2 entries"
        )
    return subset_bits


def load_metaplectic_pbp_ladder(
    partition: Iterable[int],
    chain_totals: tuple[int, ...] | None = None,
) -> PBPLadder:
    """Build Ma's bottom-up ``M/B/.../M`` packet ladder from ``Mp(0)``.

    The lift code needs an odd-length padded dual partition so that its last
    reverse-indexed seed has type M.  One zero is added for an even number of
    input rows and two zeros for an odd number of input rows.
    """

    part = validate_dual_orbit(partition)
    totals = (
        derive_metaplectic_chain_totals(part)
        if chain_totals is None
        else tuple(chain_totals)
    )
    if totals != derive_metaplectic_chain_totals(part):
        raise ValueError("chain totals do not match the metaplectic dual orbit")

    padded = part + ((0,) if len(part) % 2 == 0 else (0, 0))
    drc_stage: dict[Drc, Hashable] = {}
    packet_stage: dict[Hashable, tuple] = {}
    drc_stages: list[dict[Drc, Hashable]] = []
    packet_stages: list[dict[Hashable, tuple]] = []
    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        for index in range(len(padded) - 1, -1, -1):
            target_type = "M" if index % 2 == 0 else "B"
            drc_stage, packet_stage = lift_DRCLS(
                drc_stage,
                packet_stage,
                padded[index:],
                target_type,
                report=False,
                reportann=False,
            )
            drc_stages.append(drc_stage)
            packet_stages.append(packet_stage)

    upstream_log = captured.getvalue()
    if "Collision of keys" in upstream_log:
        raise RuntimeError("Ma's M/B ladder reported a DRC collision")
    if len(drc_stages) != len(totals) + 1:
        raise RuntimeError("the M/B ladder has an unexpected number of stages")
    if len(drc_stages[0]) != 1:
        raise RuntimeError("the Mp(0) seed is not unique")

    zero_drc = next(iter(drc_stages[0]))
    if 2 * gp_form_M(zero_drc) != 0:
        raise RuntimeError("the ladder does not start at Mp(0)")
    for stage_index, total in enumerate(totals, start=1):
        stage = drc_stages[stage_index]
        if not stage:
            raise RuntimeError(f"Ma's ladder stage {stage_index} is empty")
        if stage_index % 2:
            invalid = [
                drc
                for drc in stage
                if not verify_drc(drc, "B")
                or sum(gp_form_B_ext(drc)) != total
            ]
        else:
            invalid = [
                drc
                for drc in stage
                if not verify_drc(reg_drc(drc), "M")
                or 2 * gp_form_M(drc) != total
            ]
        if invalid:
            raise RuntimeError(
                f"Ma's ladder stage {stage_index} has invalid group data"
            )

    return PBPLadder(tuple(drc_stages), tuple(packet_stages))


def _canonical_metaplectic_subset_bits_from_edges(
    part: tuple[int, ...],
    edges: tuple[MaAncestryEdge, ...],
) -> SubsetBits | None:
    """Read raw M-boundary branches from one exact Ma ancestry."""

    row_count = len(part)
    expected_edge_count = row_count if row_count % 2 == 0 else row_count + 1
    if len(edges) != expected_edge_count:
        raise RuntimeError("Ma ancestry length is incompatible with the orbit rows")

    raw = [0] * (row_count + 1)
    if row_count % 2:
        first_form = gp_form_B_ext(edges[0].target_drc)
        if first_form == (0, 1):
            raw[1] = 0
        elif first_form == (1, 0):
            raw[1] = 1
        else:
            return None
        first_pair_edge = 2
        paired_indices = range(2, row_count, 2)
    else:
        first_pair_edge = 0
        paired_indices = range(1, row_count, 2)

    pair_lookup = {
        # (lower signature, B epsilon, following M epsilon): raw bit pair
        (True, 0, 0): (0, 0),
        (True, 1, 1): (1, 0),
        (False, 1, 0): (0, 1),
        (False, 0, 0): (1, 1),
    }
    b_edge_indices = range(first_pair_edge, len(edges), 2)
    if len(tuple(paired_indices)) != len(tuple(b_edge_indices)):
        raise RuntimeError("Ma ancestry does not contain one B/M edge pair per row pair")
    paired_indices = (
        range(2, row_count, 2)
        if row_count % 2
        else range(1, row_count, 2)
    )
    for row_index, edge_index in zip(paired_indices, b_edge_indices):
        b_edge = edges[edge_index]
        m_edge = edges[edge_index + 1]
        p, q = gp_form_B_ext(b_edge.target_drc)
        if abs(p - q) != 1:
            return None
        pair = pair_lookup.get((p < q, b_edge.epsilon, m_edge.epsilon))
        if pair is None:
            return None
        raw[row_index], raw[row_index + 1] = pair

    return _decode_metaplectic_raw_subset_bits(tuple(raw[1:]))


def canonical_metaplectic_subset_bits(
    partition: Iterable[int],
    ladder: PBPLadder,
    final_drc: Drc,
) -> SubsetBits | None:
    """Return the row-subset provenance of a final type-M ladder DRC."""

    part = validate_dual_orbit(partition)
    return _canonical_metaplectic_subset_bits_from_edges(
        part,
        trace_ma_ancestry(ladder, final_drc),
    )


def strict_metaplectic_final_local_system(
    zero_local_system: Hashable,
    chain: MetaplecticFineChain,
    history: TwistHistory,
) -> Hashable:
    """Propagate one concrete path to its final type-M local system."""

    if len(history) != len(chain.transitions):
        raise ValueError("history must select one twist at every O stage")
    local_system = zero_local_system
    for transition, twist_label in zip(chain.transitions, history):
        if twist_label not in transition.twist_labels:
            raise ValueError("history uses a twist unavailable in its VALUE chain")
        base_local_system = lift_M_B(local_system, transition.p, transition.q)
        if not base_local_system:
            raise RuntimeError("Ma's M-to-B representation lift vanished")
        twisted_local_system = char_twist_B(
            base_local_system,
            CHAR_TWISTS[twist_label],
        )
        local_system = lift_B_M(
            twisted_local_system,
            transition.mp_total // 2,
        )
        if not local_system:
            raise RuntimeError("Ma's B-to-M representation lift vanished")
    return local_system


def _validated_metaplectic_pbp(
    part: tuple[int, ...],
    final_drc: Drc,
    pbp_bijection: dict[Drc, tuple[Drc, frozenset[int]]],
    expected_shapes: dict[frozenset[int], tuple[tuple[int, ...], tuple[int, ...]]],
) -> MetaplecticPaintedBipartition:
    actual_drc = reg_drc(final_drc)
    try:
        special_reference, primitive_indices = pbp_bijection[actual_drc]
    except KeyError as error:
        raise RuntimeError("a final Ma DRC is absent from the type-M PBP bijection") from error
    special_reference = reg_drc(special_reference)
    indices = tuple(sorted(primitive_indices))

    if not verify_drc(actual_drc, "M"):
        raise RuntimeError("the selected endpoint is not a valid type-M DRC")
    if 2 * gp_form_M(actual_drc) != sum(part):
        raise RuntimeError("the selected type-M DRC has the wrong Mp rank")
    expected_shape = expected_shapes[frozenset(indices)]
    if drc_shape(actual_drc) != expected_shape:
        raise RuntimeError("the selected type-M DRC has the wrong wp-dependent shape")
    if drc_shape(special_reference) != dpart_to_bipartition(part, "M"):
        raise RuntimeError("the PBP bijection returned the wrong special reference shape")

    return MetaplecticPaintedBipartition(
        raw_drc=actual_drc,
        primitive_pair_indices=indices,
        special_reference_drc=special_reference,
    )


def build_metaplectic_endpoint_map(
    part: tuple[int, ...],
    ladder: PBPLadder,
    required_keys: set[MaEndpointKey],
) -> dict[MaEndpointKey, MetaplecticPaintedBipartition]:
    """Match strict M paths to PBPs by local system and row provenance."""

    pbp_bijection, _table = build_pbp_bijection(part, "M")
    expected_shapes = dpart2Wrepns_with_wp(part, "M")
    endpoints: dict[MaEndpointKey, MetaplecticPaintedBipartition] = {}

    for final_drc, represented_local_system in ladder.drc_stages[-1].items():
        if 2 * gp_form_M(final_drc) != sum(part):
            continue
        subset_bits = canonical_metaplectic_subset_bits(part, ladder, final_drc)
        if subset_bits is None:
            continue
        profile = row_multiplicity_profile(part, subset_bits)
        key = represented_local_system, profile
        if key not in required_keys:
            continue
        pbp = _validated_metaplectic_pbp(
            part,
            final_drc,
            pbp_bijection,
            expected_shapes,
        )
        previous = endpoints.get(key)
        if previous is not None and previous != pbp:
            raise RuntimeError(
                "distinct type-M PBPs have the same strict Ma local system "
                "and row-multiplicity profile"
            )
        endpoints[key] = pbp

    missing = required_keys - endpoints.keys()
    if missing:
        raise RuntimeError(
            "some strict metaplectic paths have no provenance-matched PBP "
            f"endpoint: {len(missing)} missing keys"
        )
    return endpoints


def _subset_mask(bits: SubsetBits) -> int:
    return sum(bit << index for index, bit in enumerate(bits))


def calculate_metaplectic(
    partition: Iterable[int],
) -> MetaplecticCalculation:
    """Calculate all ``2**r`` type-M paths, PBPs, and cycle-ready data."""

    part = validate_dual_orbit(partition)
    totals = derive_metaplectic_chain_totals(part)
    chains = enumerate_metaplectic_chains(totals)
    ladder = load_metaplectic_pbp_ladder(part, totals)
    _zero_drc, zero_local_system = next(iter(ladder.drc_stages[0].items()))

    seeds: list[_PathSeed] = []
    for chain in chains:
        for history in metaplectic_histories(chain):
            subset_bits = metaplectic_subset_bits(part, chain, history)
            selected_sum = sum(
                (row // 2) * bit
                for row, bit in zip(reversed(part), subset_bits)
            )
            final_local_system = strict_metaplectic_final_local_system(
                zero_local_system,
                chain,
                history,
            )
            endpoint_key = (
                final_local_system,
                row_multiplicity_profile(part, subset_bits),
            )
            seeds.append(
                _PathSeed(
                    chain=chain,
                    history=history,
                    subset_bits=subset_bits,
                    selected_half_row_sum=selected_sum,
                    endpoint_key=endpoint_key,
                )
            )

    expected_subsets = set(product((0, 1), repeat=len(part)))
    reached_subsets = {seed.subset_bits for seed in seeds}
    if len(seeds) != 2 ** len(part) or reached_subsets != expected_subsets:
        raise RuntimeError("the M-ending VALUE enumeration is not the full row powerset")

    endpoints = build_metaplectic_endpoint_map(
        part,
        ladder,
        {seed.endpoint_key for seed in seeds},
    )
    ordered_seeds = sorted(seeds, key=lambda seed: _subset_mask(seed.subset_bits))
    paths = tuple(
        MetaplecticConcretePath(
            number=number,
            chain=seed.chain,
            history=seed.history,
            subset_bits=seed.subset_bits,
            selected_half_row_sum=seed.selected_half_row_sum,
            painted_bipartition=endpoints[seed.endpoint_key],
        )
        for number, seed in enumerate(ordered_seeds, start=1)
    )

    tower = ("Mp(0)",) + tuple(
        f"O({total})" if stage_number % 2 else f"Mp({total})"
        for stage_number, total in enumerate(totals, start=1)
    )
    return MetaplecticCalculation(
        partition=part,
        totals=totals,
        tower=tower,
        paths=paths,
    )


__all__ = [
    "AssociatedCycle",
    "MetaplecticCalculation",
    "MetaplecticConcretePath",
    "MetaplecticFineChain",
    "MetaplecticPaintedBipartition",
    "build_metaplectic_endpoint_map",
    "calculate_metaplectic",
    "canonical_metaplectic_subset_bits",
    "derive_metaplectic_chain_totals",
    "enumerate_metaplectic_chains",
    "load_metaplectic_pbp_ladder",
    "metaplectic_histories",
    "metaplectic_subset_bits",
    "strict_metaplectic_final_local_system",
]
