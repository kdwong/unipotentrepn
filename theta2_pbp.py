#!/usr/bin/env python3
"""Fine-lowest-K-type theta chains and painted bipartitions.

Edit ``DUAL_ORBIT`` below and run this file.  The input is a good-parity
nilpotent orbit for a type-B problem, written as an all-even partition of
``2n``.  For example, ``(6, 4, 2, 2)`` is an orbit in ``Sp(14, C)``.

The program has two independent layers:

1. It reproduces the VALUE-mode fine-lowest-K-type chain enumeration in the
   user's ``theta2.py``.  At the final orthogonal stage it computes the two
   exact right-trivial O(p)xO(q)-types of exterior degrees ``k`` and ``p-k``
   separately, and only afterward groups them under the common connected
   highest-weight label ``(1^k|0)``.  These are two distinct degrees unless
   ``p=2k``.
2. It propagates every selected twist by Ma's B/M local-system lift.  Painted
   bipartitions are then attached from the exact backward Ma ancestry of the
   final DRC.  At a balanced pair this ancestry supplies the boundary
   specialization that the direct DRC lift does not retain: histories are
   grouped by their selected multiplicities inside each equal-row block, and
   no representation is discarded merely because several DRCs carry the same
   local system.  Ma's raw final DRC is retained as
   ``tau_wp in PBP(O^vee, wp)``.  The BMSZ bijection is used only to recover
   and check ``wp``; its auxiliary transport to the ``wp=empty`` shape is
   never substituted for the actual tableau.

The final real form is fixed throughout to ``O(n,n+1)``.  The retained paths
have exact right-trivial types ``(1^j 0^(n-j) | 0^(n+1))``; the output groups
``j=k`` and ``j=n-k`` only after their painted bipartitions have been
computed.  The opposite orientation ``O(n+1,n)`` is deliberately not part of
this calculator.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import os
import subprocess
import sys
from dataclasses import dataclass, field
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Hashable, Iterable


def _relaunch_in_project_venv() -> None:
    """Use the prepared project environment for direct script execution.

    VS Code may launch this file with a Microsoft Store/global Python that has
    none of the repository dependencies.  When the local ``.venv`` exists,
    transparently restart the same command with its interpreter before any
    third-party package is imported.
    """

    if __name__ != "__main__":
        return
    project_dir = Path(__file__).resolve().parent
    relative_python = (
        Path("Scripts") / "python.exe"
        if os.name == "nt"
        else Path("bin") / "python"
    )
    project_python = project_dir / ".venv" / relative_python
    if not project_python.is_file():
        return
    try:
        already_using_project_python = (
            Path(sys.executable).resolve() == project_python.resolve()
        )
    except OSError:
        already_using_project_python = False
    if already_using_project_python:
        return
    raise SystemExit(
        subprocess.call(
            [str(project_python), str(Path(__file__).resolve()), *sys.argv[1:]]
        )
    )


_relaunch_in_project_venv()


from combunipotent import test_dpart2drcLS
from combunipotent.LS import char_twist_B, lift_B_M, lift_M_B
from combunipotent.drc import gp_form_M, reg_drc, str_dgms, verify_drc
from combunipotent.drclift import (
    descent_drc,
    gp_form_B_ext,
    lift_extdrc_B_M,
    make_extdrc_B,
    split_extdrc_B,
    twist_M_nonspecial,
)
from standalone import (
    build_pbp_bijection,
    dpart2Wrepns_with_wp,
    dpart_to_bipartition,
    drc_shape,
)


# ---------------------------------------------------------------------------
# USER INPUT: change this one line to your favorite all-even dual orbit.
# ---------------------------------------------------------------------------
DUAL_ORBIT = (6, 4, 2, 2, 2)


MpWeight = tuple[Fraction, ...]
Drc = tuple[tuple[str, ...], tuple[str, ...]]
ValueList = tuple[int, ...]
TwistHistory = tuple[str, ...]
SubsetBits = tuple[int, ...]
RowMultiplicityProfile = tuple[tuple[int, int], ...]
MaEndpointKey = tuple[Hashable, RowMultiplicityProfile]

CHAR_TWISTS = {
    "tt": (1, 1),
    "dt": (-1, 1),
    "td": (1, -1),
    "dd": (-1, -1),
}

DET_TWISTS = (
    (False, False, "tt"),
    (True, False, "dt"),
    (False, True, "td"),
    (True, True, "dd"),
)


@dataclass(frozen=True)
class OrthogonalTransition:
    """One non-final O(p,q) stage followed by a metaplectic lift."""

    p: int
    q: int
    base_left: ValueList
    base_right: ValueList
    twist_labels: tuple[str, ...]
    twisted_left: ValueList
    twisted_right: ValueList
    mp_total: int
    mp_weight: MpWeight
    initial_character: bool = False

    def formatted_steps(self) -> tuple[str, ...]:
        """Return the legacy theta2.py steps used by the parity audit."""

        labels = ",".join(self.twist_labels)
        if self.initial_character:
            return (
                f"O({self.p},{self.q})[{labels}]",
                f"Mp({self.mp_total}){format_mp_weight(self.mp_weight)}",
            )
        return (
            f"O({self.p},{self.q}) base "
            f"L={format_values(self.base_left)} R={format_values(self.base_right)}",
            f"twist[{labels}] -> "
            f"L={format_values(self.twisted_left)} "
            f"R={format_values(self.twisted_right)}",
            f"Mp({self.mp_total}){format_mp_weight(self.mp_weight)}",
        )

    def display_steps(self, twist_label: str | None = None) -> tuple[str, ...]:
        """Display this stage with one concrete twist and its lowest harmonic."""

        if twist_label is not None and twist_label not in self.twist_labels:
            raise ValueError(
                f"twist {twist_label!r} is not available at O({self.p},{self.q})"
            )
        twist_text = format_twist_choice("initial twist", self.twist_labels, twist_label)
        twisted_type = format_o_type(
            self.p,
            self.q,
            self.twisted_left,
            self.twisted_right,
        )
        if self.initial_character:
            return (
                f"{twist_text} -> {twisted_type} [lowest harmonic]",
                f"Mp({self.mp_total}){format_mp_weight(self.mp_weight)}",
            )
        return (
            f"{format_o_type(self.p, self.q, self.base_left, self.base_right)} "
            "[untwisted lowest harmonic]",
            f"{format_twist_choice('twist', self.twist_labels, twist_label)} -> "
            f"{twisted_type} [lowest harmonic]",
            f"Mp({self.mp_total}){format_mp_weight(self.mp_weight)}",
        )


@dataclass(frozen=True)
class FineKChain:
    """A structured version of one chain printed by theta2.py."""

    transitions: tuple[OrthogonalTransition, ...]
    final_p: int
    final_q: int
    final_label: str
    final_base_left: ValueList = ()
    final_base_right: ValueList = ()
    final_twist_labels: tuple[str, ...] = ("tt",)
    final_left: ValueList = ()
    final_right: ValueList = ()

    @property
    def final_left_degree(self) -> int:
        """Exterior degree of the exact final O(p)-type."""

        return len(self.final_left)

    @property
    def final_exact_label(self) -> str:
        """Full-dimensional exact O(p)xO(q) label of this path."""

        degree = self.final_left_degree
        return (
            f"({'1' * degree + '0' * (self.final_p - degree)}|"
            f"{'0' * self.final_q})"
        )

    def __str__(self) -> str:
        return " -> ".join(self.display_steps())

    def display_steps(
        self,
        history: TwistHistory | None = None,
    ) -> tuple[str, ...]:
        """Show lowest harmonics with the twist selected by ``history``.

        Without a concrete history, grouped alternatives are explicitly marked
        as options.  User-facing reports pass a history, so every displayed
        twist is the actual twist used by that realization.
        """

        if history is not None and len(history) != len(self.transitions) + 1:
            raise ValueError(
                "a concrete twist history must contain one label per "
                "orthogonal stage"
            )
        steps: list[str] = []
        for index, transition in enumerate(self.transitions):
            selected = None if history is None else history[index]
            steps.extend(transition.display_steps(selected))
        final_selected = None if history is None else history[-1]
        if (
            final_selected is not None
            and final_selected not in self.final_twist_labels
        ):
            raise ValueError(
                f"final twist {final_selected!r} is not available at "
                f"O({self.final_p},{self.final_q})"
            )
        steps.extend(
            (
                f"{format_o_type(self.final_p, self.final_q, self.final_base_left, self.final_base_right)} "
                "[untwisted lowest harmonic]",
                f"{format_twist_choice('final twist', self.final_twist_labels, final_selected)} -> "
                f"{format_o_type(self.final_p, self.final_q, self.final_left, self.final_right)} "
                f"[connected label {self.final_label}]",
            )
        )
        return tuple(steps)


@dataclass(frozen=True)
class PaintedBipartition:
    """A final type-B SO parameter together with the reached O extension.

    BMSZ writes an extended painted bipartition as
    ``tau-hat = (tau_wp, wp)``.  Here ``extended_drc`` is the actual
    ``tau_wp in PBP(O^vee, wp)`` returned by Ma's lift, including its B+/B-
    tag, and ``primitive_pair_indices`` records ``wp``.  The separate
    ``outer_det_twist`` is not part of the painted bipartition: it records
    which of the two O(p,q) extensions ``pi_tau`` and ``pi_tau tensor det``
    was reached by this exact O(p)xO(q) K-type path.

    ``special_reference_extended_drc`` is only the image under the auxiliary
    bijection to ``PBP(O^vee, empty) x 2^PP``.  It is retained for auditing and
    is deliberately excluded from equality and display identity.
    """

    extended_drc: Drc
    outer_det_twist: bool
    primitive_pair_indices: tuple[int, ...] = ()
    special_reference_extended_drc: Drc | None = field(
        default=None,
        compare=False,
        hash=False,
        repr=False,
    )

    def __post_init__(self) -> None:
        if any(
            not isinstance(index, int)
            or isinstance(index, bool)
            or index < 0
            for index in self.primitive_pair_indices
        ):
            raise ValueError("Primitive-pair indices must be nonnegative integers")
        if tuple(sorted(set(self.primitive_pair_indices))) != (
            self.primitive_pair_indices
        ):
            raise ValueError("Primitive-pair indices must be sorted and unique")

    @property
    def plain_drc(self) -> Drc:
        return split_extdrc_B(self.extended_drc)[0]

    @property
    def gamma(self) -> str:
        return split_extdrc_B(self.extended_drc)[1]

    @property
    def primitive_pairs(self) -> tuple[tuple[int, int], ...]:
        """Return the type-B primitive pairs in 1-based row numbering."""

        return tuple(
            (2 * index + 2, 2 * index + 3)
            for index in self.primitive_pair_indices
        )

    @property
    def so_parameter_key(self) -> tuple[Drc, tuple[int, ...]]:
        """Identity of the SO(p,q) extended painted bipartition."""

        return self.extended_drc, self.primitive_pair_indices

    @property
    def o_parameter_key(self) -> tuple[Drc, tuple[int, ...], bool]:
        """Identity after adjoining the global O(p,q) determinant label."""

        return (*self.so_parameter_key, self.outer_det_twist)

    def painted_one_line(self) -> str:
        """Format only the SO painted-bipartition parameter."""

        p_diagram, q_diagram = self.plain_drc
        primitive_pairs = (
            "empty"
            if not self.primitive_pairs
            else "{" + ", ".join(
                f"({left},{right})" for left, right in self.primitive_pairs
            ) + "}"
        )
        return (
            f"P={p_diagram!r}, Q={q_diagram!r}, "
            f"gamma={self.gamma}, primitive_pairs={primitive_pairs}"
        )

    def one_line(self) -> str:
        outer = "det" if self.outer_det_twist else "trivial"
        return f"{self.painted_one_line()}, BMSZ_outer_epsilon={outer}"


@dataclass
class ChainPBPResult:
    """All concrete twist realizations and final PBPs of one fine-K path."""

    chain: FineKChain
    histories_by_pbp: dict[PaintedBipartition, set[TwistHistory]]


@dataclass(frozen=True)
class ConcreteThetaPath:
    """One theta path with every orthogonal twist concretely selected.

    ``FineKChain`` deliberately groups determinant twists that have the same
    VALUE output.  That grouping is useful internally, but it must never lower
    the user-facing path count: two selected twist histories can reach
    different painted bipartitions even when their lowest K-types coincide.
    """

    chain: FineKChain
    history: TwistHistory
    painted_bipartitions: tuple[PaintedBipartition, ...]


@dataclass(frozen=True)
class PBPLadder:
    """Bottom-up DRC/LS stages and the upstream lift-packet metadata."""

    drc_stages: tuple[dict[Drc, Hashable], ...]
    ls_packet_stages: tuple[dict[Hashable, tuple], ...]


@dataclass(frozen=True)
class MaAncestryEdge:
    """One exact edge in the backward Ma ancestry of a final DRC."""

    stage_index: int
    source_drc: Drc
    source_local_system: Hashable
    target_drc: Drc
    target_local_system: Hashable
    epsilon: int


def validate_dual_orbit(partition: Iterable[int]) -> tuple[int, ...]:
    """Validate a decreasing positive all-even partition."""

    part = tuple(partition)
    if not part:
        raise ValueError("DUAL_ORBIT must contain at least one row")
    if any(not isinstance(row, int) or isinstance(row, bool) for row in part):
        raise ValueError("Every row of DUAL_ORBIT must be an integer")
    if any(row <= 0 for row in part):
        raise ValueError("Every row of DUAL_ORBIT must be positive")
    if any(left < right for left, right in zip(part, part[1:])):
        raise ValueError("DUAL_ORBIT must be weakly decreasing")
    if any(row % 2 for row in part):
        raise ValueError("Good parity for this type-B problem requires all rows even")
    return part


def derive_chain_totals(partition: Iterable[int]) -> tuple[int, ...]:
    """Derive theta2's O/Mp totals from reverse cumulative row sums.

    A zero row is appended when the orbit has an even number of rows.  Reading
    upward from the bottom then gives alternating B/M stages, beginning and
    ending with B.  A B stage of suffix size S has total dimension S+1, while
    an M stage has total dimension S.
    """

    part = validate_dual_orbit(partition)
    induction_rows = part + ((0,) if len(part) % 2 == 0 else ())
    suffix_size = 0
    totals: list[int] = []
    for stage_number, row in enumerate(reversed(induction_rows), start=1):
        suffix_size += row
        if stage_number % 2 == 1:
            totals.append(suffix_size + 1)
        else:
            totals.append(suffix_size)
    if totals[-1] != sum(part) + 1:
        raise AssertionError("The derived tower does not end at SO(|O^vee|+1)")
    return tuple(totals)


def pairs_with_sum(total: int) -> tuple[tuple[int, int], ...]:
    return tuple((left, total - left) for left in range(total + 1))


def unitary_characters_o(p: int, q: int):
    """Characters of O(p)xO(q), omitting a determinant on a zero factor."""

    left_choices = (False, True) if p > 0 else (False,)
    right_choices = (False, True) if q > 0 else (False,)
    for det_left in left_choices:
        for det_right in right_choices:
            label = ("d" if det_left else "t") + ("d" if det_right else "t")
            yield det_left, det_right, label


def format_mp_weight(weight: MpWeight) -> str:
    return "(" + ",".join(str(entry) for entry in weight) + ")"


def format_values(values: ValueList) -> str:
    return "[]" if not values else "[" + ",".join(map(str, values)) + "]"


def format_o_type(
    p: int,
    q: int,
    left: ValueList,
    right: ValueList,
) -> str:
    """Format a full-dimensional lowest-harmonic O(p)xO(q) type."""

    if len(left) > p or len(right) > q:
        raise ValueError("the value lists exceed an orthogonal factor dimension")
    left_label = "".join(map(str, left)) + "0" * (p - len(left))
    right_label = "".join(map(str, right)) + "0" * (q - len(right))
    return f"O({p},{q})({left_label}|{right_label})"


def format_twist_choice(
    name: str,
    available: tuple[str, ...],
    selected: str | None,
) -> str:
    """Distinguish one selected twist from grouped path-level options."""

    if selected is not None:
        return f"{name}[{selected}]"
    labels = ",".join(available)
    if len(available) == 1:
        return f"{name}[{labels}]"
    return f"{name} options[{labels}]"


def is_all_ones(values: ValueList) -> bool:
    return all(value == 1 for value in values)


def determinant_twist_wedge(values: ValueList, dimension: int) -> ValueList:
    """Apply theta2.py's pure-wedge determinant-twist rule."""

    if not is_all_ones(values):
        raise ValueError(
            "the determinant rule is implemented only for pure exterior powers"
        )
    degree = len(values)
    if degree > dimension:
        raise ValueError("wedge degree exceeds the orthogonal dimension")
    return (1,) * (dimension - degree)


def o_character_to_mp(
    rank: int,
    p: int,
    q: int,
    det_left: bool,
    det_right: bool,
) -> MpWeight:
    shift = Fraction(p - q, 2)
    left_degree = p if det_left and p > 0 else 0
    right_degree = q if det_right and q > 0 else 0
    if left_degree + right_degree > rank:
        raise ValueError("the character has no lowest harmonic at this Mp rank")
    values = [0] * rank
    for index in range(left_degree):
        values[index] = 1
    for index in range(right_degree):
        values[rank - right_degree + index] = -1
    return tuple(Fraction(value) + shift for value in values)


def mp_to_o_values(p: int, q: int, weight: MpWeight) -> tuple[ValueList, ValueList]:
    """The VALUE map used in theta2.py, without truncation."""

    shift = Fraction(p - q, 2)
    integral_values: list[int] = []
    for coordinate in weight:
        difference = coordinate - shift
        if difference.denominator != 1:
            raise ValueError("non-integral offset in the Mp-to-O value map")
        integral_values.append(int(difference))
    left = tuple(sorted((value for value in integral_values if value > 0), reverse=True))
    right = tuple(sorted((-value for value in integral_values if value < 0), reverse=True))
    if len(left) > p or len(right) > q:
        raise ValueError("the value lists exceed an orthogonal factor dimension")
    return left, right


def o_values_to_mp(
    rank: int,
    p: int,
    q: int,
    left: ValueList,
    right: ValueList,
) -> MpWeight:
    """The lowest-harmonic O-to-Mp value map used in theta2.py."""

    if len(left) > p or len(right) > q:
        raise ValueError("the value lists exceed an orthogonal factor dimension")
    if len(left) + len(right) > rank:
        raise ValueError("the value lists do not fit in the target Mp rank")
    shift = Fraction(p - q, 2)
    values = [0] * rank
    for index, value in enumerate(left):
        values[index] = value
    for index, value in enumerate(reversed(right)):
        values[rank - len(right) + index] = -value
    return tuple(Fraction(value) + shift for value in values)


def target_labels(final_p: int, final_q: int) -> tuple[str, ...]:
    left_rank, right_rank = final_p // 2, final_q // 2
    return tuple(
        f"({'1' * degree + '0' * (left_rank - degree)}|{'0' * right_rank})"
        for degree in range(left_rank + 1)
    )


def resolve_final_form(
    final_total: int,
    final_form: tuple[int, int] | None = None,
) -> tuple[int, int]:
    """Return the sole supported final signature ``O(n,n+1)``.

    ``final_form`` remains as a compatibility check for internal callers, but
    it may only repeat this convention; the opposite orientation is rejected.
    """

    if not isinstance(final_total, int) or isinstance(final_total, bool):
        raise ValueError("final total must be an integer")
    if final_total <= 0:
        raise ValueError("final total must be positive")
    canonical = final_total // 2, final_total - final_total // 2
    if final_form is None:
        return canonical
    if (
        not isinstance(final_form, tuple)
        or len(final_form) != 2
        or any(not isinstance(value, int) or isinstance(value, bool) for value in final_form)
    ):
        raise ValueError("final_form must be a pair of nonnegative integers")
    p, q = final_form
    if p < 0 or q < 0:
        raise ValueError("final_form entries must be nonnegative")
    if p + q != final_total:
        raise ValueError(
            f"final_form must satisfy p+q={final_total}, received {(p, q)}"
        )
    if (p, q) != canonical:
        raise ValueError(
            "the calculator uses only the ordered final form "
            f"O({canonical[0]},{canonical[1]}) = O(n,n+1)"
        )
    return canonical


def o_wedge_degrees_for_connected_type(p: int, k: int) -> tuple[int, ...]:
    """Return the distinct exact degrees ``k`` and ``p-k`` for ``1^k``.

    This helper only enumerates the two right-trivial O(p)-type targets.  It
    does not infer, duplicate, or otherwise alter the resulting PBP count.
    """

    if not isinstance(p, int) or isinstance(p, bool) or p < 0:
        raise ValueError("p must be a nonnegative integer")
    if not isinstance(k, int) or isinstance(k, bool) or not 0 <= k <= p // 2:
        raise ValueError("k must satisfy 0 <= k <= floor(p/2)")
    complement = p - k
    return (k,) if complement == k else (k, complement)


def final_left_wedge_label(
    final_p: int,
    final_q: int,
    weight: MpWeight,
) -> str | None:
    """Return the untwisted final label, retained as a public helper."""

    try:
        left, right = mp_to_o_values(final_p, final_q, weight)
    except ValueError:
        return None
    return final_left_wedge_values_label(final_p, final_q, left, right)


def final_left_wedge_values_label(
    final_p: int,
    final_q: int,
    left: ValueList,
    right: ValueList,
) -> str | None:
    """Group an exact right-trivial wedge degree under ``k=min(j,p-j)``.

    This function classifies; it does not construct or duplicate a path.  In
    particular an actual degree ``p-k`` path must already have been produced
    by the VALUE/twist calculation before it can enter the degree-``k`` bin.
    """

    if right or not is_all_ones(left):
        return None
    degree = len(left)
    if degree > final_p:
        return None
    degree = min(degree, final_p - degree)
    left_rank, right_rank = final_p // 2, final_q // 2
    return (
        f"({'1' * degree + '0' * (left_rank - degree)}|"
        f"{'0' * right_rank})"
    )


def final_twist_options(
    final_p: int,
    final_q: int,
    weight: MpWeight,
) -> tuple[
    tuple[
        str,
        ValueList,
        ValueList,
        tuple[str, ...],
        ValueList,
        ValueList,
    ],
    ...,
]:
    """Return every computed right-trivial representative of a connected label.

    The VALUE convention represents ``det tensor wedge^j`` as
    ``wedge^(dimension-j)``.  Thus the output bin labelled by degree ``k``
    deliberately groups the *actual outputs* of degrees ``k`` and ``p-k``
    only after both have been constructed.  No PBP or representation count is
    obtained by multiplying a lower-degree result.
    """

    try:
        base_left, base_right = mp_to_o_values(final_p, final_q, weight)
    except ValueError:
        return ()

    grouped_twists: dict[tuple[ValueList, ValueList], list[str]] = {}
    for det_left, det_right, twist_label in DET_TWISTS:
        try:
            left = (
                determinant_twist_wedge(base_left, final_p)
                if det_left
                else base_left
            )
            right = (
                determinant_twist_wedge(base_right, final_q)
                if det_right
                else base_right
            )
        except ValueError:
            continue
        grouped_twists.setdefault((left, right), []).append(twist_label)

    options = []
    for (left, right), twist_labels in grouped_twists.items():
        label = final_left_wedge_values_label(
            final_p,
            final_q,
            left,
            right,
        )
        if label is None:
            continue
        options.append(
            (
                label,
                base_left,
                base_right,
                tuple(sorted(twist_labels)),
                left,
                right,
            )
        )
    return tuple(options)


def enumerate_fine_k_chains(
    chain_totals: tuple[int, ...],
    final_form: tuple[int, int] | None = None,
) -> dict[str, tuple[FineKChain, ...]]:
    """Enumerate theta2.py's chains, retaining structured stage data."""

    if not chain_totals or len(chain_totals) % 2 == 0:
        raise ValueError("The derived chain must begin and end with an O stage")
    if any(total % 2 for total in chain_totals[1::2]):
        raise ValueError("Every Mp total must be even")

    o_totals = chain_totals[0::2]
    mp_totals = chain_totals[1::2]
    final_total = o_totals[-1]
    final_p, final_q = resolve_final_form(final_total, final_form)
    output: dict[str, list[FineKChain]] = {
        label: [] for label in target_labels(final_p, final_q)
    }

    def record_final_options(
        weight: MpWeight,
        transitions: tuple[OrthogonalTransition, ...],
    ) -> None:
        for (
            label,
            base_left,
            base_right,
            twist_labels,
            left,
            right,
        ) in final_twist_options(final_p, final_q, weight):
            output[label].append(
                FineKChain(
                    transitions=transitions,
                    final_p=final_p,
                    final_q=final_q,
                    final_label=label,
                    final_base_left=base_left,
                    final_base_right=base_right,
                    final_twist_labels=twist_labels,
                    final_left=left,
                    final_right=right,
                )
            )

    # One-row orbit: Mp(0) lifts directly to the final balanced O group.
    if len(chain_totals) == 1:
        record_final_options((), ())
        return {label: tuple(chains) for label, chains in output.items()}

    def descend_from_mp(
        mp_index: int,
        mp_weight: MpWeight,
        transitions: tuple[OrthogonalTransition, ...],
    ) -> None:
        next_o_total = o_totals[mp_index + 1]
        is_final_o = mp_index + 1 == len(o_totals) - 1
        signatures = ((final_p, final_q),) if is_final_o else pairs_with_sum(next_o_total)

        for p, q in signatures:
            if is_final_o:
                record_final_options(mp_weight, transitions)
                continue

            try:
                base_left, base_right = mp_to_o_values(p, q, mp_weight)
            except ValueError:
                continue

            next_mp_total = mp_totals[mp_index + 1]
            next_rank = next_mp_total // 2
            grouped_twists: dict[tuple[ValueList, ValueList], list[str]] = {}

            for det_left, det_right, label in DET_TWISTS:
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
                except ValueError:
                    continue
                grouped_twists.setdefault((left, right), []).append(label)

            for (left, right), labels in grouped_twists.items():
                try:
                    next_weight = o_values_to_mp(next_rank, p, q, left, right)
                except ValueError:
                    continue
                transition = OrthogonalTransition(
                    p=p,
                    q=q,
                    base_left=base_left,
                    base_right=base_right,
                    twist_labels=tuple(sorted(labels)),
                    twisted_left=left,
                    twisted_right=right,
                    mp_total=next_mp_total,
                    mp_weight=next_weight,
                )
                descend_from_mp(
                    mp_index + 1,
                    next_weight,
                    transitions + (transition,),
                )

    first_o_total = o_totals[0]
    first_mp_total = mp_totals[0]
    first_mp_rank = first_mp_total // 2
    for p, q in pairs_with_sum(first_o_total):
        for det_left, det_right, label in unitary_characters_o(p, q):
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
                twist_labels=(label,),
                twisted_left=(1,) * (p if det_left else 0),
                twisted_right=(1,) * (q if det_right else 0),
                mp_total=first_mp_total,
                mp_weight=first_weight,
                initial_character=True,
            )
            descend_from_mp(0, first_weight, (transition,))

    # Structured equality gives an honest global deduplication, unlike a final
    # comparison of formatted strings.
    return {
        label: tuple(dict.fromkeys(chains))
        for label, chains in output.items()
    }


def load_pbp_ladder(
    partition: tuple[int, ...],
    chain_totals: tuple[int, ...],
) -> PBPLadder:
    """Load the upstream B/M ladder in bottom-up order, beginning at Mp(0)."""

    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        drc_stages, ls_packet_stages = test_dpart2drcLS(
            partition,
            "B",
            test=False,
            report=False,
            reportann=False,
        )
    upstream_log = captured.getvalue()
    if "Missing " in upstream_log:
        raise RuntimeError(
            "The upstream PBP/local-system ladder reported missing local systems:\n"
            + upstream_log
        )

    # The last entry is the empty pre-induction seed, not a group stage.
    stages = tuple(reversed(drc_stages[:-1]))
    packet_stages = tuple(reversed(ls_packet_stages[:-1]))
    expected_stage_count = len(chain_totals) + 1  # Mp(0) plus printed stages
    if len(stages) != expected_stage_count or len(packet_stages) != expected_stage_count:
        raise RuntimeError(
            f"Expected {expected_stage_count} PBP stages, received "
            f"{len(stages)} DRC stages and {len(packet_stages)} packet stages"
        )
    if len(stages[0]) != 1:
        raise RuntimeError("The Mp(0) stage is not unique")

    for index, total in enumerate(chain_totals, start=1):
        stage = stages[index]
        if index % 2 == 1:
            bad = [drc for drc in stage if sum(gp_form_B_ext(drc)) != total]
        else:
            bad = [drc for drc in stage if 2 * gp_form_M(drc) != total]
        if bad:
            raise RuntimeError(f"PBP stage {index} has an unexpected group size")
    return PBPLadder(stages, packet_stages)


def trace_ma_ancestry(
    ladder: PBPLadder,
    final_drc: Drc,
) -> tuple[MaAncestryEdge, ...]:
    """Trace one raw final DRC back to ``Mp(0)`` through exact packet edges.

    Local systems need not identify a unique DRC at a balanced pair.  Ma's
    packet records do identify the edge once both the target DRC and its local
    system are fixed, so no list-order or singleton-packet heuristic is used.
    """

    try:
        target_local_system = ladder.drc_stages[-1][final_drc]
    except KeyError as error:
        raise ValueError("the requested DRC is not in the final ladder stage") from error

    target_drc = final_drc
    reversed_edges: list[MaAncestryEdge] = []
    for stage_index in range(len(ladder.drc_stages) - 1, 0, -1):
        packet = ladder.ls_packet_stages[stage_index].get(target_local_system)
        entries = () if packet is None else packet[1]
        matches = [
            entry
            for entry in entries
            if entry[0] == target_drc and entry[1] == target_local_system
        ]
        if len(matches) != 1:
            raise RuntimeError(
                "Ma ancestry is not unique for an exact target DRC/local-system "
                f"pair at ladder stage {stage_index}: found {len(matches)} edges"
            )

        (
            recorded_target_drc,
            recorded_target_local_system,
            source_drc,
            epsilon,
            source_local_system,
        ) = matches[0][:5]
        if epsilon not in (0, 1, False, True):
            raise RuntimeError("Ma packet epsilon is not binary")
        if (
            ladder.drc_stages[stage_index - 1].get(source_drc)
            != source_local_system
        ):
            raise RuntimeError(
                f"Ma ancestry source is absent from ladder stage {stage_index - 1}"
            )
        reversed_edges.append(
            MaAncestryEdge(
                stage_index=stage_index,
                source_drc=source_drc,
                source_local_system=source_local_system,
                target_drc=recorded_target_drc,
                target_local_system=recorded_target_local_system,
                epsilon=int(epsilon),
            )
        )
        target_drc = source_drc
        target_local_system = source_local_system

    edges = tuple(reversed(reversed_edges))
    if tuple(edge.stage_index for edge in edges) != tuple(
        range(1, len(ladder.drc_stages))
    ):
        raise RuntimeError("Ma ancestry does not traverse every ladder stage")
    return edges


def _decode_raw_subset_bits(raw_bits: SubsetBits) -> SubsetBits:
    """Apply the backward XOR decoder from raw theta choices to subset bits."""

    row_count = len(raw_bits)
    if not row_count or any(bit not in (0, 1) for bit in raw_bits):
        raise ValueError("raw subset bits must be a nonempty binary tuple")

    raw = (0, *raw_bits)  # one-based indexing for the formulas
    subset = [0] * (row_count + 1)
    epsilon = raw[row_count]
    subset[row_count] = epsilon
    lower_bound = 0 if row_count % 2 else 1
    for index in range(row_count - 2, lower_bound, -2):
        subset[index] = raw[index] ^ epsilon
        subset[index + 1] = raw[index + 1] ^ epsilon
        epsilon ^= raw[index] ^ raw[index + 1]
    if row_count % 2 == 0:
        subset[1] = raw[1] ^ epsilon
    return tuple(subset[1:])


def theta_subset_bits(
    partition: Iterable[int],
    chain: FineKChain,
    history: TwistHistory,
) -> SubsetBits:
    """Decode a concrete VALUE history into its subset of orbit rows.

    Rows are numbered from the bottom: for
    ``O^vee=(2a_r,...,2a_1)``, the returned bit ``b_i`` records whether row
    ``a_i`` is selected.  The decoder remains valid when some ``a_i`` agree;
    in that case only the selected multiplicity in the equal-row block is an
    invariant of the boundary representation.
    """

    part = validate_dual_orbit(partition)
    row_count = len(part)
    expected_history_length = len(chain.transitions) + 1
    if len(history) != expected_history_length:
        raise ValueError("history must select one twist at every orthogonal stage")
    expected_orthogonal_stages = (
        (row_count + 1) // 2 if row_count % 2 else row_count // 2 + 1
    )
    if expected_history_length != expected_orthogonal_stages:
        raise ValueError("the chain length is incompatible with the orbit rows")

    raw = [0] * (row_count + 1)
    transition_index = 0
    if row_count % 2 == 0:
        first = chain.transitions[0]
        if (first.p, first.q) == (0, 1):
            raw[1] = 0
        elif (first.p, first.q) == (1, 0):
            raw[1] = 1
        else:
            raise ValueError("the even-row tower does not begin with O(0,1) or O(1,0)")
        transition_index = 1
        paired_indices = range(2, row_count - 1, 2)
    else:
        paired_indices = range(1, row_count - 1, 2)

    paired_transitions = chain.transitions[transition_index:]
    if len(paired_transitions) != len(tuple(paired_indices)):
        raise ValueError("the paired theta stages do not match the orbit rows")
    # Recreate the range after using it for the length check.
    paired_indices = (
        range(2, row_count - 1, 2)
        if row_count % 2 == 0
        else range(1, row_count - 1, 2)
    )
    for row_index, transition in zip(paired_indices, paired_transitions):
        twist_label = history[transition_index]
        if transition.p < transition.q:
            pair = {"tt": (0, 0), "dt": (1, 0)}.get(twist_label)
        elif transition.p > transition.q:
            pair = {"td": (0, 1), "tt": (1, 1)}.get(twist_label)
        else:
            pair = None
        if pair is None or abs(transition.p - transition.q) != 1:
            raise ValueError(
                "a paired theta stage has an incompatible signature/twist choice"
            )
        raw[row_index], raw[row_index + 1] = pair
        transition_index += 1

    try:
        raw[row_count] = {"tt": 0, "dt": 1}[history[-1]]
    except KeyError as error:
        raise ValueError("the final twist must be tt or dt") from error

    subset_bits = _decode_raw_subset_bits(tuple(raw[1:]))
    half_rows = tuple(row // 2 for row in reversed(part))
    weighted_degree = sum(
        half_row * bit for half_row, bit in zip(half_rows, subset_bits)
    )
    if weighted_degree != chain.final_left_degree:
        raise RuntimeError(
            "the theta subset decoder disagrees with the exact final exterior degree"
        )
    return subset_bits


def row_multiplicity_profile(
    partition: Iterable[int],
    subset_bits: SubsetBits,
) -> RowMultiplicityProfile:
    """Forget labels inside equal-row blocks but retain selected multiplicity."""

    part = validate_dual_orbit(partition)
    if len(subset_bits) != len(part) or any(bit not in (0, 1) for bit in subset_bits):
        raise ValueError("subset bits must be binary with one entry per orbit row")
    counts: dict[int, int] = {}
    for row, bit in zip(reversed(part), subset_bits):
        half_row = row // 2
        counts[half_row] = counts.get(half_row, 0) + bit
    return tuple(sorted(counts.items()))


def strict_ma_final_local_system(
    zero_local_system: Hashable,
    chain: FineKChain,
    history: TwistHistory,
) -> Hashable:
    """Propagate one concrete history by Ma's representation lift only."""

    if len(history) != len(chain.transitions) + 1:
        raise ValueError("history must select one twist at every orthogonal stage")
    local_system = zero_local_system
    for transition, twist_label in zip(chain.transitions, history):
        if twist_label not in transition.twist_labels:
            raise ValueError("history uses a twist unavailable in its VALUE chain")
        base_local_system = lift_M_B(local_system, transition.p, transition.q)
        if not base_local_system:
            raise RuntimeError("Ma's M-to-B representation lift vanished")
        input_local_system = char_twist_B(
            base_local_system,
            CHAR_TWISTS[twist_label],
        )
        local_system = lift_B_M(
            input_local_system,
            transition.mp_total // 2,
        )
        if not local_system:
            raise RuntimeError("Ma's B-to-M representation lift vanished")

    final_twist = history[-1]
    if final_twist not in chain.final_twist_labels:
        raise ValueError("history uses a final twist unavailable in its VALUE chain")
    base_final_local_system = lift_M_B(
        local_system,
        chain.final_p,
        chain.final_q,
    )
    if not base_final_local_system:
        raise RuntimeError("Ma's final M-to-B representation lift vanished")
    return char_twist_B(
        base_final_local_system,
        CHAR_TWISTS[final_twist],
    )


def _canonical_ma_subset_bits_from_edges(
    part: tuple[int, ...],
    edges: tuple[MaAncestryEdge, ...],
) -> SubsetBits | None:
    row_count = len(part)
    expected_edge_count = row_count if row_count % 2 else row_count + 1
    if len(edges) != expected_edge_count:
        raise RuntimeError("Ma ancestry length is incompatible with the orbit rows")

    raw = [0] * (row_count + 1)
    first_paired_stage = 1
    if row_count % 2 == 0:
        first_form = gp_form_B_ext(edges[0].target_drc)
        if first_form == (0, 1):
            raw[1] = 0
        elif first_form == (1, 0):
            raw[1] = 1
        else:
            return None
        first_paired_stage = 3
        paired_indices = range(2, row_count - 1, 2)
    else:
        paired_indices = range(1, row_count - 1, 2)

    pair_lookup = {
        # (lower signature, B epsilon, following M epsilon): raw bit pair
        (True, 0, 0): (0, 0),
        (True, 1, 1): (1, 0),
        (False, 1, 0): (0, 1),
        (False, 0, 0): (1, 1),
    }
    b_stage_indices = range(first_paired_stage, len(edges), 2)
    for row_index, stage_index in zip(paired_indices, b_stage_indices):
        b_edge = edges[stage_index - 1]
        m_edge = edges[stage_index]
        p, q = gp_form_B_ext(b_edge.target_drc)
        if abs(p - q) != 1:
            return None
        pair = pair_lookup.get((p < q, b_edge.epsilon, m_edge.epsilon))
        if pair is None:
            return None
        raw[row_index], raw[row_index + 1] = pair

    raw[row_count] = edges[-1].epsilon
    return _decode_raw_subset_bits(tuple(raw[1:]))


def canonical_ma_subset_bits(
    partition: Iterable[int],
    ladder: PBPLadder,
    final_drc: Drc,
) -> SubsetBits | None:
    """Read the distinct-row boundary provenance encoded by one Ma ancestry.

    ``None`` means that the final DRC ancestry is not one of the right-trivial
    VALUE histories considered by this calculator.  Exact-ancestry failures
    remain hard errors.
    """

    part = validate_dual_orbit(partition)
    edges = trace_ma_ancestry(ladder, final_drc)
    return _canonical_ma_subset_bits_from_edges(part, edges)


def resolve_type_b_pbp(
    stage: dict[Drc, Hashable],
    target_local_system: Hashable,
    form: tuple[int, int],
    source_drc: Drc,
) -> tuple[tuple[Drc, bool], ...]:
    """Legacy diagnostic for exact M-to-B DRC edges.

    Form and local system alone can collide in a large packet.  Exact DRC
    descent is preferred.  A cross-source coincidence is accepted only when
    form plus local system leaves one PBP/outer parameter.  This helper is not
    used to decide whether a concrete representation path exists.
    """

    all_matches = []
    source_matches = []
    for drc, represented_local_system in stage.items():
        if gp_form_B_ext(drc) != form:
            continue
        matches_source = reg_drc(descent_drc(drc, "B")) == reg_drc(source_drc)
        if represented_local_system == target_local_system:
            match = (drc, False)
            all_matches.append(match)
            if matches_source:
                source_matches.append(match)
        if (
            char_twist_B(represented_local_system, CHAR_TWISTS["dd"])
            == target_local_system
        ):
            match = (drc, True)
            all_matches.append(match)
            if matches_source:
                source_matches.append(match)

    exact = tuple(dict.fromkeys(source_matches))
    if exact:
        return exact
    # This singleton fallback is useful only for comparing legacy DRC output.
    fallback = tuple(dict.fromkeys(all_matches))
    return fallback if len(fallback) == 1 else ()


def resolve_type_m_lift(
    stage: dict[Drc, Hashable],
    packet_stage: dict[Hashable, tuple],
    target_local_system: Hashable,
    source_matches: tuple[tuple[Drc, bool], ...],
    added_column_length: int,
) -> tuple[Drc, ...]:
    """Legacy diagnostic for B-to-M DRC edges and singleton LS fallback."""

    exact_targets = {
        drc
        for drc, represented_local_system in stage.items()
        if represented_local_system == target_local_system
    }
    direct_matches = []
    for source_drc, outer_det in source_matches:
        try:
            candidate = lift_extdrc_B_M(source_drc, added_column_length)
            if candidate is not None and outer_det:
                candidate = twist_M_nonspecial(candidate)
        except (AssertionError, IndexError, TypeError, ValueError):
            candidate = None
        if candidate in exact_targets:
            direct_matches.append(candidate)
    packet = packet_stage.get(target_local_system)
    entries = () if packet is None else packet[1]
    source_set = {
        (reg_drc(source_drc), int(outer_det))
        for source_drc, outer_det in source_matches
    }
    matches = []
    for target_drc, represented_ls, source_drc, epsilon, _ in entries:
        if represented_ls != target_local_system:
            continue
        if (reg_drc(source_drc), epsilon) not in source_set:
            continue
        if stage.get(target_drc) != target_local_system:
            continue
        matches.append(target_drc)
    exact_matches = tuple(dict.fromkeys((*direct_matches, *matches)))
    if exact_matches:
        return exact_matches
    # As at B stages, a unique LS target certifies a cross-source coincidence;
    # a multi-PBP LS packet is never assigned heuristically.
    fallback = tuple(exact_targets)
    return fallback if len(fallback) == 1 else ()


def compute_chain_pbps(
    chain: FineKChain,
    ladder: PBPLadder,
) -> ChainPBPResult:
    """Run the former forward DRC traversal for low-level diagnostics.

    The public calculation uses :func:`compute_strict_chain_pbps`; this older
    routine must not be used as a gate at balanced pairs, where a valid Ma
    representation lift can have several DRCs with the same local system.
    """

    stages = ladder.drc_stages
    packet_stages = ladder.ls_packet_stages
    zero_drc, zero_local_system = next(iter(stages[0].items()))
    # A state is keyed by its Mp PBP and local system; its value remembers every
    # concrete sequence of character twists that reached that state.
    states: dict[tuple[Drc, Hashable], set[TwistHistory]] = {
        (zero_drc, zero_local_system): {()}
    }

    for transition_index, transition in enumerate(chain.transitions):
        b_stage = stages[1 + 2 * transition_index]
        m_stage = stages[2 + 2 * transition_index]
        m_packet_stage = packet_stages[2 + 2 * transition_index]
        next_states: dict[tuple[Drc, Hashable], set[TwistHistory]] = {}

        for (source_drc, source_local_system), histories in states.items():
            base_local_system = lift_M_B(
                source_local_system,
                transition.p,
                transition.q,
            )
            if not base_local_system:
                continue

            for twist_label in transition.twist_labels:
                input_local_system = char_twist_B(
                    base_local_system,
                    CHAR_TWISTS[twist_label],
                )
                b_matches = resolve_type_b_pbp(
                    b_stage,
                    input_local_system,
                    (transition.p, transition.q),
                    source_drc,
                )
                if not b_matches:
                    continue

                next_local_system = lift_B_M(
                    input_local_system,
                    transition.mp_total // 2,
                )
                if not next_local_system:
                    continue
                m_matches = resolve_type_m_lift(
                    m_stage,
                    m_packet_stage,
                    next_local_system,
                    b_matches,
                    (
                        transition.mp_total
                        - (transition.p + transition.q - 1)
                    )
                    // 2,
                )
                if not m_matches:
                    continue

                # The B match certifies that this twist and descent path occurs.
                # Packet source/epsilon data identify the next state even when
                # several M diagrams happen to carry the same local system.
                for next_drc in m_matches:
                    destination = (next_drc, next_local_system)
                    destination_histories = next_states.setdefault(destination, set())
                    for history in histories:
                        destination_histories.add(history + (twist_label,))

        if not next_states:
            return ChainPBPResult(chain, {})
        states = next_states

    final_stage = stages[-1]
    histories_by_pbp: dict[PaintedBipartition, set[TwistHistory]] = {}
    for (source_drc, source_local_system), histories in states.items():
        base_final_local_system = lift_M_B(
            source_local_system,
            chain.final_p,
            chain.final_q,
        )
        for final_twist_label in chain.final_twist_labels:
            final_local_system = char_twist_B(
                base_final_local_system,
                CHAR_TWISTS[final_twist_label],
            )
            matches = resolve_type_b_pbp(
                final_stage,
                final_local_system,
                (chain.final_p, chain.final_q),
                source_drc,
            )
            if not matches:
                continue
            for final_drc, outer_det in matches:
                if not verify_drc(final_drc, "B"):
                    raise RuntimeError(
                        f"Invalid final painted bipartition: {final_drc!r}"
                    )
                pbp = PaintedBipartition(final_drc, outer_det)
                destination_histories = histories_by_pbp.setdefault(pbp, set())
                for history in histories:
                    destination_histories.add(history + (final_twist_label,))

    return ChainPBPResult(chain, histories_by_pbp)


def concrete_theta_paths(
    results: Iterable[ChainPBPResult],
) -> tuple[ConcreteThetaPath, ...]:
    """Expand grouped chain objects into user-facing concrete theta paths."""

    numbered_results = tuple(enumerate(results))
    paths = []
    for _result_number, result in sorted(
        numbered_results,
        key=lambda item: (item[1].chain.final_left_degree, item[0]),
    ):
        histories = sorted(
            {
                history
                for pbp_histories in result.histories_by_pbp.values()
                for history in pbp_histories
            }
        )
        for history in histories:
            painted_bipartitions = tuple(
                sorted(
                    (
                        pbp
                        for pbp, pbp_histories in result.histories_by_pbp.items()
                        if history in pbp_histories
                    ),
                    key=lambda pbp: (
                        pbp.extended_drc,
                        pbp.primitive_pair_indices,
                        pbp.outer_det_twist,
                    ),
                )
            )
            paths.append(
                ConcreteThetaPath(
                    chain=result.chain,
                    history=history,
                    painted_bipartitions=painted_bipartitions,
                )
            )
    return tuple(paths)


def candidate_concrete_path_count(chains: Iterable[FineKChain]) -> int:
    """Count concrete histories after selecting every grouped twist option."""

    total = 0
    for chain in chains:
        history_count = len(chain.final_twist_labels)
        for transition in chain.transitions:
            history_count *= len(transition.twist_labels)
        total += history_count
    return total


def orbit_pbp_shape(partition: Iterable[int]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Return the ``wp=empty`` reference shape read from ``O^vee``.

    For type B and ``O^vee = (r_1, r_2, ...)``, these are column lengths
    ``P = (r_2/2, r_4/2, ...)`` and
    ``Q = (r_1/2, r_3/2, ...)``.

    Nonempty ``wp`` generally changes this shape; use ``orbit_pbp_shapes``
    for the complete BMSZ disjoint union.
    """

    part = validate_dual_orbit(partition)
    return dpart_to_bipartition(part, "B")


def orbit_pbp_shapes(
    partition: Iterable[int],
) -> dict[tuple[int, ...], tuple[tuple[int, ...], tuple[int, ...]]]:
    """Return the allowed raw shape of ``PBP(O^vee, wp)`` for every ``wp``."""

    part = validate_dual_orbit(partition)
    return {
        tuple(sorted(wp)): shape
        for wp, shape in dpart2Wrepns_with_wp(part, "B").items()
    }


def attach_wp_label(
    partition: tuple[int, ...],
    raw_pbp: PaintedBipartition,
    pbp_bijection: dict[Drc, tuple[Drc, frozenset[int]]],
    expected_shapes: dict[
        tuple[int, ...], tuple[tuple[int, ...], tuple[int, ...]]
    ] | None = None,
) -> PaintedBipartition:
    """Attach and check the BMSZ primitive-pair label of one raw Ma PBP."""

    shapes = orbit_pbp_shapes(partition) if expected_shapes is None else expected_shapes
    raw_plain_drc = reg_drc(raw_pbp.plain_drc)
    try:
        special_plain_drc, wp = pbp_bijection[raw_plain_drc]
    except KeyError as error:
        raise RuntimeError(
            "A final theta-lift DRC is absent from the orbit's BMSZ "
            f"PBP bijection for O^vee={partition}: {raw_plain_drc!r}"
        ) from error

    wp_indices = tuple(sorted(wp))
    expected_shape = shapes.get(wp_indices)
    actual_shape = drc_shape(raw_plain_drc)
    if expected_shape is None or actual_shape != expected_shape:
        raise RuntimeError(
            "Raw PBP shape is incompatible with its primitive-pair set: "
            f"O^vee={partition}, wp={wp_indices}, "
            f"expected={expected_shape}, actual={actual_shape}, "
            f"raw={raw_plain_drc!r}"
        )

    try:
        special_extended_drc = make_extdrc_B(
            special_plain_drc,
            raw_pbp.gamma,
        )
    except AssertionError as error:
        raise RuntimeError(
            "The auxiliary wp-empty reference is incompatible with its "
            f"B+/B- tag for O^vee={partition}: "
            f"drc={special_plain_drc!r}, gamma={raw_pbp.gamma}"
        ) from error

    if gp_form_B_ext(special_extended_drc) != gp_form_B_ext(
        raw_pbp.extended_drc
    ):
        raise RuntimeError(
            "Transport to the auxiliary wp-empty reference changed the "
            f"final real form for O^vee={partition}"
        )

    return PaintedBipartition(
        extended_drc=raw_pbp.extended_drc,
        outer_det_twist=raw_pbp.outer_det_twist,
        primitive_pair_indices=wp_indices,
        special_reference_extended_drc=special_extended_drc,
    )


def attach_wp_labels(
    partition: tuple[int, ...],
    raw_result: ChainPBPResult,
    pbp_bijection: dict[Drc, tuple[Drc, frozenset[int]]],
) -> ChainPBPResult:
    """Certify raw Ma DRCs and attach their BMSZ primitive-pair set ``wp``.

    ``test_dpart2drcLS`` deliberately contains the disjoint union over all
    primitive-pair subsets ``wp``.  A nonempty ``wp`` generally has a shifted
    shape prescribed by BMSZ (8.9).  The repository bijection identifies the
    disjoint union with ``PBP(O^vee, empty) x 2^PP``; here it is used only to
    recover/certify ``wp``.  The domain element ``tau_wp`` remains the PBP
    displayed and grouped by this program.
    """

    expected_shapes = orbit_pbp_shapes(partition)
    histories_by_pbp: dict[PaintedBipartition, set[TwistHistory]] = {}

    for raw_pbp, histories in raw_result.histories_by_pbp.items():
        pbp = attach_wp_label(
            partition,
            raw_pbp,
            pbp_bijection,
            expected_shapes,
        )
        histories_by_pbp.setdefault(pbp, set()).update(histories)

    return ChainPBPResult(raw_result.chain, histories_by_pbp)


# Backward-compatible name for callers of the earlier prototype.  The
# operation no longer canonicalizes/replaces the displayed tableau.
canonicalize_chain_pbps = attach_wp_labels


def concrete_histories(chain: FineKChain) -> tuple[TwistHistory, ...]:
    """Select one available twist at every orthogonal stage of a VALUE chain."""

    choices = [transition.twist_labels for transition in chain.transitions]
    choices.append(chain.final_twist_labels)
    return tuple(tuple(history) for history in product(*choices))


def build_ma_endpoint_map(
    partition: tuple[int, ...],
    ladder: PBPLadder,
    pbp_bijection: dict[Drc, tuple[Drc, frozenset[int]]],
    final_form: tuple[int, int],
    required_keys: set[MaEndpointKey],
) -> dict[MaEndpointKey, PaintedBipartition]:
    """Match strict Ma representations to PBPs by boundary provenance.

    The key consists of the final Ma local system and, for every distinct row
    size, the number of selected copies.  The second component is exactly the
    specialization invariant obtained by permuting equal rows.  It separates
    the two genuine PBPs that can share one local system in an all-balanced
    even tower, while identifying all strict-row branches that coalesce at the
    boundary.
    """

    expected_shapes = orbit_pbp_shapes(partition)
    endpoints: dict[MaEndpointKey, PaintedBipartition] = {}
    final_stage = ladder.drc_stages[-1]
    for final_drc, represented_local_system in final_stage.items():
        if gp_form_B_ext(final_drc) != final_form:
            continue
        edges = trace_ma_ancestry(ladder, final_drc)
        subset_bits = _canonical_ma_subset_bits_from_edges(partition, edges)
        if subset_bits is None:
            continue
        profile = row_multiplicity_profile(partition, subset_bits)

        # At the final B edge, Ma's epsilon records which O extension occurs.
        # Intermediate packet epsilons have different roles and are used only
        # by the pair table in _canonical_ma_subset_bits_from_edges.
        outer_det_twist = bool(edges[-1].epsilon)
        final_local_system = (
            char_twist_B(represented_local_system, CHAR_TWISTS["dd"])
            if outer_det_twist
            else represented_local_system
        )
        final_base = lift_M_B(
            edges[-1].source_local_system,
            *final_form,
        )
        ancestry_final_local_system = char_twist_B(
            final_base,
            CHAR_TWISTS["dt" if edges[-1].epsilon else "tt"],
        )
        if ancestry_final_local_system != final_local_system:
            raise RuntimeError(
                "The final Ma ancestry epsilon disagrees with its O extension"
            )
        key = final_local_system, profile
        if key not in required_keys:
            continue
        if not verify_drc(final_drc, "B"):
            raise RuntimeError(f"Invalid final painted bipartition: {final_drc!r}")
        pbp = attach_wp_label(
            partition,
            PaintedBipartition(final_drc, outer_det_twist),
            pbp_bijection,
            expected_shapes,
        )
        previous = endpoints.get(key)
        if previous is not None and previous != pbp:
            raise RuntimeError(
                "Distinct PBPs have the same strict Ma local system and "
                f"boundary row-multiplicity profile for O^vee={partition}"
            )
        endpoints[key] = pbp

    missing = required_keys - endpoints.keys()
    if missing:
        raise RuntimeError(
            "Some strict Ma representations have no boundary-provenance PBP "
            f"endpoint for O^vee={partition}: {len(missing)} missing keys"
        )
    return endpoints


def compute_strict_chain_pbps(
    partition: tuple[int, ...],
    chain: FineKChain,
    zero_local_system: Hashable,
    endpoints: dict[MaEndpointKey, PaintedBipartition],
) -> ChainPBPResult:
    """Attach every concrete Ma history of one VALUE chain to its endpoint."""

    histories_by_pbp: dict[PaintedBipartition, set[TwistHistory]] = {}
    for history in concrete_histories(chain):
        final_local_system = strict_ma_final_local_system(
            zero_local_system,
            chain,
            history,
        )
        subset_bits = theta_subset_bits(partition, chain, history)
        profile = row_multiplicity_profile(partition, subset_bits)
        key = final_local_system, profile
        try:
            pbp = endpoints[key]
        except KeyError as error:
            raise RuntimeError(
                "A strict Ma theta path is missing its painted-bipartition "
                f"endpoint for O^vee={partition}, history={history}"
            ) from error
        histories_by_pbp.setdefault(pbp, set()).add(history)
    return ChainPBPResult(chain, histories_by_pbp)


def concrete_history_text(
    history: TwistHistory,
    chain: FineKChain | None = None,
) -> str:
    """Format twists alone, or attach each twist to its orthogonal group."""

    if chain is None:
        return " -> ".join(history) if history else "direct Mp(0) lift"
    expected_length = len(chain.transitions) + 1
    if len(history) != expected_length:
        raise ValueError(
            "a concrete twist history must contain one label per "
            "orthogonal stage"
        )
    forms = [
        (transition.p, transition.q)
        for transition in chain.transitions
    ]
    forms.append((chain.final_p, chain.final_q))
    return " -> ".join(
        f"O({p},{q}) {twist}"
        for (p, q), twist in zip(forms, history)
    )


def print_report(
    partition: tuple[int, ...],
    chain_totals: tuple[int, ...],
    chains_by_label: dict[str, tuple[FineKChain, ...]],
    results_by_label: dict[str, tuple[ChainPBPResult, ...]],
    final_form: tuple[int, int] | None = None,
    show_diagrams: bool = False,
) -> None:
    n = sum(partition) // 2
    final_p, final_q = resolve_final_form(chain_totals[-1], final_form)
    tower_labels = ["Mp(0)"]
    for index, total in enumerate(chain_totals, start=1):
        tower_labels.append(f"O({total})" if index % 2 else f"Mp({total})")

    print(f"O^vee = {partition} in Sp({2 * n}, C)")
    print(f"theta2 totals: {chain_totals}")
    print("theta tower: " + " -> ".join(tower_labels))
    print(f"final real form: O({final_p},{final_q})")
    print("allowed raw PBP column shapes (P|Q), indexed by wp:")
    for wp, shape in sorted(orbit_pbp_shapes(partition).items()):
        pairs = tuple((2 * index + 2, 2 * index + 3) for index in wp)
        print(f"  wp={pairs or 'empty'}: {shape}")
    print()

    for k, label in enumerate(chains_by_label):
        results = results_by_label[label]
        concrete_paths = concrete_theta_paths(results)
        wedge_degrees = o_wedge_degrees_for_connected_type(final_p, k)
        print(
            f"concrete theta paths of coarse connected-K type {label}: "
            f"{len(concrete_paths)}"
        )
        print(f"  exact right-trivial O({final_p}) exterior degrees: {wedge_degrees}")
        for path_number, path in enumerate(concrete_paths, start=1):
            print(f"  path {path_number}:")
            print(
                "    concrete twist history: "
                + concrete_history_text(path.history, path.chain)
            )
            print(
                "      theta-lift chain: "
                + " -> ".join(path.chain.display_steps(path.history))
            )
            for pbp_number, pbp in enumerate(
                path.painted_bipartitions,
                start=1,
            ):
                print(f"    SO PBP hit {pbp_number}: {pbp.one_line()}")
                print(f"      raw PBP(O^vee, wp) extended DRC: {pbp.extended_drc!r}")
                if show_diagrams:
                    for line in str_dgms(pbp.extended_drc).splitlines():
                        print("      " + line)

        reached_parameters = {
            pbp
            for result in results
            for pbp in result.histories_by_pbp
        }
        so_groups: dict[
            tuple[Drc, tuple[int, ...]], list[PaintedBipartition]
        ] = {}
        for pbp in reached_parameters:
            so_groups.setdefault(pbp.so_parameter_key, []).append(pbp)
        print(f"  distinct SO painted bipartitions reached: {len(so_groups)}")
        for pbp_number, (_key, parameters) in enumerate(
            sorted(so_groups.items()),
            start=1,
        ):
            representative = parameters[0]
            reached = ",".join(
                "det" if value else "trivial"
                for value in sorted({item.outer_det_twist for item in parameters})
            )
            print(
                f"    [{pbp_number}] {representative.painted_one_line()}, "
                f"BMSZ_outer_epsilon_reached={reached}"
            )
        print()


def parse_orbit(text: str) -> tuple[int, ...]:
    try:
        rows = tuple(int(piece.strip()) for piece in text.split(",") if piece.strip())
    except ValueError as error:
        raise argparse.ArgumentTypeError("Use comma-separated integers, e.g. 6,4,2,2") from error
    try:
        return validate_dual_orbit(rows)
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def calculate(
    partition: Iterable[int],
    final_form: tuple[int, int] | None = None,
):
    """Return every strict Ma concrete history and its painted bipartition.

    Use :func:`concrete_theta_paths` to obtain the user-facing selected-twist
    path enumeration from any label's grouped results.
    """

    part = validate_dual_orbit(partition)
    totals = derive_chain_totals(part)
    resolved_final_form = resolve_final_form(totals[-1], final_form)
    candidate_chains_by_label = enumerate_fine_k_chains(
        totals,
        resolved_final_form,
    )
    ladder = load_pbp_ladder(part, totals)
    pbp_bijection, _ = build_pbp_bijection(part, "B")
    _zero_drc, zero_local_system = next(iter(ladder.drc_stages[0].items()))

    required_endpoint_keys: set[MaEndpointKey] = set()
    for candidate_chains in candidate_chains_by_label.values():
        for chain in candidate_chains:
            for history in concrete_histories(chain):
                required_endpoint_keys.add(
                    (
                        strict_ma_final_local_system(
                            zero_local_system,
                            chain,
                            history,
                        ),
                        row_multiplicity_profile(
                            part,
                            theta_subset_bits(part, chain, history),
                        ),
                    )
                )
    endpoints = build_ma_endpoint_map(
        part,
        ladder,
        pbp_bijection,
        resolved_final_form,
        required_endpoint_keys,
    )

    results_by_label = {}
    chains_by_label = {}
    for label, candidate_chains in candidate_chains_by_label.items():
        results = tuple(
            compute_strict_chain_pbps(
                part,
                chain,
                zero_local_system,
                endpoints,
            )
            for chain in candidate_chains
        )
        results_by_label[label] = results
        chains_by_label[label] = candidate_chains
    return part, totals, chains_by_label, results_by_label


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Enumerate theta2 fine-lowest-K-type chains and compute their "
            "painted bipartitions."
        )
    )
    parser.add_argument(
        "--orbit",
        type=parse_orbit,
        help=(
            "optional comma-separated orbit overriding DUAL_ORBIT in the file, "
            "for example 6,4,2,2"
        ),
    )
    parser.add_argument(
        "--show-diagrams",
        action="store_true",
        help="print the painted diagrams in addition to their exact tuple form",
    )
    args = parser.parse_args()

    selected_orbit = args.orbit if args.orbit is not None else DUAL_ORBIT
    part, totals, chains_by_label, results_by_label = calculate(selected_orbit)
    print_report(
        part,
        totals,
        chains_by_label,
        results_by_label,
        show_diagrams=args.show_diagrams,
    )


if __name__ == "__main__":
    main()
