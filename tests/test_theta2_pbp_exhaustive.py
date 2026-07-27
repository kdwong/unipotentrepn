"""Exhaustive structural audit for ``theta2_pbp.py``.

This is intentionally separate from the fast example-based regression suite.
It enumerates every all-even partition through a configurable total size,
checks the theta-chain/PBP invariants (including every wp-dependent raw shape
and primitive-pair subset), and reports (but does not assume) the
row-submultiset counting conjecture.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from collections import Counter
from collections.abc import Iterator


sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from combunipotent.drc import reg_drc, verify_drc
from combunipotent.drclift import gp_form_B_ext, split_extdrc_B
from standalone import (
    build_pbp_bijection,
    compute_PPidx,
    dpart2Wrepns_with_wp,
    drc_shape,
)
from theta2_pbp import (
    CHAR_TWISTS,
    calculate,
    determinant_twist_wedge,
    enumerate_fine_k_chains,
    final_left_wedge_values_label,
    o_character_to_mp,
    o_values_to_mp,
)


def partitions(total: int, maximum: int | None = None) -> Iterator[tuple[int, ...]]:
    """Generate decreasing positive integer partitions of ``total``."""

    if total == 0:
        yield ()
        return
    if maximum is None:
        maximum = total
    for first in range(min(total, maximum), 0, -1):
        for rest in partitions(total - first, first):
            yield (first, *rest)


def all_even_partitions(max_total: int) -> Iterator[tuple[int, ...]]:
    """Generate every nonempty all-even partition of size at most the bound."""

    for half_total in range(1, max_total // 2 + 1):
        for half_partition in partitions(half_total):
            yield tuple(2 * row for row in half_partition)


def row_submultiset_counts(
    partition: tuple[int, ...],
    max_degree: int,
) -> tuple[int, ...]:
    """Coefficient counts with equal rows indistinguishable.

    A row of size ``2s`` contributes degree ``s``.  If it occurs ``m`` times,
    its factor is ``1+x^s+...+x^(ms)``.
    """

    coefficients = [1] + [0] * max_degree
    for half_row, multiplicity in Counter(row // 2 for row in partition).items():
        updated = [0] * (max_degree + 1)
        for degree, coefficient in enumerate(coefficients):
            if coefficient == 0:
                continue
            for used in range(multiplicity + 1):
                target = degree + used * half_row
                if target > max_degree:
                    break
                updated[target] += coefficient
        coefficients = updated
    return tuple(coefficients)


def label_degree(label: str) -> int:
    left, _ = label[1:-1].split("|")
    return left.count("1")


def distinct_parameters(results) -> set:
    return {
        pbp
        for result in results
        for pbp in result.histories_by_pbp
    }


def check_final_chain(label: str, chain) -> tuple[bool, bool]:
    """Check final K-type/twist consistency.

    Returns ``(middle_alias, determinant_complement)`` for audit statistics.
    """

    degree = label_degree(label)
    exact_degree = chain.final_left_degree
    assert chain.final_label == label
    assert chain.final_left == (1,) * exact_degree
    assert chain.final_right == ()
    assert min(exact_degree, chain.final_p - exact_degree) == degree
    assert chain.final_exact_label == (
        f"{'('}{'1' * exact_degree}{'0' * (chain.final_p - exact_degree)}|"
        f"{'0' * chain.final_q})"
    )
    assert final_left_wedge_values_label(
        chain.final_p,
        chain.final_q,
        chain.final_left,
        chain.final_right,
    ) == label
    assert chain.final_twist_labels
    assert tuple(sorted(set(chain.final_twist_labels))) == chain.final_twist_labels

    for twist_label in chain.final_twist_labels:
        assert twist_label in CHAR_TWISTS
        det_left = twist_label[0] == "d"
        det_right = twist_label[1] == "d"
        left = (
            determinant_twist_wedge(chain.final_base_left, chain.final_p)
            if det_left
            else chain.final_base_left
        )
        right = (
            determinant_twist_wedge(chain.final_base_right, chain.final_q)
            if det_right
            else chain.final_base_right
        )
        assert (left, right) == (chain.final_left, chain.final_right)

    middle_alias = len(chain.final_twist_labels) > 1
    determinant_complement = exact_degree > chain.final_p // 2
    return middle_alias, determinant_complement


def check_intermediate_chain(chain, totals: tuple[int, ...]) -> None:
    """Replay every pre-final determinant twist and theta VALUE map."""

    assert len(chain.transitions) == len(totals) // 2
    for index, transition in enumerate(chain.transitions):
        assert transition.p + transition.q == totals[2 * index]
        assert transition.mp_total == totals[2 * index + 1]
        rank = transition.mp_total // 2

        if transition.initial_character:
            assert transition.base_left == transition.base_right == ()
            assert len(transition.twist_labels) == 1
            twist_label = transition.twist_labels[0]
            det_left = twist_label[0] == "d"
            det_right = twist_label[1] == "d"
            assert transition.twisted_left == (1,) * (
                transition.p if det_left else 0
            )
            assert transition.twisted_right == (1,) * (
                transition.q if det_right else 0
            )
            assert transition.mp_weight == o_character_to_mp(
                rank,
                transition.p,
                transition.q,
                det_left,
                det_right,
            )
            continue

        for twist_label in transition.twist_labels:
            det_left = twist_label[0] == "d"
            det_right = twist_label[1] == "d"
            left = (
                determinant_twist_wedge(
                    transition.base_left,
                    transition.p,
                )
                if det_left
                else transition.base_left
            )
            right = (
                determinant_twist_wedge(
                    transition.base_right,
                    transition.q,
                )
                if det_right
                else transition.base_right
            )
            assert (left, right) == (
                transition.twisted_left,
                transition.twisted_right,
            )
        assert transition.mp_weight == o_values_to_mp(
            rank,
            transition.p,
            transition.q,
            transition.twisted_left,
            transition.twisted_right,
        )


def audit_orbit(partition: tuple[int, ...]) -> dict[str, object]:
    n = sum(partition) // 2
    final_form = (n, n + 1)
    part, totals, chains_by_label, results_by_label = calculate(partition)
    assert part == partition
    assert totals[-1] == 2 * n + 1
    assert len(totals) % 2 == 1
    assert len(chains_by_label) == final_form[0] // 2 + 1
    assert tuple(chains_by_label) == tuple(results_by_label)
    candidate_chains_by_label = enumerate_fine_k_chains(totals)
    pbp_bijection, _ = build_pbp_bijection(partition, "B")
    expected_shapes = dpart2Wrepns_with_wp(partition, "B")
    valid_primitive_pair_indices = set(compute_PPidx(partition, "B"))

    row_counts = row_submultiset_counts(partition, final_form[0] // 2)
    mismatch_rows = []
    chain_count = 0
    candidate_chain_count = 0
    parameter_count = 0
    middle_alias_count = 0
    complement_count = 0
    outer_collisions = []

    for label, chains in chains_by_label.items():
        degree = label_degree(label)
        results = results_by_label[label]
        assert len(results) == len(chains)
        candidate_chains = candidate_chains_by_label[label]
        assert all(chain in candidate_chains for chain in chains)
        candidate_chain_count += len(candidate_chains)

        for chain, result in zip(chains, results):
            assert result.chain == chain
            assert (chain.final_p, chain.final_q) == final_form
            check_intermediate_chain(chain, totals)
            middle_alias, complement = check_final_chain(label, chain)
            middle_alias_count += int(middle_alias)
            complement_count += int(complement)

            for pbp, histories in result.histories_by_pbp.items():
                assert histories
                assert verify_drc(pbp.extended_drc, "B")
                assert gp_form_B_ext(pbp.extended_drc) == final_form
                assert split_extdrc_B(pbp.extended_drc)[1] in ("B+", "B-")
                assert set(pbp.primitive_pair_indices) <= (
                    valid_primitive_pair_indices
                )
                raw = reg_drc(pbp.plain_drc)
                assert raw in pbp_bijection
                _special, wp = pbp_bijection[raw]
                assert tuple(sorted(wp)) == pbp.primitive_pair_indices
                assert drc_shape(raw) == expected_shapes[wp]
                assert pbp.special_reference_extended_drc is not None
                assert verify_drc(pbp.special_reference_extended_drc, "B")
                assert gp_form_B_ext(pbp.special_reference_extended_drc) == final_form
                for history in histories:
                    assert len(history) == len(chain.transitions) + 1
                    for transition, twist_label in zip(chain.transitions, history):
                        assert twist_label in transition.twist_labels
                    assert history[-1] in chain.final_twist_labels

        parameters = distinct_parameters(results)
        so_parameters = {parameter.so_parameter_key for parameter in parameters}
        chain_count += len(chains)
        parameter_count += len(so_parameters)
        if len(parameters) != len(so_parameters):
            outer_collisions.append((label, len(so_parameters), len(parameters)))

        actual = len(so_parameters)
        predicted = row_counts[degree] * (
            1 if final_form[0] == 2 * degree else 2
        )
        if actual != predicted:
            mismatch_rows.append((degree, predicted, actual))

    return {
        "chains": chain_count,
        "candidate_chains": candidate_chain_count,
        "parameters": parameter_count,
        "middle_aliases": middle_alias_count,
        "complements": complement_count,
        "outer_collisions": tuple(outer_collisions),
        "mismatches": tuple(mismatch_rows),
    }


def run_audit(max_total: int) -> None:
    if max_total < 2 or max_total % 2:
        raise ValueError("--max-total must be a positive even integer")

    start = time.perf_counter()
    partition_count = 0
    bin_count = 0
    chain_count = 0
    candidate_chain_count = 0
    parameter_count = 0
    middle_alias_count = 0
    complement_count = 0
    mismatch_profile = {}
    outer_collisions = []

    for partition in all_even_partitions(max_total):
        summary = audit_orbit(partition)
        partition_count += 1
        n = sum(partition) // 2
        bin_count += n // 2 + 1
        chain_count += int(summary["chains"])
        candidate_chain_count += int(summary["candidate_chains"])
        parameter_count += int(summary["parameters"])
        middle_alias_count += int(summary["middle_aliases"])
        complement_count += int(summary["complements"])
        for degree, predicted, actual in summary["mismatches"]:
            mismatch_profile[(partition, degree)] = (predicted, actual)
        for collision in summary["outer_collisions"]:
            outer_collisions.append((partition, *collision))

    elapsed = time.perf_counter() - start
    print(
        "structural audit passed: "
        f"{partition_count} partitions, {bin_count} fine-K bins, "
        f"{chain_count}/{candidate_chain_count} certified/candidate chains, "
        f"{parameter_count} distinct PBP parameters"
    )
    print(
        "final determinant coverage: "
        f"{middle_alias_count} middle-degree alias chains, "
        f"{complement_count} over-middle complement chains"
    )
    print(f"PBP/outer-parameter collisions: {len(outer_collisions)}")
    print(
        "doubled row-submultiset heuristic: "
        f"{bin_count - len(mismatch_profile)}/{bin_count} bins agree"
    )
    for (partition, degree), (predicted, actual) in mismatch_profile.items():
        print(
            f"  {partition}, k={degree}: "
            f"row-submultisets={predicted}, PBPs={actual}"
        )
    print(f"elapsed: {elapsed:.1f}s")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-total",
        type=int,
        default=16,
        help="largest total orbit size to audit (default: 16)",
    )
    args = parser.parse_args()
    run_audit(args.max_total)


if __name__ == "__main__":
    main()
