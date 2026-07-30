"""Focused and bounded regressions for :mod:`mp2_pbp`."""

from __future__ import annotations

import os
import sys
from fractions import Fraction
from itertools import product
from typing import Iterator


sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from combunipotent.drc import gp_form_M, verify_drc  # noqa: E402
from mp2_pbp import (  # noqa: E402
    MetaplecticCalculation,
    calculate_metaplectic,
    canonical_metaplectic_subset_bits,
    derive_metaplectic_chain_totals,
    load_metaplectic_pbp_ladder,
)
from standalone import (  # noqa: E402
    dpart2Wrepns_with_wp,
    dpart_to_bipartition,
    drc_shape,
)


def partitions(
    total: int,
    maximum: int | None = None,
) -> Iterator[tuple[int, ...]]:
    """Yield decreasing positive partitions of ``total``."""

    if total == 0:
        yield ()
        return
    upper = total if maximum is None else min(total, maximum)
    for first in range(upper, 0, -1):
        for tail in partitions(total - first, first):
            yield (first, *tail)


def all_even_partitions(max_total: int) -> Iterator[tuple[int, ...]]:
    """Yield all nonempty all-even partitions through ``max_total``."""

    for half_total in range(1, max_total // 2 + 1):
        for half_partition in partitions(half_total):
            yield tuple(2 * row for row in half_partition)


def subset_label(bits: tuple[int, ...]) -> str:
    indices = [str(index) for index, bit in enumerate(bits, start=1) if bit]
    return "{" + ", ".join(indices) + "}"


def subset_mask(bits: tuple[int, ...]) -> int:
    return sum(bit << index for index, bit in enumerate(bits))


def assert_complete_calculation(
    calculation: MetaplecticCalculation,
    *,
    check_cycles: bool,
) -> None:
    part = calculation.partition
    row_count = len(part)
    rank = sum(part) // 2
    expected_bits = set(product((0, 1), repeat=row_count))

    assert calculation.rank == rank
    assert calculation.totals == derive_metaplectic_chain_totals(part)
    assert calculation.tower[0] == "Mp(0)"
    assert calculation.tower[-1] == f"Mp({sum(part)})"
    assert len(calculation.paths) == 2 ** row_count
    assert {path.subset_bits for path in calculation.paths} == expected_bits
    assert [subset_mask(path.subset_bits) for path in calculation.paths] == list(
        range(2 ** row_count)
    )
    assert [path.number for path in calculation.paths] == list(
        range(1, 2 ** row_count + 1)
    )

    expected_shapes = dpart2Wrepns_with_wp(part, "M")
    special_shape = dpart_to_bipartition(part, "M")
    associated_cycle_cache = {}
    checked_cycles = set()
    for path in calculation.paths:
        assert path.subset_label == subset_label(path.subset_bits)
        assert path.name == f"Path {path.subset_label}"
        assert len(path.history) == (row_count + 1) // 2
        assert len(path.display_steps) >= 2 * len(path.history)
        assert path.display_steps[-1].startswith(f"Mp({sum(part)})")
        assert path.history_text.count("O(") == len(path.history)

        assert len(path.final_weight) == rank
        assert all(abs(value) == Fraction(1, 2) for value in path.final_weight)
        assert path.fine_degree == sum(
            value == Fraction(-1, 2) for value in path.final_weight
        )
        selected_sum = sum(
            (row // 2) * bit
            for row, bit in zip(reversed(part), path.subset_bits)
        )
        assert path.selected_half_row_sum == selected_sum == path.fine_degree

        pbp = path.painted_bipartition
        assert verify_drc(pbp.raw_drc, "M")
        assert verify_drc(pbp.special_reference_drc, "M")
        assert gp_form_M(pbp.raw_drc) == rank
        assert drc_shape(pbp.raw_drc) == expected_shapes[
            frozenset(pbp.primitive_pair_indices)
        ]
        assert drc_shape(pbp.special_reference_drc) == special_shape
        assert pbp.primitive_pairs == tuple(
            (2 * index + 1, 2 * index + 2)
            for index in pbp.primitive_pair_indices
        )

        if check_cycles and pbp.parameter_key not in checked_cycles:
            cycle = pbp.associated_cycle(cache=associated_cycle_cache)
            assert cycle
            assert all(coefficient > 0 for coefficient, _ils in cycle)
            checked_cycles.add(pbp.parameter_key)


def test_known_towers_and_input_validation() -> None:
    assert derive_metaplectic_chain_totals((2,)) == (1, 2)
    assert derive_metaplectic_chain_totals((4, 2)) == (3, 6)
    assert derive_metaplectic_chain_totals((6, 4, 2, 2)) == (3, 4, 9, 14)

    for invalid in ((), (3,), (2, 4), (2, 0), (2, True)):
        try:
            derive_metaplectic_chain_totals(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError(f"expected invalid orbit to fail: {invalid!r}")


def test_focused_paths_pbps_and_cycles() -> None:
    for part in ((2,), (4, 2), (6, 4, 2, 2)):
        assert_complete_calculation(
            calculate_metaplectic(part),
            check_cycles=True,
        )

    orbit_10 = calculate_metaplectic((10,))
    paths_10 = {path.subset_label: path for path in orbit_10.paths}
    assert paths_10["{}"].fine_degree == 0
    assert paths_10["{}"].final_weight == (Fraction(1, 2),) * 5
    assert paths_10["{1}"].fine_degree == 5
    assert paths_10["{1}"].final_weight == (Fraction(-1, 2),) * 5

    orbit_42 = calculate_metaplectic((4, 2))
    expected_cycles = {
        "{}": ((1, ((0, 0), (2, 1))),),
        "{1}": ((1, ((0, 0), (2, -1))),),
        "{2}": ((1, ((0, 0), (1, -2))),),
        "{1, 2}": ((1, ((0, 0), (-1, -2))),),
    }
    assert {
        path.subset_label: path.painted_bipartition.associated_cycle()
        for path in orbit_42.paths
    } == expected_cycles


def test_ladder_padding_and_endpoint_provenance() -> None:
    for part in ((2,), (4, 2), (2, 2, 2), (6, 4, 2, 2)):
        totals = derive_metaplectic_chain_totals(part)
        ladder = load_metaplectic_pbp_ladder(part, totals)
        assert len(ladder.drc_stages) == len(totals) + 1
        assert len(ladder.drc_stages[0]) == 1
        assert all(
            2 * gp_form_M(drc) == sum(part)
            for drc in ladder.drc_stages[-1]
        )

        reached_profiles = {
            bits
            for drc in ladder.drc_stages[-1]
            if (
                bits := canonical_metaplectic_subset_bits(part, ladder, drc)
            ) is not None
        }
        assert reached_profiles
        calculation = calculate_metaplectic(part)
        assert {
            path.subset_bits for path in calculation.paths
        } == set(product((0, 1), repeat=len(part)))


def test_all_even_partitions_through_total_14() -> None:
    checked = 0
    for part in all_even_partitions(14):
        assert_complete_calculation(
            calculate_metaplectic(part),
            check_cycles=False,
        )
        checked += 1
    assert checked == 44


def main() -> None:
    tests = (
        test_known_towers_and_input_validation,
        test_focused_paths_pbps_and_cycles,
        test_ladder_padding_and_endpoint_provenance,
        test_all_even_partitions_through_total_14,
    )
    for test in tests:
        test()
    print(f"test_mp2_pbp: {len(tests)} passed, 0 failed")


if __name__ == "__main__":
    main()
