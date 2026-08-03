"""Bounded audit of exact O-path and SO-normalized associated cycles.

For every concrete theta history, this test compares two independent
computations:

1. replay Ma's representation lift along that history; and
2. compute the associated cycle recursively from the resulting extended
   painted bipartition.

The strict path result first retains the reached ``O(p,q)`` extension.  It
must equal either the painted-bipartition reference cycle or, when the reached
outer extension is nontrivial, that cycle twisted by ``(1,1)``.  Only then is
the diagonal twist undone to check path independence after restriction to the
fixed ``SO(p,q)`` parameter.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from collections import defaultdict
from collections.abc import Iterable, Iterator


sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from combunipotent.LS import char_twist_B  # noqa: E402
from multiset import FrozenMultiset  # noqa: E402
from standalone import compute_AC, str_ils  # noqa: E402
from theta2_pbp import (  # noqa: E402
    CHAR_TWISTS,
    PaintedBipartition,
    calculate,
    load_pbp_ladder,
    strict_ma_final_local_system,
)


def partitions(
    total: int,
    maximum: int | None = None,
) -> Iterator[tuple[int, ...]]:
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
    """Generate every nonempty all-even partition up to ``max_total``."""

    for half_total in range(1, max_total // 2 + 1):
        for half_partition in partitions(half_total):
            yield tuple(2 * row for row in half_partition)


def canonical_ils(ils: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    """Discard trailing zero rows, which do not change a marked diagram."""

    entries = list(ils)
    while entries and entries[-1] == (0, 0):
        entries.pop()
    return tuple(entries)


def canonical_local_system(local_system: Iterable) -> FrozenMultiset:
    """Put a path-side local system in the comparison normalization."""

    return FrozenMultiset(canonical_ils(ils) for ils in local_system)


def associated_cycle_local_system(associated_cycle: Iterable) -> FrozenMultiset:
    """Expand associated-cycle coefficients into a comparable multiset."""

    terms = []
    for coefficient, ils in associated_cycle:
        terms.extend([canonical_ils(ils)] * int(coefficient))
    return FrozenMultiset(terms)


def marked_diagram_text(local_system: Iterable) -> str:
    """Format every irreducible component for an actionable failure report."""

    return " | ".join(str_ils(ils).replace("\n", "/") for ils in local_system)


def painted_bipartition_cycle(
    pbp: PaintedBipartition,
    cache: dict,
) -> FrozenMultiset:
    """Compute the SO associated cycle from one painted bipartition."""

    reference_drc = pbp.special_reference_extended_drc or pbp.extended_drc
    associated_cycle = compute_AC(
        reference_drc,
        frozenset(pbp.primitive_pair_indices),
        pbp.gamma,
        cache=cache,
    )
    return associated_cycle_local_system(associated_cycle)


def path_cycle(
    zero_local_system,
    chain,
    history,
) -> FrozenMultiset:
    """Replay one path while retaining its exact O-extension."""

    local_system = strict_ma_final_local_system(
        zero_local_system,
        chain,
        history,
    )
    return canonical_local_system(local_system)


def diagonal_twist_cycle(local_system: Iterable) -> FrozenMultiset:
    """Twist components individually so repeated multiplicities survive."""

    terms = []
    for ils in local_system:
        twisted = char_twist_B((ils,), CHAR_TWISTS["dd"])
        if len(twisted) != 1:
            raise AssertionError("diagonal twist changed the component count")
        terms.append(canonical_ils(next(iter(twisted))))
    return FrozenMultiset(terms)


def audit_orbit(
    partition: tuple[int, ...],
    associated_cycle_cache: dict,
) -> dict[str, int]:
    """Check every concrete path/PBP occurrence for one dual orbit."""

    part, totals, _chains_by_label, results_by_label = calculate(partition)
    ladder = load_pbp_ladder(part, totals)
    _zero_drc, zero_local_system = next(iter(ladder.drc_stages[0].items()))
    statistics = defaultdict(int)

    for label, results in results_by_label.items():
        statistics["fine_k_bins"] += 1
        paths_by_so_parameter = defaultdict(list)
        for result in results:
            for pbp, histories in result.histories_by_pbp.items():
                pbp_cycle = painted_bipartition_cycle(
                    pbp,
                    associated_cycle_cache,
                )
                for history in histories:
                    exact_path_cycle = path_cycle(
                        zero_local_system,
                        result.chain,
                        history,
                    )
                    expected_exact_cycle = (
                        diagonal_twist_cycle(pbp_cycle)
                        if pbp.outer_det_twist
                        else pbp_cycle
                    )
                    statistics["path_pbp_occurrences"] += 1
                    statistics[
                        "tensor_1_1_paths"
                        if pbp.outer_det_twist
                        else "same_paths"
                    ] += 1
                    if exact_path_cycle != expected_exact_cycle:
                        raise AssertionError(
                            "Exact O-path associated-cycle mismatch: "
                            f"O^vee={partition}, K-type={label}, "
                            f"history={history}, outer_det={pbp.outer_det_twist}, "
                            f"PBP={pbp.painted_one_line()}, "
                            f"exact path={marked_diagram_text(exact_path_cycle)}, "
                            "expected="
                            f"{marked_diagram_text(expected_exact_cycle)}"
                        )

                    # The determinant character is involutive.  Removing it
                    # compares both O-extensions at their common SO parameter.
                    computed_so_cycle = (
                        diagonal_twist_cycle(exact_path_cycle)
                        if pbp.outer_det_twist
                        else exact_path_cycle
                    )
                    paths_by_so_parameter[pbp.so_parameter_key].append(
                        computed_so_cycle
                    )
                    if computed_so_cycle != pbp_cycle:
                        raise AssertionError(
                            "SO-normalized path/PBP associated-cycle mismatch: "
                            f"O^vee={partition}, K-type={label}, "
                            f"history={history}, outer_det={pbp.outer_det_twist}, "
                            f"PBP={pbp.painted_one_line()}, "
                            f"path={marked_diagram_text(computed_so_cycle)}, "
                            f"PBP descent={marked_diagram_text(pbp_cycle)}"
                        )

        statistics["so_pbp_groups"] += len(paths_by_so_parameter)
        for so_parameter, cycles in paths_by_so_parameter.items():
            if len(cycles) > 1:
                statistics["multi_path_groups"] += 1
                statistics["multi_path_occurrences"] += len(cycles)
            if len(set(cycles)) != 1:
                raise AssertionError(
                    "Paths sharing one painted bipartition have different "
                    f"associated cycles: O^vee={partition}, K-type={label}, "
                    f"SO parameter={so_parameter!r}"
                )
    return dict(statistics)


def run_audit(max_total: int) -> dict[str, int | float]:
    """Audit every all-even orbit through the requested total size."""

    if max_total < 2 or max_total % 2:
        raise ValueError("--max-total must be a positive even integer")

    started = time.perf_counter()
    totals = defaultdict(int)
    associated_cycle_cache: dict = {}
    for partition in all_even_partitions(max_total):
        totals["partitions"] += 1
        for name, value in audit_orbit(
            partition,
            associated_cycle_cache,
        ).items():
            totals[name] += value

    totals["elapsed_seconds"] = time.perf_counter() - started
    assert totals["multi_path_groups"] > 0
    assert totals["same_paths"] > 0
    assert totals["tensor_1_1_paths"] > 0
    return dict(totals)


def main() -> None:
    """Run the command-line audit and print its exact coverage."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-total",
        type=int,
        default=16,
        help="largest total orbit size to audit (default: 16)",
    )
    arguments = parser.parse_args()
    statistics = run_audit(arguments.max_total)
    print(
        "associated-cycle path-independence audit passed: "
        f"{statistics['partitions']} partitions, "
        f"{statistics['fine_k_bins']} fine-K bins, "
        f"{statistics['so_pbp_groups']} SO painted-bipartition groups"
    )
    print(
        f"all {statistics['path_pbp_occurrences']} path/PBP occurrences agree: "
        f"{statistics['same_paths']} exact cycles are displayed as-is and "
        f"{statistics['tensor_1_1_paths']} require tensor (1,1); "
        f"{statistics['multi_path_groups']} multi-path groups contain "
        f"{statistics['multi_path_occurrences']} occurrences"
    )
    print(f"elapsed: {statistics['elapsed_seconds']:.1f}s")


if __name__ == "__main__":
    main()
