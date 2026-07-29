"""Bounded audit of associated-cycle independence from the theta path.

For every concrete theta history, this test compares two independent
computations:

1. replay Ma's representation lift along that history; and
2. compute the associated cycle recursively from the resulting extended
   painted bipartition.

The path-side result is restricted from the reached ``O(p,q)`` extension to
the fixed ``SO(p,q)`` parameter before comparison.  Thus an outer determinant
twist is undone; it is not part of the painted bipartition.
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
    outer_det_twist: bool,
) -> FrozenMultiset:
    """Replay one path and restrict its reached O-extension to SO."""

    local_system = strict_ma_final_local_system(
        zero_local_system,
        chain,
        history,
    )
    if outer_det_twist:
        # The determinant character is involutive.  Removing it compares the
        # two O-extensions at their common SO painted-bipartition parameter.
        local_system = char_twist_B(local_system, CHAR_TWISTS["dd"])
    return canonical_local_system(local_system)


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
                    computed_path_cycle = path_cycle(
                        zero_local_system,
                        result.chain,
                        history,
                        pbp.outer_det_twist,
                    )
                    statistics["path_pbp_occurrences"] += 1
                    paths_by_so_parameter[pbp.so_parameter_key].append(
                        computed_path_cycle
                    )
                    if computed_path_cycle != pbp_cycle:
                        raise AssertionError(
                            "Path/PBP associated-cycle mismatch: "
                            f"O^vee={partition}, K-type={label}, "
                            f"history={history}, outer_det={pbp.outer_det_twist}, "
                            f"PBP={pbp.painted_one_line()}, "
                            f"path={marked_diagram_text(computed_path_cycle)}, "
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
        f"all {statistics['path_pbp_occurrences']} path/PBP occurrences agree; "
        f"{statistics['multi_path_groups']} multi-path groups contain "
        f"{statistics['multi_path_occurrences']} occurrences"
    )
    print(f"elapsed: {statistics['elapsed_seconds']:.1f}s")


if __name__ == "__main__":
    main()
