"""Large regression for O^vee=(8,6,4,4,2,2,2,2)."""

from __future__ import annotations

import os
import sys


sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from combunipotent.drc import reg_drc
from standalone import build_pbp_bijection, dpart2Wrepns_with_wp, drc_shape
from theta2_pbp import (
    PaintedBipartition,
    calculate,
    enumerate_fine_k_chains,
)


ORBIT = (8, 6, 4, 4, 2, 2, 2, 2)


def distinct_pbps(results):
    return {
        pbp
        for result in results
        for pbp in result.histories_by_pbp
    }


def main() -> None:
    final_form = (15, 16)
    _, totals, chains, results = calculate(ORBIT)
    candidates = enumerate_fine_k_chains(totals)

    assert totals == (1, 2, 5, 6, 9, 12, 17, 22, 31)
    assert {label: len(paths) for label, paths in chains.items()} == {
        "(0000000|00000000)": 2,
        "(1000000|00000000)": 6,
        "(1100000|00000000)": 8,
        "(1110000|00000000)": 20,
        "(1111000|00000000)": 14,
        "(1111100|00000000)": 32,
        "(1111110|00000000)": 34,
        "(1111111|00000000)": 32,
    }
    assert {
        label: len(distinct_pbps(label_results))
        for label, label_results in results.items()
    } == {
        "(0000000|00000000)": 2,
        "(1000000|00000000)": 2,
        "(1100000|00000000)": 4,
        "(1110000|00000000)": 6,
        "(1111000|00000000)": 10,
        "(1111100|00000000)": 10,
        "(1111110|00000000)": 12,
        "(1111111|00000000)": 14,
    }
    assert {
        label: len(paths)
        for label, paths in candidates.items()
    } == {
        label: len(paths)
        for label, paths in chains.items()
    }

    assert {
        PaintedBipartition(
            (("**", "*", "c"), ("**sda", "*sr", "sr", "d", "d")),
            False,
            (0, 1, 3),
        ),
        PaintedBipartition(
            (("**", "**", "*", "*"), ("**sda", "**d", "*", "*")),
            False,
            (0,),
        ),
        PaintedBipartition(
            (("***", "*", "c"), ("***da", "*d", "sd", "d", "d")),
            False,
            (1, 3),
        ),
        PaintedBipartition(
            (("***", "*c", "c"), ("***db", "*s", "r", "d", "d")),
            False,
            (3,),
        ),
        PaintedBipartition(
            (("***", "*c", "c", "c"), ("***sa", "*d", "d", "d")),
            True,
        ),
    } <= distinct_pbps(results["(1111000|00000000)"])
    assert len(distinct_pbps(results["(1111000|00000000)"])) == 10
    bijection, _ = build_pbp_bijection(ORBIT, "B")
    shapes = dpart2Wrepns_with_wp(ORBIT, "B")
    for label_results in results.values():
        for pbp in distinct_pbps(label_results):
            raw = reg_drc(pbp.plain_drc)
            assert raw in bijection
            _special, wp = bijection[raw]
            assert tuple(sorted(wp)) == pbp.primitive_pair_indices
            assert drc_shape(raw) == shapes[wp]
    print("test_theta2_pbp_86442222: 1 passed, 0 failed")


if __name__ == "__main__":
    main()
