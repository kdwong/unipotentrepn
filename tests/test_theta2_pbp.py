"""Focused regressions for explicit complementary theta paths and PBPs."""

from __future__ import annotations

import contextlib
import io
import os
import sys


sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from combunipotent.drc import reg_drc
from combunipotent.drclift import gp_form_B_ext
from standalone import build_pbp_bijection, dpart2Wrepns_with_wp, drc_shape
from theta2_pbp import (
    calculate,
    determinant_twist_wedge,
    derive_chain_totals,
    enumerate_fine_k_chains,
    load_pbp_ladder,
    o_wedge_degrees_for_connected_type,
    print_report,
    resolve_final_form,
    validate_dual_orbit,
)


def parameters(results):
    return {
        pbp
        for result in results
        for pbp in result.histories_by_pbp
    }


def so_keys(results):
    return {pbp.so_parameter_key for pbp in parameters(results)}


def assert_orbit_membership(orbit, results, final_form) -> None:
    bijection, _ = build_pbp_bijection(orbit, "B")
    shapes = dpart2Wrepns_with_wp(orbit, "B")
    for label_results in results.values():
        for pbp in parameters(label_results):
            raw = reg_drc(pbp.plain_drc)
            assert raw in bijection
            _special, wp = bijection[raw]
            assert tuple(sorted(wp)) == pbp.primitive_pair_indices
            assert drc_shape(raw) == shapes[wp]
            assert gp_form_B_ext(pbp.extended_drc) == final_form


def test_orbit_4_so23_has_two_actual_pbps() -> None:
    part, totals, chains, results = calculate((4,))
    assert part == (4,)
    assert totals == (5,)
    assert {label: len(paths) for label, paths in chains.items()} == {
        "(0|0)": 2,
        "(1|0)": 0,
    }

    k0_paths = chains["(0|0)"]
    assert [path.final_left_degree for path in k0_paths] == [0, 2]
    assert [path.final_exact_label for path in k0_paths] == [
        "(00|000)",
        "(11|000)",
    ]
    assert [path.final_twist_labels for path in k0_paths] == [
        ("tt",),
        ("dt",),
    ]
    assert so_keys(results["(0|0)"]) == {
        (((), ("sda",)), ()),
        (((), ("srb",)), ()),
    }
    reached = {
        (pbp.extended_drc, pbp.outer_det_twist)
        for pbp in parameters(results["(0|0)"])
    }
    assert reached == {
        (((), ("sda",)), False),
        (((), ("srb",)), True),
    }
    assert_orbit_membership(part, results, (2, 3))


def test_single_row_trivial_family_matches_complete_ma_packets() -> None:
    """Calibrate every direct Mp(0)->O(2n+1) lift through n=50.

    For Ovee=(2n), the only right-trivial fine type reached is k=0.  Its
    exact degrees 0 and p must recover every Ma PBP of the requested balanced
    real form ``O(n,n+1)``.  This tests actual tableaux, not a count formula.
    """

    for n in range(1, 51):
        orbit = (2 * n,)
        totals = derive_chain_totals(orbit)
        ladder = load_pbp_ladder(orbit, totals)
        final_stage = ladder.drc_stages[-1]

        final_form = (n, n + 1)
        _, calculated_totals, chains, results = calculate(orbit)
        assert calculated_totals == totals == (2 * n + 1,)
        k0_label = next(iter(chains))
        assert len(chains[k0_label]) == 2
        assert [path.final_left_degree for path in chains[k0_label]] == [0, n]
        assert [path.final_twist_labels for path in chains[k0_label]] == [
            ("tt",),
            ("dt",),
        ]
        assert all(
            not label_results
            for label, label_results in results.items()
            if label != k0_label
        )

        expected_ma_tableaux = {
            drc for drc in final_stage if gp_form_B_ext(drc) == final_form
        }
        actual_tableaux = {
            pbp.extended_drc for pbp in parameters(results[k0_label])
        }
        assert actual_tableaux == expected_ma_tableaux
        assert len(actual_tableaux) == 2
        assert all(
            not pbp.primitive_pair_indices
            for pbp in parameters(results[k0_label])
        )


def test_2x6_k1_has_ten_paths_and_two_pbps() -> None:
    orbit = (2, 2, 2, 2, 2, 2)
    _, totals, chains, results = calculate(orbit)
    assert totals == (1, 2, 5, 6, 9, 10, 13)

    k1_paths = chains["(100|000)"]
    assert len(k1_paths) == 10
    assert {path.final_left_degree for path in k1_paths} == {1, 5}
    assert {path.final_exact_label for path in k1_paths} == {
        "(100000|0000000)",
        "(111110|0000000)",
    }
    assert {path.final_twist_labels for path in k1_paths} == {("tt",), ("dt",)}
    assert so_keys(results["(100|000)"]) == {
        ((("*", "*", "c"), ("*b", "*", "d")), ()),
        ((("*", "c", "c"), ("*b", "d", "d")), ()),
    }
    assert len(so_keys(results["(100|000)"])) == 2
    assert_orbit_membership(orbit, results, (6, 7))


def test_6422_complements_are_computed_not_multiplied() -> None:
    orbit = (6, 4, 2, 2)
    _, totals, chains, results = calculate(orbit, final_form=(7, 8))
    assert totals == (1, 2, 5, 8, 15)
    assert {label: len(paths) for label, paths in chains.items()} == {
        "(000|0000)": 2,
        "(100|0000)": 2,
        "(110|0000)": 4,
        "(111|0000)": 4,
    }
    assert {
        label: len(so_keys(label_results))
        for label, label_results in results.items()
    } == {
        "(000|0000)": 2,
        "(100|0000)": 2,
        "(110|0000)": 4,
        "(111|0000)": 4,
    }
    assert {path.final_left_degree for path in chains["(100|0000)"]} == {1, 6}
    assert so_keys(results["(100|0000)"]) == {
        ((("*", "*"), ("*sda", "*d")), (0,)),
        ((("*", "*"), ("*srb", "*d")), (0,)),
    }
    assert_orbit_membership(orbit, results, (7, 8))


def test_64222_raw_wp_shape_and_middle_degree() -> None:
    orbit = (6, 4, 2, 2, 2)
    _, _, chains, results = calculate(orbit, final_form=(8, 9))
    assert {
        label: len(so_keys(label_results))
        for label, label_results in results.items()
    } == {
        "(0000|0000)": 2,
        "(1000|0000)": 2,
        "(1100|0000)": 4,
        "(1110|0000)": 6,
        "(1111|0000)": 2,
    }

    # The nontrivial-wp raw shape is retained as (1,1)|(3,2,1), i.e. the
    # user's row-notation (2 x 321), rather than replaced by the wp-empty one.
    k1_parameters = parameters(results["(1000|0000)"])
    assert len({pbp.so_parameter_key for pbp in k1_parameters}) == 2
    assert all(pbp.primitive_pairs == ((2, 3),) for pbp in k1_parameters)
    assert all(drc_shape(reg_drc(pbp.plain_drc)) == ((1, 1), (3, 2, 1)) for pbp in k1_parameters)

    middle_paths = chains["(1111|0000)"]
    assert {path.final_left_degree for path in middle_paths} == {4}
    assert all(path.final_twist_labels == ("dt", "tt") for path in middle_paths)
    assert len(so_keys(results["(1111|0000)"])) == 2
    assert_orbit_membership(orbit, results, (8, 9))


def test_concrete_twists_and_lowest_harmonics_are_printed() -> None:
    orbit, totals, chains, results = calculate((6, 4, 2))
    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        print_report(orbit, totals, chains, results)
    report = captured.getvalue()

    assert "concrete twist history: O(2,1) tt -> O(6,7) tt" in report
    assert "concrete twist history: O(2,1) tt -> O(6,7) dt" in report
    assert "final twist[tt] -> O(6,7)(111000|0000000)" in report
    assert "final twist[dt] -> O(6,7)(111000|0000000)" in report
    assert "O(6,7)(111000|0000000) [untwisted lowest harmonic]" in report
    assert "base L=" not in report
    assert "final twist options[dt,tt]" not in report


def test_cross_path_dedup_and_empty_bin() -> None:
    orbit = (4, 2, 2)
    _, totals, chains, results = calculate(orbit, final_form=(4, 5))
    assert totals == (3, 4, 9)
    assert len(chains["(10|00)"]) == 4
    assert len(so_keys(results["(10|00)"])) == 2
    assert_orbit_membership(orbit, results, (4, 5))

    _, totals_44, chains_44, results_44 = calculate((4, 4), final_form=(4, 5))
    assert totals_44 == (1, 4, 9)
    assert chains_44["(10|00)"] == ()
    assert results_44["(10|00)"] == ()


def test_final_form_and_input_validation() -> None:
    for partition in ((), (4, 6), (4, 3), (4, 0), (4, -2)):
        try:
            validate_dual_orbit(partition)
        except ValueError:
            pass
        else:
            raise AssertionError(f"Expected invalid orbit to fail: {partition}")

    assert resolve_final_form(5) == (2, 3)
    assert resolve_final_form(5, (2, 3)) == (2, 3)
    for bad_form in ((3, 2), (3, 3), (-1, 6), (5,), [2, 3]):
        try:
            resolve_final_form(5, bad_form)  # type: ignore[arg-type]
        except ValueError:
            pass
        else:
            raise AssertionError(f"Expected invalid final form: {bad_form!r}")

    assert derive_chain_totals((6, 4, 2)) == (3, 6, 13)
    assert o_wedge_degrees_for_connected_type(7, 1) == (1, 6)
    assert o_wedge_degrees_for_connected_type(8, 4) == (4,)
    assert [
        path.final_left_degree
        for path in enumerate_fine_k_chains((5,))["(0|0)"]
    ] == [0, 2]


def test_unsupported_determinant_twist_is_not_silently_accepted() -> None:
    try:
        determinant_twist_wedge((2,), 4)
    except ValueError:
        pass
    else:
        raise AssertionError("a non-exterior determinant twist must be rejected")


def main() -> None:
    tests = (
        test_orbit_4_so23_has_two_actual_pbps,
        test_single_row_trivial_family_matches_complete_ma_packets,
        test_2x6_k1_has_ten_paths_and_two_pbps,
        test_6422_complements_are_computed_not_multiplied,
        test_64222_raw_wp_shape_and_middle_degree,
        test_concrete_twists_and_lowest_harmonics_are_printed,
        test_cross_path_dedup_and_empty_bin,
        test_final_form_and_input_validation,
        test_unsupported_determinant_twist_is_not_silently_accepted,
    )
    for test in tests:
        test()
    print(f"test_theta2_pbp: {len(tests)} passed, 0 failed")


if __name__ == "__main__":
    main()
