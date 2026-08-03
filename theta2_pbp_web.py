#!/usr/bin/env python3
"""Web API for the SO(n,n+1) and Mp(2n,R) theta-path calculator.

Run this file, then open the printed address in a browser.  The server uses
only Python's standard library; all mathematical work is delegated to
the project calculators so the command-line and browser interfaces share the
same strict-Ma concrete-path computations.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import threading
import time
import webbrowser
from collections import OrderedDict
from fractions import Fraction
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import urlsplit
from urllib.request import urlopen


def _relaunch_in_project_venv() -> None:
    """Restart a direct invocation with the prepared project environment."""

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


from combunipotent.LS import char_twist_B  # noqa: E402
from theta2_pbp import (  # noqa: E402 (the venv relaunch must happen first)
    CHAR_TWISTS,
    ConcreteThetaPath,
    FineKChain,
    PaintedBipartition,
    calculate,
    concrete_history_text,
    concrete_theta_paths,
    load_pbp_ladder,
    o_wedge_degrees_for_connected_type,
    orbit_pbp_shape,
    orbit_pbp_shapes,
    resolve_final_form,
    strict_ma_final_local_system,
    theta_subset_bits,
    validate_dual_orbit,
)
from mp2_pbp import (  # noqa: E402
    MetaplecticConcretePath,
    MetaplecticFineChain,
    MetaplecticPaintedBipartition,
    calculate_metaplectic,
)
from standalone import (  # noqa: E402
    compute_AC,
    dpart2Wrepns_with_wp,
    dpart_to_bipartition,
    str_ils,
)


PROJECT_DIR = Path(__file__).resolve().parent
WEB_DIR = PROJECT_DIR / "web"
MAX_REQUEST_BYTES = 32 * 1024
DEFAULT_ALLOWED_ORIGINS = frozenset({"https://kdwong.github.io"})
ALLOWED_ORIGINS = DEFAULT_ALLOWED_ORIGINS | frozenset(
    origin.strip().rstrip("/")
    for origin in os.environ.get("THETA_PBP_ALLOWED_ORIGINS", "").split(",")
    if origin.strip()
)

# The calculator temporarily redirects stdout while loading upstream packet
# data.  Serializing calculations prevents concurrent requests from
# interfering with that process-global redirection.
CALCULATION_LOCK = threading.Lock()
CALCULATION_CACHE_SIZE = 8
CALCULATION_CACHE: OrderedDict[
    tuple[str, tuple[int, ...]],
    dict[str, Any],
] = OrderedDict()

STATIC_FILES = {
    "/": ("index.html", "text/html; charset=utf-8"),
    "/index.html": ("index.html", "text/html; charset=utf-8"),
    "/app.js": ("app.js", "text/javascript; charset=utf-8"),
    "/styles.css": ("styles.css", "text/css; charset=utf-8"),
}


class RequestError(ValueError):
    """A client error with a specific HTTP response status."""

    def __init__(self, message: str, status: HTTPStatus = HTTPStatus.BAD_REQUEST):
        super().__init__(message)
        self.status = status


def _fraction_text(value: Fraction) -> str:
    """Represent weights exactly; JSON floating point would lose fractions."""

    return str(value)


def _chain_steps(
    chain: FineKChain,
    history: tuple[str, ...] | None = None,
) -> list[str]:
    """Return display steps for one concrete twist realization."""

    return list(chain.display_steps(history))


def _serialize_chain_structure(chain: FineKChain) -> dict[str, Any]:
    """Expose structured path data as well as its ready-to-display steps."""

    transitions = []
    for transition in chain.transitions:
        transitions.append(
            {
                "orthogonal_form": {"p": transition.p, "q": transition.q},
                "initial_character": transition.initial_character,
                "base_left": list(transition.base_left),
                "base_right": list(transition.base_right),
                "twist_labels": list(transition.twist_labels),
                "twisted_left": list(transition.twisted_left),
                "twisted_right": list(transition.twisted_right),
                "mp_total": transition.mp_total,
                "mp_weight": [
                    _fraction_text(value) for value in transition.mp_weight
                ],
            }
        )
    return {
        "transitions": transitions,
        "final": {
            "orthogonal_form": {"p": chain.final_p, "q": chain.final_q},
            "label": chain.final_exact_label,
            "exact_label": chain.final_exact_label,
            "connected_label": chain.final_label,
            "left_degree": chain.final_left_degree,
            "base_left": list(chain.final_base_left),
            "base_right": list(chain.final_base_right),
            "twist_labels": list(chain.final_twist_labels),
            "left": list(chain.final_left),
            "right": list(chain.final_right),
        },
    }


def _serialize_pbp(pbp: PaintedBipartition) -> dict[str, Any]:
    """Serialize the actual raw ``tau_wp`` and its primitive-pair set."""

    p_columns, q_columns = pbp.plain_drc
    extended_p, extended_q = pbp.extended_drc
    special_reference = (
        pbp.special_reference_extended_drc or pbp.extended_drc
    )
    special_p, special_q = special_reference
    primitive_pairs = [list(pair) for pair in pbp.primitive_pairs]
    return {
        # Entries in these arrays are the diagram's column strings.  Keeping
        # them as separate strings lets the browser draw true tableau cells
        # without reverse-engineering the console's ASCII rendering.
        "p": list(p_columns),
        "q": list(q_columns),
        "gamma": pbp.gamma,
        "primitive_pair_indices": list(pbp.primitive_pair_indices),
        "primitive_pairs": primitive_pairs,
        "primitive_pairs_label": (
            "∅"
            if not primitive_pairs
            else "{" + ", ".join(
                f"({left},{right})" for left, right in pbp.primitive_pairs
            ) + "}"
        ),
        "extended_drc": [list(extended_p), list(extended_q)],
        "raw_extended_drc": [list(extended_p), list(extended_q)],
        "special_reference_extended_drc": [list(special_p), list(special_q)],
        "one_line": pbp.painted_one_line(),
    }


def _pbp_sort_key(pbp: PaintedBipartition) -> tuple[Any, ...]:
    """Return the deterministic display order for painted bipartitions."""

    return (
        pbp.extended_drc,
        pbp.primitive_pair_indices,
        pbp.outer_det_twist,
    )


def _subset_label(indices: tuple[int, ...]) -> str:
    """Format a subset of orbit-row indices with mathematical braces."""

    return "{" + ", ".join(map(str, indices)) + "}"


Ils = tuple[tuple[int, int], ...]
CycleCoefficients = dict[Ils, int]


def _canonical_ils(ils: Iterable[tuple[int, int]]) -> Ils:
    """Remove trailing zero rows, which do not change a marked diagram."""

    entries = [(int(p_value), int(q_value)) for p_value, q_value in ils]
    while entries and entries[-1] == (0, 0):
        entries.pop()
    return tuple(entries)


def _associated_cycle_coefficients(
    raw_terms: Iterable[tuple[int, Iterable[tuple[int, int]]]],
) -> CycleCoefficients:
    """Combine equal marked diagrams in first-occurrence order."""

    coefficients: CycleCoefficients = {}
    for coefficient, ils in raw_terms:
        canonical = _canonical_ils(ils)
        coefficients[canonical] = coefficients.get(canonical, 0) + int(
            coefficient
        )
        if coefficients[canonical] == 0:
            del coefficients[canonical]
    return coefficients


def _serialize_cycle_coefficients(
    coefficients: CycleCoefficients,
) -> dict[str, Any]:
    """Serialize consolidated marked associated-cycle components."""

    terms = []
    for ils, coefficient in coefficients.items():
        rendered = str_ils(ils)
        marked_rows = [] if rendered == "(trivial)" else rendered.splitlines()
        entries = [
            {"row_length": row_length, "p": p_value, "q": q_value}
            for row_length, (p_value, q_value) in enumerate(ils, start=1)
            if (p_value, q_value) != (0, 0)
        ]
        terms.append(
            {
                "coefficient": coefficient,
                "marked_rows": marked_rows,
                "underlying_partition": sorted(
                    (len(row) for row in marked_rows),
                    reverse=True,
                ),
                "ils_entries": entries,
            }
        )

    return {
        "term_count": len(terms),
        "total_multiplicity": sum(term["coefficient"] for term in terms),
        "terms": terms,
    }


def _diagonal_twist_cycle(
    coefficients: CycleCoefficients,
) -> CycleCoefficients:
    """Twist every type-B marked diagram by the global ``(1,1)`` character."""

    twisted_terms = []
    for ils, coefficient in coefficients.items():
        twisted_local_system = char_twist_B((ils,), CHAR_TWISTS["dd"])
        if len(twisted_local_system) != 1:
            raise RuntimeError("a diagonal twist changed the component count")
        twisted_ils = next(iter(twisted_local_system))
        twisted_terms.append((coefficient, twisted_ils))
    return _associated_cycle_coefficients(twisted_terms)


def _path_cycle_relation(
    path_coefficients: CycleCoefficients,
    reference_coefficients: CycleCoefficients,
    outer_det_twist: bool,
) -> str:
    """Verify and label the exact O-extension relative to its SO reference."""

    expected_coefficients = (
        _diagonal_twist_cycle(reference_coefficients)
        if outer_det_twist
        else reference_coefficients
    )
    if path_coefficients != expected_coefficients:
        expected_relation = (
            "its diagonal (1,1) twist"
            if outer_det_twist
            else "the untwisted reference"
        )
        raise RuntimeError(
            "an exact O-path associated cycle does not equal "
            f"{expected_relation} of its painted-bipartition cycle"
        )
    # Use the O-extension label even if a marked cycle happens to be fixed by
    # the diagonal twist, since equality alone cannot distinguish that case.
    return "tensor_1_1" if outer_det_twist else "same"


def _serialize_associated_cycle_data(
    reference_drc: tuple[tuple[str, ...], tuple[str, ...]],
    primitive_pair_indices: Iterable[int],
    parameter_type: str,
    cache: dict[Any, Any],
) -> dict[str, Any]:
    """Serialize the associated-cycle terms as marked Young diagrams.

    The associated-cycle descent is defined on the special reference DRC.
    The shifted diagram displayed by the painted-bipartition panel is not a
    valid substitute in this calculation.
    """

    raw_terms = compute_AC(
        reference_drc,
        frozenset(primitive_pair_indices),
        parameter_type,
        cache=cache,
    )

    # Iterated theta lifting can reach the same marked diagram through more
    # than one source term.  Combining those coefficients produces the actual
    # associated-cycle multiplicity while preserving first-occurrence order.
    return _serialize_cycle_coefficients(
        _associated_cycle_coefficients(raw_terms)
    )


def _painted_bipartition_cycle_coefficients(
    pbp: PaintedBipartition,
    cache: dict[Any, Any],
) -> CycleCoefficients:
    """Compute the SO-normalized cycle attached to one painted bipartition."""

    reference_drc = pbp.special_reference_extended_drc or pbp.extended_drc
    raw_terms = compute_AC(
        reference_drc,
        frozenset(pbp.primitive_pair_indices),
        pbp.gamma,
        cache=cache,
    )
    return _associated_cycle_coefficients(raw_terms)


def _serialize_group_path(
    part: tuple[int, ...],
    path: ConcreteThetaPath,
    path_number: int,
    k: int,
    pbp: PaintedBipartition,
    zero_local_system: Any,
    reference_cycle: CycleCoefficients,
) -> dict[str, Any]:
    """Serialize one path and its exact O-extension associated cycle."""

    history = path.history
    subset_bits = theta_subset_bits(part, path.chain, history)
    subset_indices = tuple(
        index for index, bit in enumerate(subset_bits, start=1) if bit
    )
    subset_label = _subset_label(subset_indices)
    subset_slug = "-".join(map(str, subset_indices)) or "empty"
    exact_local_system = strict_ma_final_local_system(
        zero_local_system,
        path.chain,
        history,
    )
    exact_cycle = _associated_cycle_coefficients(
        (1, ils) for ils in exact_local_system
    )
    cycle_relation = _path_cycle_relation(
        exact_cycle,
        reference_cycle,
        pbp.outer_det_twist,
    )
    concrete_realizations = [
        {
            "twist_history": list(history),
            "twist_history_text": concrete_history_text(history, path.chain),
            "steps": _chain_steps(path.chain, history),
            "exact_associated_cycle": _serialize_cycle_coefficients(
                exact_cycle
            ),
            "pbp_cycle_relation": cycle_relation,
            "pbp_cycle_twist": (
                [0, 0] if cycle_relation == "same" else [1, 1]
            ),
            "pbp_cycle_relation_label": (
                "same as painted-bipartition reference"
                if cycle_relation == "same"
                else "painted-bipartition reference twisted by (1,1)"
            ),
        }
    ]
    preview_steps = concrete_realizations[0]["steps"]
    serialized = {
        "id": path_number,
        "stable_id": f"k{k}-path-subset-{subset_slug}",
        "name": f"Path {subset_label}",
        "number": path_number,
        "subset_bits": list(subset_bits),
        "subset_indices": list(subset_indices),
        "subset_label": subset_label,
        "target_label": path.chain.final_exact_label,
        "left_degree": path.chain.final_left_degree,
        "text": " -> ".join(preview_steps),
        # Retain one concrete chain here for older clients.  New clients use
        # concrete_realizations so a displayed final twist can never be a
        # grouped list of alternatives.
        "steps": preview_steps,
        "concrete_realizations": concrete_realizations,
        "twist_histories": [list(history)],
        "twist_history_texts": [concrete_history_text(history, path.chain)],
        "outer_epsilon": pbp.outer_det_twist,
        "outer_epsilon_label": (
            "det" if pbp.outer_det_twist else "trivial"
        ),
    }
    serialized.update(_serialize_chain_structure(path.chain))
    return serialized


def _serialize_metaplectic_chain_structure(
    chain: MetaplecticFineChain,
) -> dict[str, Any]:
    """Expose the structured M-ending theta chain for one fine weight."""

    transitions = []
    for transition in chain.transitions:
        transitions.append(
            {
                "orthogonal_form": {"p": transition.p, "q": transition.q},
                "initial_character": transition.initial_character,
                "base_left": list(transition.base_left),
                "base_right": list(transition.base_right),
                "twist_labels": list(transition.twist_labels),
                "twisted_left": list(transition.twisted_left),
                "twisted_right": list(transition.twisted_right),
                "mp_total": transition.mp_total,
                "mp_weight": [
                    _fraction_text(value) for value in transition.mp_weight
                ],
            }
        )
    final_label = chain.display_steps()[-1].removesuffix(" [lowest harmonic]")
    return {
        "transitions": transitions,
        "final": {
            "group": f"Mp({chain.final_mp_total},R)",
            "mp_total": chain.final_mp_total,
            "mp_weight": [
                _fraction_text(value) for value in chain.final_mp_weight
            ],
            "label": final_label,
            "exact_label": final_label,
            "left_degree": chain.fine_degree,
        },
    }


def _serialize_metaplectic_pbp(
    pbp: MetaplecticPaintedBipartition,
) -> dict[str, Any]:
    """Serialize one actual type-M painted bipartition and its reference."""

    raw_p, raw_q = pbp.raw_drc
    special_p, special_q = pbp.special_reference_drc
    primitive_pairs = [list(pair) for pair in pbp.primitive_pairs]
    return {
        "p": list(raw_p),
        "q": list(raw_q),
        "parameter_type": "M",
        "primitive_pair_indices": list(pbp.primitive_pair_indices),
        "primitive_pairs": primitive_pairs,
        "primitive_pairs_label": (
            "∅"
            if not primitive_pairs
            else "{" + ", ".join(
                f"({left},{right})" for left, right in pbp.primitive_pairs
            ) + "}"
        ),
        "extended_drc": [list(raw_p), list(raw_q)],
        "raw_extended_drc": [list(raw_p), list(raw_q)],
        "special_reference_extended_drc": [
            list(special_p),
            list(special_q),
        ],
        "one_line": pbp.one_line(),
    }


def _serialize_metaplectic_path(
    path: MetaplecticConcretePath,
) -> dict[str, Any]:
    """Serialize one subset-labelled M-ending concrete twist history."""

    subset_slug = "-".join(map(str, path.subset_indices)) or "empty"
    steps = list(path.display_steps)
    realization = {
        "twist_history": list(path.history),
        "twist_history_text": path.history_text,
        "steps": steps,
    }
    serialized = {
        "id": path.number,
        "stable_id": f"mp-path-subset-{subset_slug}",
        "name": path.name,
        "number": path.number,
        "subset_bits": list(path.subset_bits),
        "subset_indices": list(path.subset_indices),
        "subset_label": path.subset_label,
        "selected_half_row_sum": path.selected_half_row_sum,
        "target_label": path.final_weight_text,
        "left_degree": path.fine_degree,
        "text": " -> ".join(steps),
        "steps": steps,
        "concrete_realizations": [realization],
        "twist_histories": [list(path.history)],
        "twist_history_texts": [path.history_text],
    }
    serialized.update(_serialize_metaplectic_chain_structure(path.chain))
    return serialized


def _serialize_so_uncached(
    part: tuple[int, ...],
) -> dict[str, Any]:
    """Perform one calculation and build the complete browser data model."""

    final_form = resolve_final_form(sum(part) + 1)
    part, totals, _chains_by_label, results_by_label = calculate(part)
    ladder = load_pbp_ladder(part, totals)
    _zero_drc, zero_local_system = next(iter(ladder.drc_stages[0].items()))
    associated_cycle_cache: dict[Any, Any] = {}

    n = sum(part) // 2
    empty_wp_p, empty_wp_q = orbit_pbp_shape(part)
    allowed_shapes = []
    for wp, (p_shape, q_shape) in sorted(orbit_pbp_shapes(part).items()):
        primitive_pairs = [
            [2 * index + 2, 2 * index + 3]
            for index in wp
        ]
        allowed_shapes.append(
            {
                "primitive_pair_indices": list(wp),
                "primitive_pairs": primitive_pairs,
                "primitive_pairs_label": (
                    "∅"
                    if not primitive_pairs
                    else "{" + ", ".join(
                        f"({left},{right})" for left, right in primitive_pairs
                    ) + "}"
                ),
                "p": list(p_shape),
                "q": list(q_shape),
            }
        )
    tower = ["Mp(0)"]
    for stage_number, total in enumerate(totals, start=1):
        tower.append(f"O({total})" if stage_number % 2 else f"Mp({total})")

    labels = []
    for k, (label, results) in enumerate(results_by_label.items()):
        concrete_paths = concrete_theta_paths(results)
        grouped: dict[
            tuple[Any, ...],
            list[tuple[int, ConcreteThetaPath, PaintedBipartition]],
        ] = {}
        for path_number, path in enumerate(concrete_paths, start=1):
            for pbp in path.painted_bipartitions:
                grouped.setdefault(pbp.so_parameter_key, []).append(
                    (path_number, path, pbp)
                )

        groups = []
        for group_number, so_key in enumerate(sorted(grouped), start=1):
            # Within a fixed connected K-type, display the target with fewer
            # 1's first: degree k precedes its complementary degree p-k.
            occurrences = sorted(
                grouped[so_key],
                key=lambda occurrence: (
                    occurrence[1].chain.final_left_degree,
                    occurrence[0],
                    _pbp_sort_key(occurrence[2]),
                ),
            )
            pbp = occurrences[0][2]
            reference_cycle = _painted_bipartition_cycle_coefficients(
                pbp,
                associated_cycle_cache,
            )
            paths = [
                _serialize_group_path(
                    part,
                    path,
                    path_number,
                    k,
                    parameter,
                    zero_local_system,
                    reference_cycle,
                )
                for path_number, path, parameter in occurrences
            ]
            reached_extensions = sorted(
                {parameter.outer_det_twist for _, _, parameter in occurrences}
            )
            groups.append(
                {
                    "id": f"k{k}-pbp-{group_number}",
                    "number": group_number,
                    "pbp": _serialize_pbp(pbp),
                    "associated_cycle": _serialize_cycle_coefficients(
                        reference_cycle
                    ),
                    "path_count": len({path_number for path_number, _, _ in occurrences}),
                    "path_occurrence_count": len(paths),
                    "outer_epsilons_reached": reached_extensions,
                    "outer_epsilon_labels": [
                        "det" if value else "trivial"
                        for value in reached_extensions
                    ],
                    "paths": paths,
                }
            )

        final_p, final_q = final_form
        wedge_degrees = o_wedge_degrees_for_connected_type(final_p, k)
        labels.append(
            {
                "k": k,
                "label": label,
                "concrete_path_count": len(concrete_paths),
                "distinct_pbp_count": len(groups),
                "distinct_so_pbp_count": len(groups),
                "o_wedge_degrees": list(wedge_degrees),
                "o_wedge_degree_count": len(wedge_degrees),
                "middle_degree": final_p == 2 * k,
                "groups": groups,
            }
        )

    return {
        "group": "so",
        "group_label": "SO(n,n+1)",
        "orbit": list(part),
        "orbit_text": "(" + ", ".join(map(str, part)) + ")",
        "ambient_group": f"Sp({2 * n}, C)",
        # Kept as a compatibility alias.  It is only the wp=empty reference,
        # not a required shape for nonempty wp.
        "expected_bipartition": {
            "p": list(empty_wp_p),
            "q": list(empty_wp_q),
        },
        "empty_wp_bipartition": {
            "p": list(empty_wp_p),
            "q": list(empty_wp_q),
        },
        "allowed_bipartition_shapes": allowed_shapes,
        "totals": list(totals),
        "tower": tower,
        "final_form": {"p": final_form[0], "q": final_form[1]},
        "final_group_label": f"O({final_form[0]},{final_form[1]})",
        "path_index_set_size": len(part),
        "expected_path_count": 2 ** len(part),
        "results": labels,
    }


def _serialize_metaplectic_uncached(
    part: tuple[int, ...],
) -> dict[str, Any]:
    """Perform one type-M calculation and build the browser data model."""

    calculation = calculate_metaplectic(part)
    associated_cycle_cache: dict[Any, Any] = {}
    expected_path_count = 2 ** len(part)
    if len(calculation.paths) != expected_path_count:
        raise RuntimeError(
            "the metaplectic calculation did not produce one path per row subset"
        )

    subset_bits = {path.subset_bits for path in calculation.paths}
    if len(subset_bits) != expected_path_count:
        raise RuntimeError("the metaplectic path labels do not form the full powerset")

    empty_wp_p, empty_wp_q = dpart_to_bipartition(part, "M")
    allowed_shapes = []
    allowed_shape_map = dpart2Wrepns_with_wp(part, "M")
    for wp, (p_shape, q_shape) in sorted(
        allowed_shape_map.items(),
        key=lambda item: tuple(sorted(item[0])),
    ):
        primitive_pair_indices = tuple(sorted(wp))
        primitive_pairs = [
            [2 * index + 1, 2 * index + 2]
            for index in primitive_pair_indices
        ]
        allowed_shapes.append(
            {
                "primitive_pair_indices": list(primitive_pair_indices),
                "primitive_pairs": primitive_pairs,
                "primitive_pairs_label": (
                    "∅"
                    if not primitive_pairs
                    else "{" + ", ".join(
                        f"({left},{right})" for left, right in primitive_pairs
                    ) + "}"
                ),
                "p": list(p_shape),
                "q": list(q_shape),
            }
        )

    paths_by_degree: dict[int, list[MetaplecticConcretePath]] = {}
    for path in calculation.paths:
        paths_by_degree.setdefault(path.fine_degree, []).append(path)

    results = []
    for degree in sorted(paths_by_degree):
        degree_paths = sorted(
            paths_by_degree[degree],
            key=lambda path: path.number,
        )
        grouped: dict[
            tuple[Any, ...],
            list[MetaplecticConcretePath],
        ] = {}
        for path in degree_paths:
            grouped.setdefault(
                path.painted_bipartition.parameter_key,
                [],
            ).append(path)

        groups = []
        for group_number, parameter_key in enumerate(sorted(grouped), start=1):
            group_paths = grouped[parameter_key]
            pbp = group_paths[0].painted_bipartition
            groups.append(
                {
                    "id": f"mp-k{degree}-pbp-{group_number}",
                    "number": group_number,
                    "pbp": _serialize_metaplectic_pbp(pbp),
                    "associated_cycle": _serialize_associated_cycle_data(
                        pbp.special_reference_drc,
                        pbp.primitive_pair_indices,
                        "M",
                        associated_cycle_cache,
                    ),
                    "path_count": len(group_paths),
                    "path_occurrence_count": len(group_paths),
                    "paths": [
                        _serialize_metaplectic_path(path)
                        for path in group_paths
                    ],
                }
            )

        final_weight_label = degree_paths[0].final_weight_text
        results.append(
            {
                "k": degree,
                "label": final_weight_label,
                "title": "Fine genuine U(n)-type",
                "degree_description": (
                    f"The final weight has {calculation.rank - degree} "
                    f"entries +1/2 and {degree} entries -1/2."
                ),
                "concrete_path_count": len(degree_paths),
                "distinct_pbp_count": len(groups),
                "groups": groups,
            }
        )

    return {
        "group": "mp",
        "group_label": "Mp(2n,R)",
        "orbit": list(part),
        "orbit_text": "(" + ", ".join(map(str, part)) + ")",
        "ambient_group": f"Sp({sum(part)}, C)",
        "expected_bipartition": {
            "p": list(empty_wp_p),
            "q": list(empty_wp_q),
        },
        "empty_wp_bipartition": {
            "p": list(empty_wp_p),
            "q": list(empty_wp_q),
        },
        "allowed_bipartition_shapes": allowed_shapes,
        "totals": list(calculation.totals),
        "tower": list(calculation.tower),
        "final_group_label": f"Mp({sum(part)},R)",
        "rank": calculation.rank,
        "path_index_set_size": len(part),
        "expected_path_count": expected_path_count,
        "results": results,
    }


def parse_group_value(value: Any) -> str:
    """Normalize the two public calculator choices."""

    if value is None:
        return "so"
    if not isinstance(value, str):
        raise RequestError("The 'group' field must be 'so' or 'mp'.")
    group = value.strip().casefold()
    aliases = {
        "so": "so",
        "so(n,n+1)": "so",
        "mp": "mp",
        "mp(2n,r)": "mp",
        "mp(2n,ℝ)": "mp",
    }
    try:
        return aliases[group]
    except KeyError as error:
        raise RequestError("The 'group' field must be 'so' or 'mp'.") from error


def serialize_calculation(
    partition: Iterable[int],
    group: str = "so",
) -> dict[str, Any]:
    """Return a JSON-ready model for the selected real-group convention."""

    part = validate_dual_orbit(partition)
    group_kind = parse_group_value(group)
    cache_key = (group_kind, part)
    # Keep the lookup and calculation under one lock.  In addition to guarding
    # upstream stdout redirection, this prevents two simultaneous requests for
    # a large orbit from both doing the same expensive work.
    with CALCULATION_LOCK:
        cached = CALCULATION_CACHE.get(cache_key)
        if cached is not None:
            CALCULATION_CACHE.move_to_end(cache_key)
            return cached
        payload = (
            _serialize_metaplectic_uncached(part)
            if group_kind == "mp"
            else _serialize_so_uncached(part)
        )
        CALCULATION_CACHE[cache_key] = payload
        CALCULATION_CACHE.move_to_end(cache_key)
        while len(CALCULATION_CACHE) > CALCULATION_CACHE_SIZE:
            CALCULATION_CACHE.popitem(last=False)
        return payload


def parse_orbit_value(value: Any) -> tuple[int, ...]:
    """Accept the UI's integer array and a convenient comma-separated form."""

    if isinstance(value, str):
        compact = value.strip()
        if compact.isdigit() and "," not in compact:
            # The user's customary notation ``64222`` means (6,4,2,2,2).
            # Multi-digit row sizes remain available in comma notation.
            pieces = list(compact)
        else:
            pieces = [piece.strip() for piece in value.split(",")]
        if not pieces or any(not piece for piece in pieces):
            raise RequestError("Enter rows as comma-separated even integers.")
        try:
            value = [int(piece) for piece in pieces]
        except ValueError as error:
            raise RequestError(
                "Enter rows as comma-separated even integers, for example 6,4,2,2."
            ) from error
    if not isinstance(value, list):
        raise RequestError("The 'orbit' field must be an array of even integers.")
    try:
        return validate_dual_orbit(value)
    except ValueError as error:
        raise RequestError(str(error)) from error


def parse_request_body(
    body: bytes,
) -> tuple[str, tuple[int, ...]]:
    try:
        document = json.loads(body.decode("utf-8"))
    except UnicodeDecodeError as error:
        raise RequestError("The request body must be UTF-8 JSON.") from error
    except json.JSONDecodeError as error:
        raise RequestError("The request body is not valid JSON.") from error
    if not isinstance(document, dict):
        raise RequestError("The JSON body must be an object with an 'orbit' field.")
    if "orbit" not in document:
        raise RequestError("The JSON body is missing the 'orbit' field.")
    group = parse_group_value(document.get("group"))
    part = parse_orbit_value(document["orbit"])
    # Accept the now-obsolete field only when it repeats the fixed convention,
    # so an already-open copy of the old page continues to work after refresh.
    final_form_value = document.get("final_form")
    if final_form_value is None:
        return group, part
    if group != "so":
        raise RequestError("The 'final_form' field applies only to SO(n,n+1).")
    if (
        not isinstance(final_form_value, list)
        or len(final_form_value) != 2
        or any(
            not isinstance(value, int) or isinstance(value, bool)
            for value in final_form_value
        )
    ):
        raise RequestError("The 'final_form' field must be [p,q].")
    final_form = tuple(final_form_value)
    try:
        resolve_final_form(sum(part) + 1, final_form)
    except ValueError as error:
        raise RequestError(str(error)) from error
    return group, part


class ThetaPBPRequestHandler(BaseHTTPRequestHandler):
    """Serve the allow-listed UI assets and the calculation endpoint."""

    server_version = "ThetaPBP/1.0"

    def do_OPTIONS(self) -> None:  # noqa: N802 (BaseHTTPRequestHandler API)
        if urlsplit(self.path).path != "/api/calculate":
            self._send_error_json(HTTPStatus.NOT_FOUND, "Not found")
            return
        if not self._origin_is_allowed():
            self._send_error_json(HTTPStatus.FORBIDDEN, "Origin not allowed")
            return
        self.send_response(HTTPStatus.NO_CONTENT)
        self._common_headers("text/plain; charset=utf-8", 0)
        self.send_header("Access-Control-Allow-Methods", "POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type, Accept")
        self.send_header("Access-Control-Max-Age", "86400")
        self.end_headers()

    def do_HEAD(self) -> None:  # noqa: N802 (BaseHTTPRequestHandler API)
        path = urlsplit(self.path).path
        if path in STATIC_FILES:
            self._serve_static(path, include_body=False)
            return
        if path == "/health":
            self._send_json(HTTPStatus.OK, {"status": "ok"}, include_body=False)
            return
        self._send_error_json(HTTPStatus.NOT_FOUND, "Not found", include_body=False)

    def do_GET(self) -> None:  # noqa: N802 (BaseHTTPRequestHandler API)
        path = urlsplit(self.path).path
        if path in STATIC_FILES:
            self._serve_static(path)
            return
        if path == "/health":
            self._send_json(HTTPStatus.OK, {"status": "ok"})
            return
        if path == "/favicon.ico":
            self.send_response(HTTPStatus.NO_CONTENT)
            self._common_headers("image/x-icon", 0)
            self.end_headers()
            return
        self._send_error_json(HTTPStatus.NOT_FOUND, "Not found")

    def do_POST(self) -> None:  # noqa: N802 (BaseHTTPRequestHandler API)
        if urlsplit(self.path).path != "/api/calculate":
            self._send_error_json(HTTPStatus.NOT_FOUND, "Not found")
            return
        if not self._origin_is_allowed():
            self._send_error_json(HTTPStatus.FORBIDDEN, "Origin not allowed")
            return

        try:
            content_length_text = self.headers.get("Content-Length")
            if content_length_text is None:
                raise RequestError(
                    "Content-Length is required.", HTTPStatus.LENGTH_REQUIRED
                )
            try:
                content_length = int(content_length_text)
            except ValueError as error:
                raise RequestError("Invalid Content-Length header.") from error
            if content_length < 0:
                raise RequestError("Invalid Content-Length header.")
            if content_length > MAX_REQUEST_BYTES:
                raise RequestError(
                    f"Request body exceeds {MAX_REQUEST_BYTES} bytes.",
                    HTTPStatus.REQUEST_ENTITY_TOO_LARGE,
                )
            group, orbit = parse_request_body(self.rfile.read(content_length))
            payload = serialize_calculation(orbit, group=group)
        except RequestError as error:
            self._send_error_json(error.status, str(error))
            return
        except Exception as error:  # keep implementation details out of JSON
            self.log_error("calculation failed: %r", error)
            self._send_error_json(
                HTTPStatus.INTERNAL_SERVER_ERROR,
                "The calculation failed. See the server window for details.",
            )
            return

        self._send_json(HTTPStatus.OK, payload)

    def _serve_static(self, request_path: str, include_body: bool = True) -> None:
        filename, content_type = STATIC_FILES[request_path]
        file_path = WEB_DIR / filename
        try:
            content = file_path.read_bytes()
        except OSError as error:
            self.log_error("could not read UI asset %s: %r", file_path, error)
            self._send_error_json(
                HTTPStatus.INTERNAL_SERVER_ERROR,
                f"The web interface asset {filename!r} is unavailable.",
                include_body=include_body,
            )
            return
        self.send_response(HTTPStatus.OK)
        self._common_headers(content_type, len(content))
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        if include_body:
            self.wfile.write(content)

    def _send_json(
        self,
        status: HTTPStatus,
        document: dict[str, Any],
        include_body: bool = True,
    ) -> None:
        content = json.dumps(
            document,
            ensure_ascii=False,
            separators=(",", ":"),
        ).encode("utf-8")
        self.send_response(status)
        self._common_headers("application/json; charset=utf-8", len(content))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        if include_body:
            self.wfile.write(content)

    def _send_error_json(
        self,
        status: HTTPStatus,
        message: str,
        include_body: bool = True,
    ) -> None:
        self._send_json(
            status,
            {"error": message, "status": int(status)},
            include_body=include_body,
        )

    def _common_headers(self, content_type: str, content_length: int) -> None:
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(content_length))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        origin = self.headers.get("Origin", "").rstrip("/")
        if origin and self._origin_is_allowed():
            self.send_header("Access-Control-Allow-Origin", origin)
            self.send_header("Vary", "Origin")

    def _origin_is_allowed(self) -> bool:
        """Allow same-origin/no-Origin clients and explicitly trusted webpages."""

        origin = self.headers.get("Origin")
        if origin is None:
            return True
        normalized_origin = origin.rstrip("/")
        if normalized_origin in ALLOWED_ORIGINS:
            return True

        parsed_origin = urlsplit(normalized_origin)
        request_host = self.headers.get("Host", "").casefold()
        return (
            parsed_origin.scheme in {"http", "https"}
            and parsed_origin.netloc.casefold() == request_host
            and parsed_origin.path in {"", "/"}
            and not parsed_origin.query
            and not parsed_origin.fragment
            and parsed_origin.username is None
            and parsed_origin.password is None
        )


class ThetaPBPServer(ThreadingHTTPServer):
    daemon_threads = True
    allow_reuse_address = True


def create_server(host: str = "127.0.0.1", port: int = 8000) -> ThetaPBPServer:
    return ThetaPBPServer((host, port), ThetaPBPRequestHandler)


def _open_browser_when_ready(url: str) -> None:
    """Open the page only after the server has answered a health request."""

    health_url = url + "health"
    for _attempt in range(50):
        try:
            with urlopen(health_url, timeout=1) as response:
                if response.status == HTTPStatus.OK:
                    webbrowser.open(url, new=2)
                    return
        except OSError:
            pass
        time.sleep(0.1)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the local painted-bipartition theta-chain webpage."
    )
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="interface to listen on (default: 127.0.0.1, local computer only)",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8000,
        help="TCP port to listen on (default: 8000; use 0 for an available port)",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="do not open the webpage automatically after the server starts",
    )
    args = parser.parse_args()
    if not 0 <= args.port <= 65535:
        parser.error("--port must be between 0 and 65535")

    try:
        server = create_server(args.host, args.port)
    except OSError as error:
        parser.error(f"could not start the server: {error}")

    actual_host, actual_port = server.server_address[:2]
    display_host = "127.0.0.1" if actual_host in ("0.0.0.0", "::") else actual_host
    url = f"http://{display_host}:{actual_port}/"
    print(f"Painted-bipartition webpage: {url}", flush=True)
    print("Press Ctrl+C to stop the server.", flush=True)
    if not args.no_browser:
        # The helper waits for a successful health response, avoiding a race
        # between the browser's first request and serve_forever().
        threading.Thread(
            target=_open_browser_when_ready,
            args=(url,),
            daemon=True,
        ).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nServer stopped.", flush=True)
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
