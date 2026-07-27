#!/usr/bin/env python3
"""Local web interface for the painted-bipartition theta-chain calculator.

Run this file, then open the printed address in a browser.  The server uses
only Python's standard library; all mathematical work is delegated to
``theta2_pbp.calculate`` so the command-line and browser interfaces share the
same certified-chain computation.
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


from theta2_pbp import (  # noqa: E402 (the venv relaunch must happen first)
    ChainPBPResult,
    FineKChain,
    PaintedBipartition,
    calculate,
    concrete_history_text,
    enumerate_fine_k_chains,
    o_wedge_degrees_for_connected_type,
    orbit_pbp_shape,
    orbit_pbp_shapes,
    resolve_final_form,
    validate_dual_orbit,
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
    tuple[int, ...],
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
    return (
        pbp.extended_drc,
        pbp.primitive_pair_indices,
        pbp.outer_det_twist,
    )


def _serialize_group_path(
    result: ChainPBPResult,
    path_number: int,
    k: int,
    pbp: PaintedBipartition,
) -> dict[str, Any]:
    histories = sorted(result.histories_by_pbp[pbp])
    concrete_realizations = [
        {
            "twist_history": list(history),
            "twist_history_text": concrete_history_text(history, result.chain),
            "steps": _chain_steps(result.chain, history),
        }
        for history in histories
    ]
    preview_steps = (
        concrete_realizations[0]["steps"]
        if concrete_realizations
        else _chain_steps(result.chain)
    )
    serialized = {
        "id": path_number,
        "stable_id": f"k{k}-path-{path_number}",
        "name": f"Path {path_number}",
        "number": path_number,
        "target_label": result.chain.final_exact_label,
        "left_degree": result.chain.final_left_degree,
        "text": " -> ".join(preview_steps),
        # Retain one concrete chain here for older clients.  New clients use
        # concrete_realizations so a displayed final twist can never be a
        # grouped list of alternatives.
        "steps": preview_steps,
        "concrete_realizations": concrete_realizations,
        "twist_histories": [list(history) for history in histories],
        "twist_history_texts": [
            concrete_history_text(history, result.chain)
            for history in histories
        ],
        "outer_epsilon": pbp.outer_det_twist,
        "outer_epsilon_label": (
            "det" if pbp.outer_det_twist else "trivial"
        ),
    }
    serialized.update(_serialize_chain_structure(result.chain))
    return serialized


def _serialize_uncached(
    part: tuple[int, ...],
) -> dict[str, Any]:
    """Perform one calculation and build the complete browser data model."""

    final_form = resolve_final_form(sum(part) + 1)
    part, totals, _chains_by_label, results_by_label = calculate(part)
    candidate_chains = enumerate_fine_k_chains(totals)

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
        grouped: dict[
            tuple[Any, ...],
            list[tuple[int, ChainPBPResult, PaintedBipartition]],
        ] = {}
        for path_number, result in enumerate(results, start=1):
            for pbp in result.histories_by_pbp:
                grouped.setdefault(pbp.so_parameter_key, []).append(
                    (path_number, result, pbp)
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
            paths = [
                _serialize_group_path(result, path_number, k, parameter)
                for path_number, result, parameter in occurrences
            ]
            reached_extensions = sorted(
                {parameter.outer_det_twist for _, _, parameter in occurrences}
            )
            groups.append(
                {
                    "id": f"k{k}-pbp-{group_number}",
                    "number": group_number,
                    "pbp": _serialize_pbp(pbp),
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
                "candidate_count": len(candidate_chains[label]),
                "certified_count": len(results),
                "distinct_pbp_count": len(groups),
                "distinct_so_pbp_count": len(groups),
                "o_wedge_degrees": list(wedge_degrees),
                "o_wedge_degree_count": len(wedge_degrees),
                "middle_degree": final_p == 2 * k,
                "groups": groups,
            }
        )

    return {
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
        "results": labels,
    }


def serialize_calculation(
    partition: Iterable[int],
) -> dict[str, Any]:
    """Return a JSON-ready model for the fixed ``O(n,n+1)`` convention."""

    part = validate_dual_orbit(partition)
    # Keep the lookup and calculation under one lock.  In addition to guarding
    # upstream stdout redirection, this prevents two simultaneous requests for
    # a large orbit from both doing the same expensive work.
    with CALCULATION_LOCK:
        cached = CALCULATION_CACHE.get(part)
        if cached is not None:
            CALCULATION_CACHE.move_to_end(part)
            return cached
        payload = _serialize_uncached(part)
        CALCULATION_CACHE[part] = payload
        CALCULATION_CACHE.move_to_end(part)
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
) -> tuple[int, ...]:
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
    part = parse_orbit_value(document["orbit"])
    # Accept the now-obsolete field only when it repeats the fixed convention,
    # so an already-open copy of the old page continues to work after refresh.
    final_form_value = document.get("final_form")
    if final_form_value is None:
        return part
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
    return part


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
            orbit = parse_request_body(self.rfile.read(content_length))
            payload = serialize_calculation(orbit)
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
