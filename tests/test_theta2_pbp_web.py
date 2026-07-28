"""Regression tests for the local theta2 painted-bipartition web API."""

from __future__ import annotations

import http.client
import json
import os
import sys
import threading


sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from theta2_pbp_web import (  # noqa: E402
    CALCULATION_CACHE,
    RequestError,
    create_server,
    parse_orbit_value,
    serialize_calculation,
)


def test_orbit_parser() -> None:
    assert parse_orbit_value([6, 4, 2, 2]) == (6, 4, 2, 2)
    assert parse_orbit_value("6, 4, 2, 2") == (6, 4, 2, 2)
    assert parse_orbit_value("64222") == (6, 4, 2, 2, 2)
    for invalid in ("", "643", [4, 6], [4, True], {"rows": [4, 2]}):
        try:
            parse_orbit_value(invalid)
        except RequestError:
            pass
        else:
            raise AssertionError(f"expected invalid orbit input to fail: {invalid!r}")


def test_explicit_paths_grouping_fixed_form_and_cache() -> None:
    CALCULATION_CACHE.clear()
    payload = serialize_calculation((4,))
    assert payload["orbit"] == [4]
    assert payload["totals"] == [5]
    assert payload["tower"] == ["Mp(0)", "O(5)"]
    assert payload["final_form"] == {"p": 2, "q": 3}

    k0 = payload["results"][0]
    assert k0["k"] == 0
    assert k0["candidate_concrete_path_count"] == 2
    assert k0["concrete_path_count"] == 2
    assert "candidate_count" not in k0
    assert "certified_count" not in k0
    assert k0["distinct_pbp_count"] == 2
    assert k0["o_wedge_degrees"] == [0, 2]
    assert "full_o_extension_count" not in k0
    assert {
        group["pbp"]["extended_drc"][1][0]
        for group in k0["groups"]
    } == {"sda", "srb"}
    assert {
        path["final"]["exact_label"]
        for group in k0["groups"]
        for path in group["paths"]
    } == {"(00|000)", "(11|000)"}
    assert {
        path["final"]["left_degree"]
        for group in k0["groups"]
        for path in group["paths"]
    } == {0, 2}
    assert all("all_o_extensions" not in group for group in k0["groups"])

    # The single orbit-keyed cache reuses the fixed O(n,n+1) calculation.
    assert serialize_calculation((4,)) is payload

    repeated = serialize_calculation((2, 2, 2, 2, 2, 2))
    k1 = repeated["results"][1]
    assert k1["distinct_pbp_count"] == 2
    assert k1["concrete_path_count"] == 12
    assert k1["o_wedge_degrees"] == [1, 5]
    assert {
        path["final"]["exact_label"]
        for group in k1["groups"]
        for path in group["paths"]
    } == {"(100000|0000000)", "(111110|0000000)"}
    assert all(
        [path["left_degree"] for path in group["paths"]]
        == sorted(path["left_degree"] for path in group["paths"])
        for group in k1["groups"]
    )
    assert all(
        path["target_label"] == path["final"]["exact_label"]
        for group in k1["groups"]
        for path in group["paths"]
    )

    raw_shape = serialize_calculation((6, 4, 2, 2, 2))
    allowed = {
        tuple(shape["primitive_pair_indices"]): (shape["p"], shape["q"])
        for shape in raw_shape["allowed_bipartition_shapes"]
    }
    assert allowed == {
        (): ([2, 1], [3, 1, 1]),
        (0,): ([1, 1], [3, 2, 1]),
    }

    middle = serialize_calculation((6, 4, 2))["results"][3]
    realizations = [
        realization
        for group in middle["groups"]
        for path in group["paths"]
        for realization in path["concrete_realizations"]
    ]
    assert {
        tuple(realization["twist_history"])
        for realization in realizations
    } == {("tt", "tt"), ("tt", "dt")}
    assert {
        realization["twist_history_text"]
        for realization in realizations
    } == {
        "O(2,1) tt -> O(6,7) tt",
        "O(2,1) tt -> O(6,7) dt",
    }
    for realization in realizations:
        final_twist = realization["twist_history"][-1]
        steps = realization["steps"]
        assert any(
            step.startswith(f"final twist[{final_twist}] -> O(6,7)")
            for step in steps
        )
        assert any(
            step.startswith("O(6,7)(111000|0000000)")
            and "untwisted lowest harmonic" in step
            for step in steps
        )
        assert not any("base L=" in step or "twist options[" in step for step in steps)

    boundary = serialize_calculation((4, 2, 2))
    assert [section["concrete_path_count"] for section in boundary["results"]] == [
        2,
        4,
        2,
    ]
    assert [
        section["candidate_concrete_path_count"]
        for section in boundary["results"]
    ] == [2, 4, 2]
    assert sum(
        section["concrete_path_count"] for section in boundary["results"]
    ) == 8
    assert sum(
        section["distinct_so_pbp_count"] for section in boundary["results"]
    ) == 6
    boundary_middle = boundary["results"][2]
    middle_paths = [
        path
        for group in boundary_middle["groups"]
        for path in group["paths"]
    ]
    assert len(middle_paths) == 2
    assert {path["number"] for path in middle_paths} == {1, 2}
    assert {
        tuple(path["concrete_realizations"][0]["twist_history"])
        for path in middle_paths
    } == {("tt", "tt"), ("tt", "dt")}
    assert all(len(path["concrete_realizations"]) == 1 for path in middle_paths)

    multistage = serialize_calculation((6, 4, 2, 2))
    first_multistage_realization = (
        multistage["results"][0]["groups"][0]["paths"][0]
        ["concrete_realizations"][0]
    )
    assert first_multistage_realization["twist_history_text"] == (
        "O(0,1) tt -> O(2,3) tt -> O(7,8) tt"
    )


def test_http_api() -> None:
    server = create_server("127.0.0.1", 0)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        connection = http.client.HTTPConnection(
            "127.0.0.1", server.server_address[1], timeout=30
        )
        body = json.dumps({"orbit": [4]}).encode("utf-8")
        connection.request(
            "POST",
            "/api/calculate",
            body=body,
            headers={"Content-Type": "application/json"},
        )
        response = connection.getresponse()
        document = json.loads(response.read())
        assert response.status == 200
        assert response.getheader("Content-Type") == "application/json; charset=utf-8"
        assert document["final_form"] == {"p": 2, "q": 3}
        assert document["results"][0]["distinct_pbp_count"] == 2

        local_origin = f"http://127.0.0.1:{server.server_address[1]}"
        connection.request(
            "POST",
            "/api/calculate",
            body=body,
            headers={
                "Content-Type": "application/json",
                "Origin": local_origin,
            },
        )
        response = connection.getresponse()
        document = json.loads(response.read())
        assert response.status == 200
        assert response.getheader("Access-Control-Allow-Origin") == local_origin
        assert document["final_form"] == {"p": 2, "q": 3}

        connection.request(
            "OPTIONS",
            "/api/calculate",
            headers={
                "Origin": "https://kdwong.github.io",
                "Access-Control-Request-Method": "POST",
                "Access-Control-Request-Headers": "content-type",
            },
        )
        response = connection.getresponse()
        response.read()
        assert response.status == 204
        assert response.getheader("Access-Control-Allow-Origin") == (
            "https://kdwong.github.io"
        )
        assert "POST" in response.getheader("Access-Control-Allow-Methods")

        connection.request(
            "POST",
            "/api/calculate",
            body=body,
            headers={
                "Content-Type": "application/json",
                "Origin": "https://example.com",
            },
        )
        response = connection.getresponse()
        response.read()
        assert response.status == 403

        bad_body = json.dumps({"orbit": [4], "final_form": [3, 2]}).encode("utf-8")
        connection.request(
            "POST",
            "/api/calculate",
            body=bad_body,
            headers={"Content-Type": "application/json"},
        )
        response = connection.getresponse()
        error = json.loads(response.read())
        assert response.status == 400
        assert "O(n,n+1)" in error["error"]

        connection.request("GET", "/health")
        response = connection.getresponse()
        assert response.status == 200
        assert json.loads(response.read()) == {"status": "ok"}

        connection.request("GET", "/")
        response = connection.getresponse()
        html = response.read().decode("utf-8")
        assert response.status == 200
        assert 'name="final-orientation"' not in html
        assert "SO(n,n+1)" in html
        assert "Jia-Jun Ma" in html
        assert "ems.press/journals/jems/articles/14298688" in html
        assert "doi.org/10.1090/jams/1082" in html
        assert 'id="empty-view"' not in html
        assert "<footer>" not in html

        connection.request("GET", "/app.js")
        response = connection.getresponse()
        javascript = response.read().decode("utf-8")
        assert response.status == 200
        assert 'body: JSON.stringify({ orbit })' in javascript
        assert "final-orientation" not in javascript
        assert "pathHistoryPreview(path.realizations)" in javascript
        assert 'makeStat(totalPaths, "concrete theta paths")' in javascript
        connection.close()
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def main() -> None:
    tests = (
        test_orbit_parser,
        test_explicit_paths_grouping_fixed_form_and_cache,
        test_http_api,
    )
    for test in tests:
        test()
    print(f"test_theta2_pbp_web: {len(tests)} passed, 0 failed")


if __name__ == "__main__":
    main()
