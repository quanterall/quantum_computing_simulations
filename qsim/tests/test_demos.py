"""The ten demos: their exact outcome probabilities, the response fields the demo page reads, and the API."""
import numpy as np
import pytest
from fastapi.testclient import TestClient

from qsim import demos
from qsim.app import app

SHOTS = 20_000


def share(counts, key):
    return counts.get(key, 0) / sum(counts.values())


def test_bell():
    r = demos.bell(SHOTS)
    assert set(r["counts"]) <= {"00", "11"}
    assert abs(share(r["counts"], "00") - 0.5) < 0.02


def test_teleport_random_state_arrives():
    # a random ψ teleported and undone on Bob's side: q2 (the leftmost character) reads 0 every time, Alice's two
    # bits are uniformly random
    for _ in range(3):
        r = demos.teleport(SHOTS)
        assert set(r["counts"]) == {"000", "001", "010", "011"} and r["bob_zero_fraction"] == 1.0
        assert all(abs(share(r["counts"], k) - 0.25) < 0.02 for k in r["counts"])


def test_grover_finds_11():
    assert demos.grover(1000)["counts"] == {"11": 1000}


def test_qml_matches_its_exact_distribution():
    from qsim.engine import Circuit
    p = Circuit(4).rx(0, 0.5).rx(1, 0.3).cx(0, 2).cx(1, 3).rx(2, 0.7).rx(3, 0.7).cx(2, 3).probabilities()
    r = demos.qml(SHOTS)["counts"]
    for i, pi in enumerate(p):
        assert abs(share(r, format(i, "04b")) - pi) < 0.02


def test_shor15_finds_the_period_and_the_factors():
    r = demos.shor15(SHOTS)
    assert set(r["counts"]) == {"000", "010", "100", "110"}                     # k/8 with period 4
    assert all(abs(share(r["counts"], k) - 0.25) < 0.02 for k in r["counts"])
    by_k = {c["measured"]: c for c in r["classical"]}
    assert by_k["000"]["factors"] == []
    assert all(by_k[k]["period"] == 4 and by_k[k]["factors"] == [3, 5] for k in ("010", "100", "110"))


def test_cswap_and_cp():
    from qsim.engine import Circuit
    assert Circuit(3).x(0).x(1).cswap(0, 1, 2).sample(10) == {"101": 10}       # control 1: q1, q2 swapped
    assert Circuit(3).x(1).cswap(0, 1, 2).sample(10) == {"010": 10}            # control 0: nothing
    c = Circuit(2).h(0).h(1).cp(0, 1, np.pi).h(0).h(1)                          # cp(pi) = cz: H.CZ.H on |00>
    assert np.allclose(c.probabilities(), Circuit(2).h(0).h(1).cz(0, 1).h(0).h(1).probabilities())


@pytest.mark.parametrize("name", ["bitflip", "phaseflip", "shor9", "steane7", "surface3"])
def test_qec_response_fields(name):
    r = demos.QEC[name](2000, "none")
    for k in ("counts_before", "syndrome_counts", "corrected_counts", "logical_fidelity", "error_type",
              "error_qubit", "plots", "circuit_qasm", "execution_time_ms", "qubits", "shots"):
        assert k in r
    assert set(r["plots"]) == {"counts_before", "corrected_counts", "syndrome"}
    assert sum(r["corrected_counts"].values()) == 2000


def test_bitflip_no_error_and_error_on_q1():
    r = demos.bitflip(1000, "none")
    assert r["counts_before"] == {"000": 1000} and r["syndrome_counts"] == {"00000": 1000}
    assert r["logical_fidelity"] == 1.0
    r = demos.bitflip(1000, "1")
    assert r["counts_before"] == {"010": 1000}
    assert r["syndrome_counts"] == {"11010": 1000}          # both checks fire: s0 (q3) and s1 (q4)
    assert r["logical_fidelity"] == 1.0                     # one flipped bit of three: the vote corrects it


def test_error_qubit_out_of_range_means_none():
    assert demos.bitflip(10, "7")["counts_before"] == {"000": 10}
    assert demos.bitflip(10, "x")["counts_before"] == {"000": 10}


def test_api_paths_shots_and_cors():
    client = TestClient(app)
    r = client.get("/quantum/bell?shots=50")
    assert r.status_code == 200 and r.json()["shots"] == 50 and r.json()["example"] == "bell"
    assert r.headers["access-control-allow-origin"] == "*"
    assert client.get("/quantum/bell?shots=0").json()["shots"] == 1024
    assert client.get("/quantum/bell?shots=99999").json()["shots"] == 1024
    assert client.get("/quantum/qec/bitflip?shots=10&error_qubit=2").json()["counts_before"] == {"100": 10}
    assert client.get("/quantum/nope").status_code == 404
    assert client.get("/quantum/qec/nope").status_code == 404
    assert client.get("/health").json()["status"] == "ok"


def test_all_ten_on_the_api():
    client = TestClient(app)
    for name in demos.AI:
        assert client.get(f"/quantum/{name}?shots=100").json()["status"] == "ok"
    for name in demos.QEC:
        assert client.get(f"/quantum/qec/{name}?shots=100&error_qubit=1").json()["status"] == "ok"


@pytest.mark.parametrize("name,k", [("phaseflip", 3), ("shor9", 9), ("steane7", 7)])
def test_every_error_position_is_corrected(name, k):
    for q in ["none"] + [str(i) for i in range(k)]:
        r = demos.QEC[name](200, q)
        assert r["logical_fidelity"] == 1.0, (name, q, r["syndrome"], r["correction"])
        if q != "none" and name != "phaseflip":
            assert r["logical_before"] == {"1": 200}, (name, q)


def test_surface3_corrects_every_single_error():
    for q in ["none"] + [str(i) for i in range(9)]:
        r = demos.QEC["surface3"](200, q)
        assert r["logical_fidelity"] == 1.0, (q, r["syndrome"], r["correction"])
        if q in ("0", "1", "2"):
            assert r["logical_before"] == {"1": 200}
