#!/usr/bin/env bash
# test_all_quantum_demos.sh — Quanterall Quantum Demo Suite
# Tests all the endpoints and prints pass/fail for each. On the host: bash test_all_quantum_demos.sh http://127.0.0.1:8100

BASE="${1:-http://localhost:8080}"
PASS=0; FAIL=0

run_test() {
  local name="$1"
  local url="$2"
  local response
  sleep 0.25; response=$(curl -sf "$url" 2>/dev/null)   # under nginx's 5 requests a second
  if echo "$response" | grep -q '"status":"ok"'; then
    echo "  ✓ $name"
    ((PASS++))
  else
    echo "  ✗ $name  →  $(echo "$response" | head -c 120)"
    ((FAIL++))
  fi
}

echo ""
echo "Quanterall Quantum Demo Suite — Integration Test"
echo "Target: $BASE"
echo "────────────────────────────────────────────────"
# a correction cycle must leave the logical qubit at 0 in every shot
run_fixed() {
  local name="$1" url="$2" response
  sleep 0.25; response=$(curl -sf "$url" 2>/dev/null)   # under nginx's 5 requests a second
  if echo "$response" | grep -q '"logical_fidelity":1.0'; then
    echo "  ✓ $name"; ((PASS++))
  else
    echo "  ✗ $name  →  $(echo "$response" | grep -o '"logical_fidelity":[0-9.]*')"; ((FAIL++))
  fi
}

echo ""
echo "Algorithms, protocols and QML:"
run_test "Bell State"           "$BASE/quantum/bell?shots=100"
run_test "Teleportation"        "$BASE/quantum/teleport?shots=100"
run_test "Grover Search"        "$BASE/quantum/grover?shots=100"
run_test "QML Classifier"       "$BASE/quantum/qml?shots=100"
run_test "Shor-15"              "$BASE/quantum/shor15?shots=100"

# a check on the JSON answer, in Python: the expression gets d, the parsed response
run_check() {
  local name="$1" url="$2" expr="$3" response
  sleep 0.25; response=$(curl -sf "$url" 2>/dev/null)   # under nginx's 5 requests a second
  if echo "$response" | python3 -c "import sys,json; d=json.load(sys.stdin); sys.exit(0 if ($expr) else 1)" 2>/dev/null; then
    echo "  ✓ $name"; ((PASS++))
  else
    echo "  ✗ $name"; ((FAIL++))
  fi
}
run_check "Shor-15 finds 3 x 5"     "$BASE/quantum/shor15?shots=200" \
  'set(d["counts"]) <= {"000","010","100","110"} and any(c["factors"] == [3, 5] for c in d["classical"])'
run_check "Teleport random ψ: Bob 0" "$BASE/quantum/teleport?shots=300" 'all(k[0] == "0" for k in d["counts"]) and d["bob_zero_fraction"] == 1.0 and len(d["counts"]) == 4'
run_check "Shor general: 21 with a=2" "$BASE/quantum/shor?N=21&a=2&shots=2000&seed=4" 'd["period_true"] == 6 and d["success"] > 0.7 and any(c["factors"] == [3, 7] for c in d["classical"])'
run_check "Shor general: 15 with a=14 fails" "$BASE/quantum/shor?N=15&a=14&shots=500" 'd["success"] == 0'
run_check "Shor general: gcd shortcut"  "$BASE/quantum/shor?N=21&a=6" 'd["circuit_run"] is False and d["factors"] == [3, 7]'
run_check "Shor general: 21, a=4 odd r" "$BASE/quantum/shor?N=21&a=4&shots=500&seed=2" 'd["success"] == 0 and "period odd" in d["failures"] and d["success_exact"] == 0'
run_check "Shor: success fractions"   "$BASE/quantum/shor15?shots=4000&seed=5" '0.45 < d["success_direct"] < 0.55 and 0.70 < d["success_with_check"] < 0.80'
run_check "Seed repeats a run"      "$BASE/quantum/qml?shots=300&seed=12345" 'd["seed"] == 12345'
A=$(curl -sf "$BASE/quantum/qml?shots=300&seed=12345" | python3 -c "import sys,json; print(sorted(json.load(sys.stdin)['counts'].items()))")
B=$(curl -sf "$BASE/quantum/qml?shots=300&seed=12345" | python3 -c "import sys,json; print(sorted(json.load(sys.stdin)['counts'].items()))")
if [ -n "$A" ] && [ "$A" = "$B" ]; then echo "  ✓ Same seed, same counts"; ((PASS++)); else echo "  ✗ Same seed, same counts"; ((FAIL++)); fi

run_check "Grover n=8 at its optimum"  "$BASE/quantum/grover/search?n=8&marked=77&shots=300" \
  'd["optimal_iterations"] == 12 and d["success"] > 0.97 and abs([c for c in d["curve"] if c["k"] == 12][0]["p"] - [c for c in d["curve"] if c["k"] == 12][0]["theory"]) < 1e-6'
run_check "Surface-3 corrects Y on q4" "$BASE/quantum/qec/surface3?shots=100&error_qubit=4" 'd["logical_fidelity"] == 1.0 and d["correction"] == ["x q[4]", "z q[4]"]'
run_check "Stabilizer: GHZ on 500"     "$BASE/quantum/stabilizer/ghz?n=500&shots=500" 'd["counts"]["other"] == 0 and d["counts"]["all 0"] > 150 and d["counts"]["all 1"] > 150'
run_check "Stabilizer: noiseless memory" "$BASE/quantum/stabilizer/memory?code=surface&d=5&p=0&shots=500" 'd["shots_with_detections"] == 0 and d["observable_flip_rate_undecoded"] == 0'
run_check "Stabilizer: noisy memory"   "$BASE/quantum/stabilizer/memory?code=surface&d=3&p=0.01&shots=2000" '0.02 < d["detection_rate"] < 0.2'

run_check "Optics: HOM at gt=1"      "$BASE/quantum/optics/hom?gt=1" 'abs(d["point"]["coincidence"] - d["point"]["formula"]) < 1e-12'
run_check "Optics: NS 1/4, fidelity 1" "$BASE/quantum/optics/ns?seed=4" 'abs(d["success"] - 0.25) < 1e-6 and d["fidelity"] > 1 - 1e-12'
run_check "Optics: CZ exact"          "$BASE/quantum/optics/cz?input=11&shots=2000&seed=3" 'd["exact"]["process_fidelity"] > 1 - 1e-12 and abs(d["exact"]["success"] - 1/16) < 1e-6 and set(d["heralded_outputs"]) == {"11"}'
run_check "Optics: type-I fusion"     "$BASE/quantum/optics/fusion?shots=500&seed=2" 'abs(d["success_exact"] - d["hh_plus_vv_weight"]) < 1e-12'
run_check "Optics: CZ power check"    "$BASE/quantum/optics/cz?input=11&shots=100&theta2_offset=0.5" 'd["exact"]["process_fidelity"] < 0.9999'
run_check "Optics: GKP Z readout"     "$BASE/quantum/optics/gkp?db=10&sigma=0.4&shots=10000&seed=3" 'd["shots"] == 10000 and abs(d["error_sampled"] - d["error_exact"]) < 0.01 and abs(d["error_exact"] - d["error_closed_form"]) < 1e-3'
run_check "Optics: GKP X readout"     "$BASE/quantum/optics/gkp?basis=x&db=10&sigma=0.4&shots=10000&seed=3" 'abs(d["error_sampled"] - d["error_exact"]) < 0.01 and d["basis"] == "x"'
run_check "Optics: cat parity, loss"  "$BASE/quantum/optics/cat?beta=2&parity=odd&loss=0.1&shots=5000&seed=4" 'abs(d["parity_exact"] - d["parity_closed_form"]) < 1e-12 and abs(d["parity_flip_sampled"] - d["parity_flip_exact"]) < 0.03'
run_check "Optics: cat, no loss"      "$BASE/quantum/optics/cat?beta=2&loss=0&shots=2000&seed=4" 'd["parity_flip_sampled"] == 0'
run_check "Optics: heralded comb"     "$BASE/quantum/optics/comb?y=0" 'all(abs(g - 1.7725) < 0.02 for g in d["gaps"])'

echo ""
echo "Process models (BEAM):"
run_test "Entanglement (CHSH)"  "$BASE/quantum/entangle?shots=200"
run_test "Teleportation, 2 nodes" "$BASE/quantum/migrate?shots=20"

echo ""
echo "Quantum Error Correction:"
run_test "Bit-Flip (no error)"  "$BASE/quantum/qec/bitflip?shots=100"
run_test "Bit-Flip (Q1 error)"  "$BASE/quantum/qec/bitflip?shots=100&error_qubit=1"
run_fixed "Phase-Flip (no err)"  "$BASE/quantum/qec/phaseflip?shots=100"
run_fixed "Phase-Flip (Z on Q2)" "$BASE/quantum/qec/phaseflip?shots=100&error_qubit=2"
run_fixed "Shor-9 (Y on Q3)"     "$BASE/quantum/qec/shor9?shots=100&error_qubit=3"
run_fixed "Steane-7 (no err)"    "$BASE/quantum/qec/steane7?shots=100"
run_fixed "Steane-7 (Y on Q5)"   "$BASE/quantum/qec/steane7?shots=100&error_qubit=5"
run_fixed "Surface-3 (Y on Q0)"    "$BASE/quantum/qec/surface3?shots=100&error_qubit=0"

echo ""
echo "────────────────────────────────────────────────"
echo "Results: $PASS passed, $FAIL failed"
echo ""
[ "$FAIL" -eq 0 ] && exit 0 || exit 1
