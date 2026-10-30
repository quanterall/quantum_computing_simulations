"""Phase 0: KLM's NS-gate unitary (Nature 409, 46, Eq. 9) gives the NS gate with probability 1/4, and which reading
of the review's beamsplitter angles (Romero–Milburn Fig. 8) reproduces it. Run: .venv-tools/bin/python tools/optics/ns_check.py"""
import itertools
import math

import numpy as np

s2 = math.sqrt(2)
U = np.array([[1 - s2, 2 ** -0.25, math.sqrt(3 / s2 - 2)],
              [2 ** -0.25, 0.5, 0.5 - 1 / s2],
              [math.sqrt(3 / s2 - 2), 0.5 - 1 / s2, s2 - 0.5]])


def perm(M):
    return sum(np.prod([M[i, p[i]] for i in range(M.shape[0])]) for p in itertools.permutations(range(M.shape[0]))) if M.size else 1.0


def amp(M, tin, tout):
    rows = [i for i, c in enumerate(tout) for _ in range(c)]
    cols = [i for i, c in enumerate(tin) for _ in range(c)]
    return perm(M[np.ix_(rows, cols)]) / math.sqrt(np.prod([math.factorial(c) for c in tin + tout]))


def B(i, j, t):
    M = np.eye(3); c, s = math.cos(t), math.sin(t)
    M[i, i], M[i, j], M[j, i], M[j, j] = c, s, -s, c
    return M


assert np.allclose(U @ U.T, np.eye(3))
print("KLM Eq. 9, heralded amplitudes for 0, 1, 2 photons:", [round(amp(U, (n, 1, 0), (n, 1, 0)), 9) for n in range(3)])
d = math.pi / 180
for t3 in (-22.5, 22.5):
    M = B(1, 2, t3 * d) @ B(0, 1, 65.5302 * d) @ B(1, 2, 22.5 * d) @ np.diag([-1, 1, 1])
    a = [amp(M, (n, 1, 0), (n, 1, 0)) for n in range(3)]
    print(f"theta3 = {t3:+}: amplitudes {[round(x, 4) for x in a]} ->", "NS" if np.allclose(a, [a[0], a[0], -a[0]], atol=1e-3) else "not NS")
