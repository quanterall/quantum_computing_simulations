"""The review paper's "2.7 % single-photon loss" (p. 27), traced to Bartolucci et al., Nat. Commun. 14, 912 (2023),
arXiv 2101.09310, Sec. VI.2 and Appendix B — reproduced from their formulas, not from a simulation of ours.

  p0   = 1 - (1 - pfail/2) * eta**(1/pfail)                  (B1) erasure of one fusion outcome; 1/pfail photons per fusion
  penc = ([1 - (1-p0)**2]**2 + 1 - (1 - p0**2)**2) / 2        (B2) with the (2,2)-Shor-encoded resource states
  loss threshold: penc = p*_erasure, with their simulated marginal erasure thresholds 11.98 % (6-ring), 6.90 % (4-star)

The thresholds are taken from the paper (we did not simulate the fusion networks).
    ~/quanterall/.venv-tools/bin/python tools/optics/fbqc_loss.py"""
from scipy.optimize import brentq

p0 = lambda pl, pf: 1 - (1 - pf / 2) * (1 - pl) ** (1 / pf)
penc = lambda q: ((1 - (1 - q) ** 2) ** 2 + 1 - (1 - q * q) ** 2) / 2

assert abs(penc(1 / 8) - 0.043) < 5e-4                          # their baseline at eta = 1, pfail = 1/4
for net, th in (("6-ring", 0.1198), ("4-star", 0.0690)):
    for pf in (1 / 2, 1 / 4, 1 / 8):
        for enc, f in (("encoded", lambda pl: penc(p0(pl, pf))), ("unencoded", lambda pl: p0(pl, pf))):
            try:
                pl = brentq(lambda x: f(x) - th, 0, 0.5)
                print(f"{net} {enc:9} pfail {pf:.3f}: loss per photon {pl:.3%}, at least one of {round(1 / pf)} photons lost {1 - (1 - pl) ** (1 / pf):.1%}")
            except ValueError:
                print(f"{net} {enc:9} pfail {pf:.3f}: no loss tolerance (above threshold at zero loss)")
pf_star = brentq(lambda pf: penc(pf / 2) - 0.1198, 0.01, 0.99)
print(f"6-ring encoded, marginal failure threshold at zero loss: {pf_star:.1%}")
pl = brentq(lambda x: penc(p0(x, 0.25)) - 0.1198, 0, 0.5)
assert abs(pl - 0.027) < 5e-4 and abs(1 - (1 - pl) ** 4 - 0.104) < 5e-4 and abs(pf_star - 0.432) < 5e-4
print("reproduced: 2.7 % per photon = 10.4 % per four-photon fusion; 43.2 % failure threshold")
