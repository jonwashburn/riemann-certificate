# Certification Plan (Path B)

This folder contains minimal, certified computations needed to close the Final Certificate inequality without prime sums. We certify:

1. C_Gamma(psi): Archimedean windowed constant on adaptive intervals
2. C_H(psi): Hilbert-window constant via direct pairing and adaptive L^1 bound

We use interval arithmetic (via mpmath/arb/python-flint if available) to produce rigorous enclosures.

Outline:
- `cert_cgamma.py`: computes an upper bound for C_Gamma(psi) on a grid of T with L(T)=1/(1+log(2+|T|)) and takes a supremum enclosure.
- `cert_chilbert.py`: computes an upper bound for sup_{z\in[-1,1]} |(1/π) PV\int_{-1}^1 ψ(u)/(z-u) du| and multiplies by a uniform L^1 bound U0 for u' on an adaptive window.

No prime sums are used.

## Current certification (canonical bump ψ)

- C_Gamma(ψ) ≤ 2.392 (sup over T in [-50,50], margin = 0.02)
- C_H(ψ) ≤ 0.671 (sup over z in [-1,1], margin = 0.01), with U0 = 1.0 (placeholder)

Slack calculation:
- π/2 ≈ 1.571; C_Gamma + C_H ≈ 3.063 > π/2 → no slack for κ at present.

## Next steps to close via Path B

- Reduce C_H(ψ): CH scales linearly with U0. If the uniform-eps theorem yields a certified U0 ≤ 0.1 on adaptive windows, then C_H(ψ) ≤ 0.067, giving C_Gamma + C_H ≈ 2.459.
- Reduce C_Gamma(ψ): two options
  1) Analytic: incorporate the affine calibration (subtract a + b t) of the archimedean integrand justified by the affine-embedding lemma (main text), then certify the residual. This should significantly lower C_Gamma.
  2) Smoother ψ: test a family of C^∞ windows to seek a lower envelope for sup_T (1/L)∫ φ_I(t)·g(t) dt.

Once C_Gamma(ψ) + C_H(ψ) < π/2, choose κ ≤ 0.5 (π/2 − C_Gamma − C_H) to close PSC ⇒ (P+) ⇒ Schur/PSD ⇒ RH.
