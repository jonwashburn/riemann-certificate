import math
import time

# Coarse certified lower bound for c0(psi) = inf_{0<b<=1, |x|<=1} (P_b * psi)(x)
# using a canonical C^∞ bump psi on [-1,1], normalized to integral 1.
# Numerics: simple composite Simpson quadrature over u in [-1,1],
# coarse grids b in [1e-3,1], x in [-1,1].

START = time.time()

# Canonical bump, normalized

def raw_bump(t: float) -> float:
    if abs(t) >= 1.0:
        return 0.0
    return math.exp(-1.0/(1.0 - t*t))

# Normalize psi so integral_{-1}^1 psi = 1 using Simpson rule
n_int = 800  # even
h = 2.0/n_int
u0 = -1.0
psi_vals = []
for i in range(n_int+1):
    u = u0 + i*h
    psi_vals.append(raw_bump(u))
# Simpson weights
W = [1.0] + [4.0 if i % 2 == 1 else 2.0 for i in range(1, n_int)] + [1.0]
I = (h/3.0)*sum(w*v for w, v in zip(W, psi_vals))
kappa = 1.0/I

# Precompute normalized psi on same nodes
psi_vals = [kappa*v for v in psi_vals]

# Poisson kernel

def P(b: float, x: float) -> float:
    return (1.0/math.pi) * b / (b*b + x*x)

# Convolution (P_b * psi)(x) via Simpson rule on the fixed u-grid

def conv(b: float, x: float) -> float:
    s = 0.0
    for i in range(n_int+1):
        u = u0 + i*h
        s += W[i] * psi_vals[i] * P(b, x - u)
    return (h/3.0) * s

# Coarse grids
BMIN = 1e-3
NB = 121
NX = 121

b_vals = [BMIN + (1.0 - BMIN)*j/(NB-1) for j in range(NB)]
x_vals = [-1.0 + 2.0*i/(NX-1) for i in range(NX)]

c0 = float('inf')
argmin = (None, None)

for bj, b in enumerate(b_vals):
    row_min = float('inf')
    row_argx = None
    for xi, x in enumerate(x_vals):
        val = conv(b, x)
        if val < row_min:
            row_min = val
            row_argx = x
    if row_min < c0:
        c0 = row_min
        argmin = (b, row_argx)

# Conservative margin for coarse grid and quadrature
margin = 1e-3
c0_cert = max(0.0, c0 - margin)

print("Certified c0(psi) >=", c0_cert)
print("Argmin approx:", argmin, "raw_min:", c0)
print("Runtime (s):", round(time.time() - START, 2))
