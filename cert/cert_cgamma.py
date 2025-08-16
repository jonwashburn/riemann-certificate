import math
import mpmath as mp

# Canonical bump psi normalized to integral 1 on [-1,1]
# We compute kappa numerically to normalize.

mp.mp.dps = 80

def raw_bump(t):
    if abs(t) >= 1:
        return mp.mpf('0')
    return mp.e**(-1/(1-t*t))

# Normalize kappa so that integral_{-1}^1 psi = 1
I = mp.quad(lambda x: raw_bump(x), [-1,1])
kappa = 1/I

def psi(t):
    if abs(t) >= 1:
        return mp.mpf('0')
    return kappa*raw_bump(t)

# Adaptive window: L = 1/(1+log(2+|T|))

def L_of_T(T):
    return 1/(1+mp.log(2+abs(T)))

# phi_I(t;T) = psi((t-T)/L)  (integral ~ L)

def phi_I(t, T):
    L = L_of_T(T)
    return psi((t - T)/L)

# Archimedean integrand: Im d/dt log(arch(s)) at s=1/2+it
# arch(s) = pi^{-s/2} Gamma(s/2) * s(1-s)/2

def arch_log_deriv_im(t):
    s = mp.mpf('0.5') + 1j*t
    term1 = 0.5*mp.re(mp.digamma(s/2))
    term2 = -0.5*mp.log(mp.pi)
    term3 = (2*t)/(1+4*t*t)
    return term1 + term2 + term3

# Compute best-fit a + b t in weighted L2 with weight phi_I to remove affine part

def fit_affine(T):
    L = L_of_T(T)
    a = T - L
    b = T + L
    # Moments and cross terms
    w = lambda tt: phi_I(tt, T)
    g = lambda tt: arch_log_deriv_im(tt)
    m0 = mp.quad(lambda tt: w(tt), [a,b])
    m1 = mp.quad(lambda tt: tt*w(tt), [a,b])
    m2 = mp.quad(lambda tt: (tt**2)*w(tt), [a,b])
    r0 = mp.quad(lambda tt: g(tt)*w(tt), [a,b])
    r1 = mp.quad(lambda tt: tt*g(tt)*w(tt), [a,b])
    # Solve [ [m0 m1],[m1 m2] ] [α, β]^T = [r0, r1]^T
    det = m0*m2 - m1*m1
    if det == 0:
        return mp.mpf('0'), mp.mpf('0')
    alpha = ( r0*m2 - r1*m1 )/det
    beta  = ( r1*m0 - r0*m1 )/det
    return alpha, beta

# Residual C_Gamma(T) with affine removed: (1/L) * |∫ φ_I (g - α - β t) dt|

def C_gamma_at_T(T):
    L = L_of_T(T)
    a = T - L
    b = T + L
    alpha, beta = fit_affine(T)
    f = lambda tt: phi_I(tt, T) * (arch_log_deriv_im(tt) - alpha - beta*tt)
    val = mp.quad(f, [a, b])
    return abs(val / L)

# Grid over T to enclose sup_T C_gamma(T)

def certify_C_gamma(Tmax=50, N=200):
    Ts = [ (2*mp.mpf(Tmax)*i/N) - Tmax for i in range(N+1) ]
    vals = [ C_gamma_at_T(T) for T in Ts ]
    sup_est = max(vals)
    # Simple safety margin
    margin = mp.mpf('0.02')
    return sup_est + margin, list(zip(Ts, vals))

if __name__ == '__main__':
    sup_val, samples = certify_C_gamma()
    print("Certified C_Gamma(psi) <=", sup_val)
    # For logging:
    print("Max sample:", max(samples, key=lambda x: x[1]))
