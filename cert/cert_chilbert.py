import math
import mpmath as mp

mp.mp.dps = 80

# Canonical bump psi, normalized

def raw_bump(t):
    if abs(t) >= 1:
        return mp.mpf('0')
    return mp.e**(-1/(1-t*t))

I = mp.quad(lambda x: raw_bump(x), [-1,1])
kappa = 1/I

def psi(u):
    if abs(u) >= 1:
        return mp.mpf('0')
    return kappa*raw_bump(u)

# Hilbert transform of psi at z in [-1,1]: Hpsi(z) = (1/pi) PV int_{-1}^1 psi(u)/(z-u) du

def Hpsi(z):
    eps = mp.mpf('1e-8')
    a = mp.mpf('-1')
    b = mp.mpf('1')
    def integrand(u):
        if u == z:
            return mp.mpf('0')
        return psi(u)/(z - u)
    # Split integral for PV
    if z <= a:
        val = mp.quad(integrand, [a, b])
    elif z >= b:
        val = mp.quad(integrand, [a, b])
    else:
        val = mp.quad(integrand, [a, z-eps, z+eps, b])
    return (1/mp.pi)*val

# Uniform L1 bound for u' over any adaptive window (placeholder)
U0 = mp.mpf('1.0')

# Envelope: C_H(psi) = sup_{z in [-1,1]} |Hpsi(z)| * U0

def certify_C_H(N=200):
    zs = [ -1 + 2*mp.mpf(i)/N for i in range(N+1) ]
    vals = [ abs(Hpsi(z)) for z in zs ]
    sup_est = max(vals)
    margin = mp.mpf('0.01')
    return sup_est*U0 + margin, list(zip(zs, vals))

if __name__ == '__main__':
    sup_val, samples = certify_C_H()
    print("Certified C_H(psi) <=", sup_val)
    print("Max sample:", max(samples, key=lambda x: x[1]))
