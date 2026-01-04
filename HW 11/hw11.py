import numpy as np

def dft(f_vals, x, w):
    dx = x[1] - x[0]
    F = np.zeros_like(w, dtype=np.complex128)
    for k, wk in enumerate(w):
        integrand = f_vals * np.exp(-1j * wk * x)
        F[k] = np.trapz(integrand, x)
    return F

def idft(F_vals, w, x):
    f_rec = np.zeros_like(x, dtype=np.complex128)
    for n, xn in enumerate(x):
        integrand = F_vals * np.exp(1j * w * xn)
        f_rec[n] = (1 / (2 * np.pi)) * np.trapz(integrand, w)
    return f_rec

if __name__ == "__main__":
    L = 10.0
    Nx = 2001
    x = np.linspace(-L, L, Nx)

    def f(x):
        return np.exp(-x**2)

    f_vals = f(x)

    W = 20.0
    Nw = 2001
    w = np.linspace(-W, W, Nw)

    F_vals = dft(f_vals, x, w)
    f_back = idft(F_vals, w, x)

    err = np.max(np.abs(f_back.real - f_vals))
    print("max abs error:", err)

    for i in [0, Nx//4, Nx//2, 3*Nx//4, Nx-1]:
        print(f"x={x[i]: .3f}  f={f_vals[i]: .6f}  f_back={f_back.real[i]: .6f}")
