import cmath

def root2(a, b, c):
    if a == 0:
        raise ValueError("a must not be zero for quadratic polynomial")

    D = b*b - 4*a*c               # discriminant
    sqrtD = cmath.sqrt(D)        # bisa real atau kompleks

    x1 = (-b + sqrtD) / (2*a)
    x2 = (-b - sqrtD) / (2*a)

    # fungsi polinomial
    f = lambda x: a*x*x + b*x + c

    # verifikasi numerik
    assert cmath.isclose(f(x1), 0, abs_tol=1e-9), f"Root1 invalid: f(x1)={f(x1)}"
    assert cmath.isclose(f(x2), 0, abs_tol=1e-9), f"Root2 invalid: f(x2)={f(x2)}"

    return x1, x2
