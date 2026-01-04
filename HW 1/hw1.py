import math
import random

def df(f, x, h=1e-6):
    # central difference: (f(x+h)-f(x-h)) / (2h)
    return (f(x + h) - f(x - h)) / (2 * h)

def integral(f, a, b, n=20000):
    # composite trapezoidal rule
    if a == b:
        return 0.0
    if b < a:
        return -integral(f, b, a, n)

    h = (b - a) / n
    s = 0.5 * (f(a) + f(b))
    x = a
    for _ in range(1, n):
        x += h
        s += f(x)
    return s * h

def theorem1(f, x, tol=1e-4):
    g = lambda t: integral(f, 0.0, t)  # g(x) = ∫_0^x f(t) dt
    lhs = df(g, x)
    rhs = f(x)
    if abs(lhs - rhs) > tol:
        raise AssertionError(f"Fail at x={x}: df(∫f)={lhs}, f(x)={rhs}, diff={abs(lhs-rhs)}")
    return lhs, rhs

if __name__ == "__main__":
    tests = [
        ("sin", lambda x: math.sin(x), (-2.0, 2.0)),
        ("x^2", lambda x: x*x, (-2.0, 2.0)),
        ("exp", lambda x: math.exp(x), (-1.0, 1.0)),
        ("1/(1+x^2)", lambda x: 1.0/(1.0 + x*x), (-2.0, 2.0)),
    ]

    random.seed(0)

    for name, f, (L, R) in tests:
        for _ in range(20):
            x = random.uniform(L, R)
            lhs, rhs = theorem1(f, x, tol=2e-4)
        print(f"[OK] {name} passed (20 random points in [{L}, {R}])")
