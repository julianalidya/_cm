import cmath

def root3(a, b, c, d):
    if a == 0:
        raise ValueError("a must not be zero for cubic polynomial")

    # 1. Normalize
    A = b / a
    B = c / a
    C = d / a

    # 2. Depressed cubic: t^3 + pt + q = 0
    p = B - (A*A)/3
    q = (2*A*A*A)/27 - (A*B)/3 + C

    # 3. Discriminant
    delta = (q/2)**2 + (p/3)**3

    # 4. Cardano parts
    s = cmath.sqrt(delta)
    u = (-q/2 + s) ** (1/3)
    v = (-q/2 - s) ** (1/3)

    # 5. Cube roots of unity
    omega = cmath.exp(2j * cmath.pi / 3)
    omega2 = omega**2

    # 6. Three roots
    t1 = u + v
    t2 = u*omega + v*omega2
    t3 = u*omega2 + v*omega

    # shift back
    x1 = t1 - A/3
    x2 = t2 - A/3
    x3 = t3 - A/3

    # optional verification
    f = lambda x: a*x**3 + b*x**2 + c*x + d
    for x in (x1, x2, x3):
        assert cmath.isclose(f(x), 0, abs_tol=1e-8)

    return x1, x2, x3
