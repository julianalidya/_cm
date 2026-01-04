import numpy as np

def solve_ode_general(coefficients, tol=1e-6):
    """
    Solve: a0*y^(n) + a1*y^(n-1) + ... + a_{n-1}*y' + a_n*y = 0
    coefficients: [a0, a1, ..., a_n]  (a0 != 0)
    Return: a readable general solution string
    """
    coefficients = np.array(coefficients, dtype=float)
    if len(coefficients) < 2:
        return "Invalid input: need at least 2 coefficients."

    if abs(coefficients[0]) < 1e-14:
        return "Invalid input: leading coefficient cannot be 0."

    # Characteristic polynomial: a0*r^n + a1*r^(n-1) + ... + a_n = 0
    roots = np.roots(coefficients)

    # Group roots by approximate equality (to detect repeated roots)
    used = [False] * len(roots)
    groups = []  # each item: (representative_root, multiplicity)

    for i in range(len(roots)):
        if used[i]:
            continue
        r = roots[i]
        count = 1
        used[i] = True
        for j in range(i + 1, len(roots)):
            if not used[j] and abs(roots[j] - r) < tol:
                used[j] = True
                count += 1
        # Use the average of grouped roots for nicer printing
        members = [roots[k] for k in range(len(roots)) if abs(roots[k] - r) < tol]
        rep = sum(members) / len(members)
        groups.append((rep, count))

    # Sort: real first, then by real part
    groups.sort(key=lambda x: (abs(x[0].imag) > tol, x[0].real, x[0].imag))

    terms = []
    c_index = 1

    def fmt_real(x):
        # nice formatting for near-integers
        if abs(x - round(x)) < 1e-6:
            return str(int(round(x)))
        return f"{x:.6g}"

    def add_term(term):
        nonlocal c_index
        terms.append(f"C_{c_index}{term}")
        c_index += 1

    for r, mult in groups:
        if abs(r.imag) < tol:
            # real root
            a = r.real
            for k in range(mult):
                # e^{ax}, x e^{ax}, x^2 e^{ax}, ...
                if k == 0:
                    add_term(f" * exp({fmt_real(a)}x)")
                else:
                    add_term(f" * x^{k} * exp({fmt_real(a)}x)")
        else:
            # complex root a ± bi -> exp(ax)(C cos(bx) + C sin(bx))
            a = r.real
            b = abs(r.imag)
            for k in range(mult):
                # repeated complex pair -> x^k * exp(ax) * (C cos(bx) + C sin(bx))
                base = f" * exp({fmt_real(a)}x)"
                if k > 0:
                    base = f" * x^{k}" + base

                # cos term
                add_term(base + f" * cos({fmt_real(b)}x)")
                # sin term
                add_term(base + f" * sin({fmt_real(b)}x)")

    if not terms:
        return "No solution form generated."

    return "y(x) = " + " + ".join(terms)


# ---------- Test Examples ----------
if __name__ == "__main__":
    # (1) distinct real roots: y'' - 3y' + 2y = 0 -> r=1,2
    print("--- Example 1: distinct real roots ---")
    coeffs1 = [1, -3, 2]
    print("Coefficients:", coeffs1)
    print(solve_ode_general(coeffs1))
    print()

    # (2) repeated real root: y'' - 4y' + 4y = 0 -> r=2 (double)
    print("--- Example 2: repeated real root ---")
    coeffs2 = [1, -4, 4]
    print("Coefficients:", coeffs2)
    print(solve_ode_general(coeffs2))
    print()

    # (3) complex roots: y'' + y = 0 -> r=±i
    print("--- Example 3: complex roots ---")
    coeffs3 = [1, 0, 1]
    print("Coefficients:", coeffs3)
    print(solve_ode_general(coeffs3))
    print()
