import numpy as np
from collections import Counter

def solve_ode_general(coefficients):
    coeffs = np.array(coefficients, dtype=float)

    # ---- helpers ----
    EPS = 1e-7          # for deciding if imag/real ~ 0
    KEY_ROUND = 6       # rounding digits for grouping roots

    def clean0(x):
        return 0.0 if abs(x) < EPS else x

    def fmt_num(x):
        x = clean0(float(x))
        # force 1 decimal like 2.0, -0.0 becomes 0.0
        x = 0.0 if abs(x) < EPS else x
        return f"{x:.1f}"

    def x_pow(j):
        if j == 0:
            return ""
        if j == 1:
            return "x"
        return f"x^{j}"

    def exp_part(alpha):
        return f"e^({fmt_num(alpha)}x)"

    def join_mul(parts):
        parts = [p for p in parts if p]
        if not parts:
            return ""
        return "".join(parts)

    # ---- 1) characteristic roots ----
    roots = np.roots(coeffs)

    # ---- 2) group roots with tolerance (rounding key) ----
    # key uses rounded real/imag so near-equal roots group together
    root_keys = []
    for r in roots:
        a = clean0(r.real)
        b = clean0(r.imag)
        ar = round(a, KEY_ROUND)
        br = round(b, KEY_ROUND)
        root_keys.append((ar, br))

    counts = Counter(root_keys)

    # ---- 3) build solution terms ----
    terms = []
    c_index = 1
    used = set()

    # sort keys for stable output
    keys_sorted = sorted(counts.keys(), key=lambda k: (k[0], k[1]))

    for (ar, br) in keys_sorted:
        if (ar, br) in used:
            continue

        # real root
        if abs(br) < EPS:
            m = counts[(ar, br)]
            for j in range(m):
                poly = x_pow(j)
                piece = join_mul([f"C_{c_index}", poly, exp_part(ar)])
                terms.append(piece)
                c_index += 1
            used.add((ar, br))
            continue

        # complex root: handle only imag > 0 as the representative, pair with conjugate
        if br < 0:
            continue

        conj_key = (ar, round(-br, KEY_ROUND))
        if conj_key not in counts:
            # should not happen for real-coefficient ODE, but just in case:
            # treat it as a single complex exponential term
            m = counts[(ar, br)]
            for j in range(m):
                poly = x_pow(j)
                piece = join_mul([f"C_{c_index}", poly, exp_part(ar)]) + f"(cos({fmt_num(br)}x) + i sin({fmt_num(br)}x))"
                terms.append(piece)
                c_index += 1
            used.add((ar, br))
            continue

        m = min(counts[(ar, br)], counts[conj_key])

        # For each multiplicity level j: x^j e^(ax)cos(bx) and x^j e^(ax)sin(bx)
        for j in range(m):
            poly = x_pow(j)
            common = join_mul([poly, exp_part(ar)])
            cos_term = f"C_{c_index}" + common + f"cos({fmt_num(br)}x)"
            c_index += 1
            sin_term = f"C_{c_index}" + common + f"sin({fmt_num(br)}x)"
            c_index += 1
            terms.append(cos_term)
            terms.append(sin_term)

        used.add((ar, br))
        used.add(conj_key)

    # ---- final string ----
    return "y(x) = " + " + ".join(terms)

