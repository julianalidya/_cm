import numpy as np

def root(coeffs, eps=1e-14):
    coeffs = list(coeffs)

    # remove trailing zeros (highest degree terms)
    while len(coeffs) > 1 and abs(coeffs[-1]) < eps:
        coeffs.pop()

    n = len(coeffs) - 1  # degree

    if n == 0:
        return []
    if n == 1:
        return [-(coeffs[0] / coeffs[1])]

    # normalize to monic: x^n + a_{n-1} x^{n-1} + ... + a0
    lead = coeffs[-1]
    a = np.array(coeffs[:-1], dtype=np.complex128) / lead  # a0..a_{n-1}

    # companion matrix (standard form)
    C = np.zeros((n, n), dtype=np.complex128)
    C[1:, :-1] = np.eye(n - 1, dtype=np.complex128)
    C[:, -1] = -a  # last column = -[a0, a1, ..., a_{n-1}]^T

    roots = np.linalg.eigvals(C)

    # verification: plug back (Horner)
    def eval_poly(x):
        y = 0j
        for k in range(len(coeffs) - 1, -1, -1):
            y = y * x + coeffs[k]
        return y

    for r in roots:
        assert abs(eval_poly(r)) < 1e-7

    return roots

# demo
poly = [-8, 14, -7, 1]   # -8 + 14x - 7x^2 + x^3
print(root(poly))
