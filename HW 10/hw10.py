import numpy as np

# ------------------------------------------------------------
# 1) Determinant by recursion (Laplace expansion)
# ------------------------------------------------------------
def det_recursive(A, eps=1e-12):
    A = np.array(A, dtype=float)
    n, m = A.shape
    if n != m:
        raise ValueError("det_recursive: matrix must be square")
    if n == 0:
        return 1.0
    if n == 1:
        return float(A[0, 0])
    if n == 2:
        return float(A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0])

    # expand along row 0 (can be expensive for large n)
    d = 0.0
    for j in range(n):
        if abs(A[0, j]) < eps:
            continue
        minor = np.delete(np.delete(A, 0, axis=0), j, axis=1)
        cofactor = ((-1) ** j) * A[0, j]
        d += cofactor * det_recursive(minor, eps=eps)
    return float(d)


# ------------------------------------------------------------
# 2) LU decomposition (with partial pivoting): P @ A = L @ U
#    and determinant using LU
# ------------------------------------------------------------
def lu_decompose(A, eps=1e-12):
    A = np.array(A, dtype=float)
    n, m = A.shape
    if n != m:
        raise ValueError("lu_decompose: matrix must be square")

    U = A.copy()
    L = np.eye(n, dtype=float)
    P = np.eye(n, dtype=float)
    swap_count = 0

    for k in range(n):
        # pivot
        pivot_row = k + np.argmax(np.abs(U[k:, k]))
        if abs(U[pivot_row, k]) < eps:
            raise ValueError("lu_decompose: singular or near-singular matrix")

        if pivot_row != k:
            U[[k, pivot_row], :] = U[[pivot_row, k], :]
            P[[k, pivot_row], :] = P[[pivot_row, k], :]
            if k > 0:
                L[[k, pivot_row], :k] = L[[pivot_row, k], :k]
            swap_count += 1

        # elimination
        for i in range(k + 1, n):
            L[i, k] = U[i, k] / U[k, k]
            U[i, k:] = U[i, k:] - L[i, k] * U[k, k:]
            U[i, k] = 0.0

    return P, L, U, swap_count


def det_via_lu(A, eps=1e-12):
    P, L, U, swap_count = lu_decompose(A, eps=eps)
    detP = -1.0 if (swap_count % 2 == 1) else 1.0
    detU = float(np.prod(np.diag(U)))
    return detP * detU


# ------------------------------------------------------------
# 3) Verify reconstructions: LU / eigen / SVD
# ------------------------------------------------------------
def verify_lu(A):
    P, L, U, _ = lu_decompose(A)
    A_rec = np.linalg.inv(P) @ (L @ U)  # since P @ A = L @ U => A = P^{-1} L U
    err = np.linalg.norm(A - A_rec, ord="fro")
    return A_rec, err


def verify_eigendecomp(A, symmetric_only=True):
    A = np.array(A, dtype=float)
    if symmetric_only:
        if not np.allclose(A, A.T, atol=1e-10):
            raise ValueError("verify_eigendecomp: matrix must be symmetric if symmetric_only=True")
        w, V = np.linalg.eigh(A)  # stable for symmetric
        A_rec = V @ np.diag(w) @ V.T
    else:
        w, V = np.linalg.eig(A)
        # For general matrices, reconstruction uses V^{-1} (if diagonalizable)
        Vinv = np.linalg.inv(V)
        A_rec = V @ np.diag(w) @ Vinv

    err = np.linalg.norm(A - A_rec, ord="fro")
    return A_rec, err


def verify_svd(A):
    U, s, Vt = np.linalg.svd(A, full_matrices=False)
    A_rec = U @ np.diag(s) @ Vt
    err = np.linalg.norm(A - A_rec, ord="fro")
    return A_rec, err


# ------------------------------------------------------------
# 4) Build SVD using eigen decomposition (works for any A)
#    - eigen of (A^T A) gives V and singular values
#    - U computed as U = A V Sigma^{-1}
# ------------------------------------------------------------
def svd_via_eigen(A, eps=1e-12):
    A = np.array(A, dtype=float)
    m, n = A.shape

    AtA = A.T @ A
    # symmetric PSD -> use eigh
    vals, V = np.linalg.eigh(AtA)

    # sort by descending eigenvalue
    idx = np.argsort(vals)[::-1]
    vals = vals[idx]
    V = V[:, idx]

    # singular values are sqrt(eigenvalues) (clip for numerical stability)
    s = np.sqrt(np.clip(vals, 0.0, None))

    # keep only non-zero singular values to avoid divide-by-zero
    r = np.sum(s > eps)
    s_r = s[:r]
    V_r = V[:, :r]

    # U = A V Sigma^{-1}
    U_r = A @ V_r
    U_r = U_r / s_r  # broadcast division by each singular value

    # Orthonormalize U_r (numerical cleanup)
    # QR gives orthonormal columns
    U_r, _ = np.linalg.qr(U_r)

    # Build Sigma (m x n, but we usually use compact)
    S_r = np.diag(s_r)
    Vt_r = V_r.T
    return U_r, s_r, Vt_r


# ------------------------------------------------------------
# 5) PCA using SVD
#    - center data
#    - SVD of centered X: Xc = U S V^T
#    - principal directions = columns of V
# ------------------------------------------------------------
def pca_svd(X, k=None, center=True):
    X = np.array(X, dtype=float)
    if center:
        mu = X.mean(axis=0, keepdims=True)
        Xc = X - mu
    else:
        mu = np.zeros((1, X.shape[1]))
        Xc = X

    U, s, Vt = np.linalg.svd(Xc, full_matrices=False)

    # explained variance (using sample covariance scaling)
    # variance along PC_i = (s_i^2) / (n_samples - 1)
    n_samples = X.shape[0]
    var = (s ** 2) / max(n_samples - 1, 1)
    var_ratio = var / (var.sum() if var.sum() > 0 else 1.0)

    if k is None:
        k = X.shape[1]
    k = int(k)

    components = Vt[:k, :]          # shape (k, n_features)
    scores = Xc @ components.T      # projected data (n_samples, k)

    return {
        "mean": mu.reshape(-1),
        "components": components,
        "scores": scores,
        "singular_values": s,
        "explained_variance": var,
        "explained_variance_ratio": var_ratio
    }


# ------------------------------------------------------------
# Demo / quick self-test
# ------------------------------------------------------------
if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)

    # Square matrix for det, LU, eigen
    A = np.array([
        [2, 1, 3],
        [4, 1, 6],
        [1, 0, 1]
    ], dtype=float)

    print("A=\n", A)

    print("\n[1] det_recursive(A) =", det_recursive(A))
    print("[1] numpy det(A)      =", float(np.linalg.det(A)))

    print("\n[2] det_via_lu(A)    =", det_via_lu(A))

    A_lu_rec, lu_err = verify_lu(A)
    print("\n[3] LU recon error (Frobenius) =", lu_err)

    # Eigen verification needs symmetric matrix (use a symmetric one)
    S = np.array([
        [4, 1, 2],
        [1, 3, 0],
        [2, 0, 2]
    ], dtype=float)
    S_rec, eig_err = verify_eigendecomp(S, symmetric_only=True)
    print("\n[3] Symmetric eigen recon error =", eig_err)

    # SVD verification works for any matrix (even non-square)
    B = np.array([
        [1, 2, 3],
        [4, 5, 6]
    ], dtype=float)
    B_rec, svd_err = verify_svd(B)
    print("\n[3] SVD recon error =", svd_err)

    # [4] SVD via eigen
    Ue, se, Vte = svd_via_eigen(B)
    B_rec2 = Ue @ np.diag(se) @ Vte
    err2 = np.linalg.norm(B - B_rec2, ord="fro")
    print("\n[4] SVD via eigen recon error =", err2)

    # [5] PCA demo
    X = np.array([
        [2.5, 2.4],
        [0.5, 0.7],
        [2.2, 2.9],
        [1.9, 2.2],
        [3.1, 3.0],
        [2.3, 2.7],
        [2.0, 1.6],
        [1.0, 1.1],
        [1.5, 1.6],
        [1.1, 0.9],
    ])
    pca = pca_svd(X, k=1, center=True)
    print("\n[5] PCA mean =", pca["mean"])
    print("[5] first PC (components) =", pca["components"])
    print("[5] explained variance ratio =", pca["explained_variance_ratio"][:3])
    print("[5] projected scores (first 5) =\n", pca["scores"][:5])
