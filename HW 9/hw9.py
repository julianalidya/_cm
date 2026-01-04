import math
import random
from typing import List, Tuple, Dict

EPS = 1e-15


# -------------------------
# (1) Probability p^10000
# -------------------------
def prob_all_heads(n: int = 10000, p: float = 0.5) -> float:
    return p ** n


# -------------------------
# (2) log(p^n) = n log(p)
# -------------------------
def log_prob_all_heads(n: int = 10000, p: float = 0.5, base: str = "e") -> float:
    if base == "e":
        return n * math.log(p)
    elif base == "2":
        return n * math.log(p, 2)
    elif base == "10":
        return n * math.log10(p)
    else:
        raise ValueError("base must be 'e', '2', or '10'")


# -------------------------
# Helpers for distributions
# -------------------------
def normalize(dist: List[float]) -> List[float]:
    s = sum(dist)
    if s <= 0:
        raise ValueError("Distribution sum must be > 0")
    return [x / s for x in dist]


def safe_log(x: float, base: float = 2.0) -> float:
    # avoid log(0)
    x = max(x, EPS)
    return math.log(x, base)


# -------------------------
# (3) Entropy, Cross-Entropy, KL
# -------------------------
def entropy(p: List[float], base: float = 2.0) -> float:
    p = normalize(p)
    return -sum(pi * safe_log(pi, base) for pi in p)


def cross_entropy(p: List[float], q: List[float], base: float = 2.0) -> float:
    p = normalize(p)
    q = normalize(q)
    return -sum(pi * safe_log(qi, base) for pi, qi in zip(p, q))


def kl_divergence(p: List[float], q: List[float], base: float = 2.0) -> float:
    p = normalize(p)
    q = normalize(q)
    return sum(pi * (safe_log(pi, base) - safe_log(qi, base)) for pi, qi in zip(p, q))


# -------------------------
# (3) Mutual Information
# -------------------------
def mutual_information(joint: List[List[float]], base: float = 2.0) -> float:
    # joint is a matrix P(X=i, Y=j)
    total = sum(sum(row) for row in joint)
    if total <= 0:
        raise ValueError("Joint distribution sum must be > 0")

    Pxy = [[v / total for v in row] for row in joint]
    Px = [sum(row) for row in Pxy]
    Py = [sum(Pxy[i][j] for i in range(len(Pxy))) for j in range(len(Pxy[0]))]

    mi = 0.0
    for i in range(len(Pxy)):
        for j in range(len(Pxy[0])):
            pxy = Pxy[i][j]
            if pxy > 0:
                mi += pxy * (safe_log(pxy, base) - safe_log(Px[i] * Py[j], base))
    return mi


# -------------------------
# (4) Verify cross_entropy(p,p) <= cross_entropy(p,q)
# -------------------------
def verify_cross_entropy_property(trials: int = 5) -> None:
    for t in range(trials):
        # create a random p and q over 4 symbols
        p = normalize([random.random() for _ in range(4)])
        q = normalize([random.random() for _ in range(4)])

        cpp = cross_entropy(p, p)
        cpq = cross_entropy(p, q)

        print(f"Trial {t+1}")
        print("p =", [round(x, 4) for x in p])
        print("q =", [round(x, 4) for x in q])
        print("cross_entropy(p,p) =", cpp)
        print("cross_entropy(p,q) =", cpq)
        print("Is cross_entropy(p,p) <= cross_entropy(p,q)?", cpp <= cpq + 1e-12)
        print("-" * 50)


# -------------------------
# (5) (7,4) Hamming Code
# Positions (1-indexed): 1 2 3 4 5 6 7
# Bits:                 p1 p2 d1 p4 d2 d3 d4
# Parity checks:
# p1 covers positions 1,3,5,7
# p2 covers positions 2,3,6,7
# p4 covers positions 4,5,6,7
# Using even parity
# -------------------------
def hamming74_encode(data4: str) -> str:
    if len(data4) != 4 or any(c not in "01" for c in data4):
        raise ValueError("data4 must be a 4-bit string like '1010'")

    d1, d2, d3, d4 = map(int, data4)

    # parity bits (even parity)
    p1 = (d1 + d2 + d4) % 2
    p2 = (d1 + d3 + d4) % 2
    p4 = (d2 + d3 + d4) % 2

    code7 = [p1, p2, d1, p4, d2, d3, d4]
    return "".join(str(b) for b in code7)


def hamming74_syndrome(code7: List[int]) -> int:
    # returns error position (1..7), or 0 if no error
    # s1 checks parity over 1,3,5,7
    s1 = (code7[0] + code7[2] + code7[4] + code7[6]) % 2
    # s2 checks parity over 2,3,6,7
    s2 = (code7[1] + code7[2] + code7[5] + code7[6]) % 2
    # s4 checks parity over 4,5,6,7
    s4 = (code7[3] + code7[4] + code7[5] + code7[6]) % 2

    # syndrome bits form a binary number: s4 s2 s1
    return s1 + 2 * s2 + 4 * s4


def hamming74_decode(code7_str: str) -> Tuple[str, Dict[str, int]]:
    if len(code7_str) != 7 or any(c not in "01" for c in code7_str):
        raise ValueError("code7_str must be a 7-bit string like '0110011'")

    code7 = [int(c) for c in code7_str]
    err_pos = hamming74_syndrome(code7)

    corrected = code7[:]
    if err_pos != 0:
        corrected[err_pos - 1] ^= 1  # flip the bit at err_pos

    # extract data bits positions 3,5,6,7 (1-indexed)
    d1 = corrected[2]
    d2 = corrected[4]
    d3 = corrected[5]
    d4 = corrected[6]
    data4 = f"{d1}{d2}{d3}{d4}"

    info = {
        "error_position": err_pos,
        "corrected": int(err_pos != 0)
    }
    return data4, info


def flip_one_bit(code7: str, pos: int) -> str:
    # pos is 1..7
    if pos < 1 or pos > 7:
        raise ValueError("pos must be 1..7")
    bits = list(code7)
    bits[pos - 1] = "1" if bits[pos - 1] == "0" else "0"
    return "".join(bits)


# -------------------------
# Demo / main
# -------------------------
if __name__ == "__main__":
    print("=== (1) p^10000 ===")
    val = prob_all_heads(10000, 0.5)
    print("0.5^10000 =", val)  # likely prints 0.0 due to underflow in float

    print("\n=== (2) log(0.5^10000) ===")
    print("ln:", log_prob_all_heads(10000, 0.5, base="e"))
    print("log2:", log_prob_all_heads(10000, 0.5, base="2"))
    print("log10:", log_prob_all_heads(10000, 0.5, base="10"))

    print("\n=== (3) entropy / cross-entropy / KL ===")
    p = [0.5, 0.25, 0.25]
    q = [0.4, 0.3, 0.3]
    print("H(p) =", entropy(p))
    print("H(p,q) =", cross_entropy(p, q))
    print("KL(p||q) =", kl_divergence(p, q))

    print("\n=== (3) mutual information ===")
    joint = [
        [0.25, 0.25],
        [0.10, 0.40]
    ]
    print("I(X;Y) =", mutual_information(joint))

    print("\n=== (4) verify cross-entropy property ===")
    verify_cross_entropy_property(trials=5)

    print("\n=== (5) Hamming (7,4) encode/decode ===")
    data = "1011"
    code = hamming74_encode(data)
    print("data:", data)
    print("encoded:", code)

    # simulate a 1-bit error
    corrupted = flip_one_bit(code, pos=6)
    print("corrupted (bit 6 flipped):", corrupted)

    decoded, info = hamming74_decode(corrupted)
    print("decoded:", decoded)
    print("info:", info)
