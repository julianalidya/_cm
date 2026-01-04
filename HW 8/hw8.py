import numpy as np
from math import sqrt
from scipy.stats import norm, t


def one_sample_z_test(x, mu0, sigma, alternative="two-sided"):
    """
    One-sample z-test: population std sigma is known.
    alternative: "two-sided", "less", "greater"
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    xbar = x.mean()

    z = (xbar - mu0) / (sigma / sqrt(n))

    if alternative == "two-sided":
        p = 2 * (1 - norm.cdf(abs(z)))
    elif alternative == "greater":
        p = 1 - norm.cdf(z)
    elif alternative == "less":
        p = norm.cdf(z)
    else:
        raise ValueError("alternative must be 'two-sided', 'less', or 'greater'")

    return {"stat": z, "p_value": p, "n": n, "mean": xbar}


def one_sample_t_test(x, mu0, alternative="two-sided"):
    """
    One-sample t-test: population std unknown (use sample std).
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    xbar = x.mean()
    s = x.std(ddof=1)

    t_stat = (xbar - mu0) / (s / sqrt(n))
    df = n - 1

    if alternative == "two-sided":
        p = 2 * (1 - t.cdf(abs(t_stat), df))
    elif alternative == "greater":
        p = 1 - t.cdf(t_stat, df)
    elif alternative == "less":
        p = t.cdf(t_stat, df)
    else:
        raise ValueError("alternative must be 'two-sided', 'less', or 'greater'")

    return {"stat": t_stat, "df": df, "p_value": p, "n": n, "mean": xbar, "s": s}


def two_sample_t_test(x1, x2, alternative="two-sided", equal_var=False):
    """
    Independent two-sample t-test.
    equal_var=False -> Welch's t-test (recommended)
    equal_var=True  -> pooled-variance t-test
    """
    x1 = np.asarray(x1, dtype=float)
    x2 = np.asarray(x2, dtype=float)

    n1, n2 = x1.size, x2.size
    m1, m2 = x1.mean(), x2.mean()
    s1, s2 = x1.std(ddof=1), x2.std(ddof=1)

    if equal_var:
        # pooled
        sp2 = ((n1 - 1) * s1**2 + (n2 - 1) * s2**2) / (n1 + n2 - 2)
        se = sqrt(sp2 * (1 / n1 + 1 / n2))
        t_stat = (m1 - m2) / se
        df = n1 + n2 - 2
    else:
        # Welch
        se2 = s1**2 / n1 + s2**2 / n2
        se = sqrt(se2)
        t_stat = (m1 - m2) / se
        df = (se2**2) / ((s1**2 / n1) ** 2 / (n1 - 1) + (s2**2 / n2) ** 2 / (n2 - 1))

    if alternative == "two-sided":
        p = 2 * (1 - t.cdf(abs(t_stat), df))
    elif alternative == "greater":
        p = 1 - t.cdf(t_stat, df)
    elif alternative == "less":
        p = t.cdf(t_stat, df)
    else:
        raise ValueError("alternative must be 'two-sided', 'less', or 'greater'")

    return {
        "stat": t_stat,
        "df": df,
        "p_value": p,
        "n1": n1,
        "n2": n2,
        "mean1": m1,
        "mean2": m2,
        "s1": s1,
        "s2": s2,
        "equal_var": equal_var,
    }


def paired_t_test(before, after, alternative="two-sided"):
    """
    Paired t-test = one-sample t-test on differences (after - before).
    """
    before = np.asarray(before, dtype=float)
    after = np.asarray(after, dtype=float)
    if before.size != after.size:
        raise ValueError("before and after must have the same length")

    d = after - before
    res = one_sample_t_test(d, mu0=0.0, alternative=alternative)
    res["mean_diff"] = d.mean()
    return res


if __name__ == "__main__":
    # Example usage:
    x = [52, 50, 49, 51, 53, 52, 48, 50]
    print("One-sample Z-test:", one_sample_z_test(x, mu0=50, sigma=2.0))
    print("One-sample T-test:", one_sample_t_test(x, mu0=50))

    a = [10, 12, 9, 11, 13]
    b = [8, 7, 9, 10, 6]
    print("Two-sample T-test (Welch):", two_sample_t_test(a, b, equal_var=False))
    print("Two-sample T-test (pooled):", two_sample_t_test(a, b, equal_var=True))

    before = [70, 72, 68, 75, 71]
    after = [74, 73, 70, 78, 72]
    print("Paired T-test:", paired_t_test(before, after))
