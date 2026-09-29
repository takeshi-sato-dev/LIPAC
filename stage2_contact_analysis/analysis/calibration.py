"""Calibration of Stage 2 for time-correlated trajectories.

Frames of a molecular dynamics trajectory are not independent observations. A model
that treats them as independent returns credible intervals that are too narrow, and a
comparison of a linear model with a mixture model prefers the mixture whenever the
contact numbers are not normally distributed, whether or not binding has any effect.
This module provides the two corrections that the Stage 2 results need before they are
reported.

1. Block bootstrap. Every trajectory (protein copy) is cut into consecutive time blocks,
   and whole blocks are resampled. The interval of the effect of binding on a contact
   number then carries the correlation between frames.
2. Time-shift null. The binding-state series of every copy is shifted circularly in time
   against the contact-number series of the same copy. A shift keeps the autocorrelation
   and the distribution of both series and removes any relation between them. The
   statistics of the observed data are compared with their distribution over many shifts.

Both statistics are computed by maximum likelihood (least squares and an EM fit of the
mixture of mixture_causal_analysis.py), which makes hundreds of null replicates cheap.
"""
import numpy as np
from scipy.stats import norm


# ----------------------------------------------------------------------------- models
def fit_linear(Y, T, P, nc):
    """Y = alpha_p + beta_p T + e, common sigma. Returns log-likelihood and beta_p."""
    n = len(Y)
    X = np.zeros((n, 2 * nc))
    X[np.arange(n), P] = 1.0
    X[np.arange(n), nc + P] = T
    b, *_ = np.linalg.lstsq(X, Y, rcond=None)
    r = Y - X @ b
    s = max(r.std(), 1e-6)
    return norm.logpdf(r, 0, s).sum(), b[nc:]


def fit_mixture(Y, T, P, nc, iters=300):
    """Mixture model of mixture_causal_analysis.py, fitted by EM.
    Unbound frames: N(alpha_p, s). Bound frames: (1-pi) N(alpha_p + bd_p, s) + pi N(alpha_p + bd_p + bc, s).
    Returns log-likelihood, beta_coop and pi."""
    a = np.array([Y[(P == c) & (T == 0)].mean() if ((P == c) & (T == 0)).any() else Y[P == c].mean() for c in range(nc)])
    bd = np.array([Y[(P == c) & (T == 1)].mean() - a[c] if ((P == c) & (T == 1)).any() else 0.0 for c in range(nc)])
    s = max(Y.std(), 1e-6); bc = s; pi = 0.3
    for _ in range(iters):
        mu = a[P] + bd[P] * T
        l0 = np.log(1 - pi) + norm.logpdf(Y, mu, s)
        l1 = np.log(pi) + norm.logpdf(Y, mu + bc, s)
        w = np.where(T == 1, 1.0 / (1.0 + np.exp(np.clip(l0 - l1, -50, 50))), 0.0)
        if (T == 1).any():
            pi = float(np.clip(w[T == 1].mean(), 1e-4, 1 - 1e-4))
        R = Y - w * bc
        for c in range(nc):
            u = (P == c) & (T == 0); b = (P == c) & (T == 1)
            if u.any():
                a[c] = R[u].mean()
            if b.any():
                bd[c] = R[b].mean() - a[c]
        mu = a[P] + bd[P] * T
        bc = float((w * (Y - mu)).sum() / max(w.sum(), 1e-9))
        s = max(np.sqrt(((1 - w) * (Y - mu) ** 2 + w * (Y - mu - bc) ** 2).mean()), 1e-6)
    mu = a[P] + bd[P] * T
    ll = np.where(T == 1,
                  np.logaddexp(np.log(1 - pi) + norm.logpdf(Y, mu, s), np.log(pi) + norm.logpdf(Y, mu + bc, s)),
                  norm.logpdf(Y, mu, s)).sum()
    return ll, bc, pi


def statistics(Y, T, P, nc):
    """Mean effect over copies (linear model) and the AIC advantage of the mixture model."""
    l1, beta = fit_linear(Y, T, P, nc)
    l2, bc, pi = fit_mixture(Y, T, P, nc)
    d_aic = 2 * (l2 - l1) - 2 * 2          # the mixture has two more parameters (beta_coop, pi)
    return dict(beta_mean=float(beta.mean()), beta=beta, d_aic=float(d_aic), beta_coop=bc, pi=pi)


# ----------------------------------------------------------------------------- helpers
def autocorrelation_time(y, max_lag=None):
    """Integrated autocorrelation time in frames (initial positive sequence)."""
    y = np.asarray(y, float) - np.mean(y)
    n = len(y); max_lag = max_lag or n // 4
    f = np.fft.rfft(y, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:max_lag] / (y.var() * n)
    tau = 1.0
    for k in range(1, max_lag):
        if ac[k] <= 0:
            break
        tau += 2 * ac[k]
    return tau


def _split(Y, T, P, nc):
    return [(Y[P == c], T[P == c]) for c in range(nc)]


def _join(parts):
    Y = np.concatenate([p[0] for p in parts]); T = np.concatenate([p[1] for p in parts])
    P = np.concatenate([np.full(len(p[0]), c) for c, p in enumerate(parts)])
    return Y, T, P


# ----------------------------------------------------------------------------- 1. block bootstrap
def choose_n_blocks(Y, P, nc, lo=5, hi=20):
    """Number of blocks per copy: each block about two autocorrelation times long of the
    most correlated copy, between lo and hi blocks."""
    L = min((P == c).sum() for c in range(nc))
    tau = max(autocorrelation_time(Y[P == c]) for c in range(nc))
    return int(np.clip(L // max(2 * tau, 1), lo, hi))


def block_bootstrap(Y, T, P, nc, n_blocks=None, n_boot=1000, seed=0):
    """Per-copy and mean effect of binding (bound minus unbound mean contact number),
    with 95% intervals from resampling whole time blocks within each copy.
    n_blocks=None chooses the number from the autocorrelation time of the contact numbers."""
    rng = np.random.default_rng(seed)
    if n_blocks is None:
        n_blocks = choose_n_blocks(Y, P, nc)
    parts = _split(Y, T, P, nc)
    blocks = [np.array_split(np.arange(len(y)), n_blocks) for y, _ in parts]
    def effect(y, t):
        return y[t == 1].mean() - y[t == 0].mean() if (t == 1).any() and (t == 0).any() else np.nan
    obs = np.array([effect(y, t) for y, t in parts])
    boot = np.full((n_boot, nc), np.nan)
    for i in range(n_boot):
        for c, (y, t) in enumerate(parts):
            idx = np.concatenate([blocks[c][j] for j in rng.integers(0, n_blocks, n_blocks)])
            boot[i, c] = effect(y[idx], t[idx])
    lo, hi = np.nanpercentile(boot, [2.5, 97.5], axis=0)
    m = np.nanmean(boot, axis=1)
    return dict(effect=obs, ci_low=lo, ci_high=hi, mean_effect=float(np.nanmean(obs)),
                mean_ci=tuple(np.nanpercentile(m, [2.5, 97.5])),
                block_frames=[len(b[0]) for b in blocks],
                tau_frames=[autocorrelation_time(y) for y, _ in parts])


# ----------------------------------------------------------------------------- 2. time-shift null
def time_shift_null(Y, T, P, nc, n_null=200, min_frac=0.1, seed=0):
    """Distribution of the statistics when the binding series of every copy is shifted
    circularly by a random amount between min_frac and 1 - min_frac of its length."""
    rng = np.random.default_rng(seed)
    parts = _split(Y, T, P, nc)
    obs = statistics(Y, T, P, nc)
    null_beta, null_aic = [], []
    for _ in range(n_null):
        shifted = []
        for y, t in parts:
            L = len(t); k = rng.integers(int(min_frac * L), int((1 - min_frac) * L))
            shifted.append((y, np.roll(t, k)))
        st = statistics(*_join(shifted), nc)
        null_beta.append(st['beta_mean']); null_aic.append(st['d_aic'])
    null_beta = np.array(null_beta); null_aic = np.array(null_aic)
    p_beta = (1 + np.sum(np.abs(null_beta) >= abs(obs['beta_mean']))) / (n_null + 1)
    p_aic = (1 + np.sum(null_aic >= obs['d_aic'])) / (n_null + 1)
    return dict(observed=obs, null_beta=null_beta, null_d_aic=null_aic, p_effect=float(p_beta),
                p_mixture=float(p_aic), d_aic_null_95=float(np.percentile(null_aic, 95)))


def calibrated_classification(Y, T, P, nc, n_null=200, n_blocks=None, alpha=0.05, seed=0):
    """Classification of one lipid type that holds against time-correlated frames.

    effect      : the mean effect over copies lies outside the time-shift null (p < alpha)
                  and its block-bootstrap interval excludes zero.
    cooperative : in addition, the mixture model beats the linear model by more than it does
                  in the time-shift null (p < alpha).
    """
    Y = np.asarray(Y, float); T = np.asarray(T, int); P = np.asarray(P, int)
    bb = block_bootstrap(Y, T, P, nc, n_blocks=n_blocks, seed=seed)
    ns = time_shift_null(Y, T, P, nc, n_null=n_null, seed=seed)
    lo, hi = bb['mean_ci']
    effect = ns['p_effect'] < alpha and not (lo <= 0 <= hi)
    coop = effect and ns['p_mixture'] < alpha
    label = 'cooperative' if coop else ('linear' if effect else 'no detectable effect')
    return dict(classification=label, block_bootstrap=bb, time_shift_null=ns)
