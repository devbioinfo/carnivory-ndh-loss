#!/usr/bin/env python3
"""Marginal likelihood of each clock model by 8-point Gauss-Legendre quadrature of E_beta[lnL] over beta in (0,1),
with batch-means standard errors; Bayes factor and posterior model probabilities (equal prior odds)."""
import glob, math, os
W = [0.1012285363, 0.2223810345, 0.3137066459, 0.3626837834, 0.3626837834, 0.3137066459, 0.2223810345, 0.1012285363]
def lnl(d):
    v = []
    for f in sorted(glob.glob(f'{d}/mcmc_part*.txt')) + [f'{d}/mcmc.txt']:
        if not os.path.exists(f): continue
        with open(f) as fh:
            h = fh.readline().rstrip('\n').split('\t'); k = h.index('lnL')
            for l in fh:
                p = l.rstrip('\n').split('\t')
                if len(p) == len(h):
                    try: v.append(float(p[k]))
                    except ValueError: pass
    return v
def mean_se(v, nb=20):
    n = len(v); m = sum(v) / n; b = n // nb
    bm = [sum(v[i*b:(i+1)*b]) / b for i in range(nb)]
    return m, math.sqrt(sum((x - m) ** 2 for x in bm) / (nb - 1) / nb)
res = {}
for name in ('ILN', 'GBM'):
    tot = var = 0.0; rows = []
    for i in range(8):
        v = lnl(f'runs_bf/{name}_b{i}')
        if len(v) < 100: print(f'{name}_b{i}: only {len(v)} samples'); break
        m, se = mean_se(v); tot += W[i] / 2 * m; var += (W[i] / 2 * se) ** 2; rows.append((i, len(v), m, se))
    else:
        res[name] = (tot, math.sqrt(var))
        for i, n, m, se in rows: print(f'{name} beta-step {i}: n={n} E[lnL]={m:.2f} (SE {se:.2f})')
for k, (z, se) in res.items(): print(f'log marginal likelihood {k}: {z:.2f} (SE {se:.2f})')
if len(res) == 2:
    d = res['ILN'][0] - res['GBM'][0]; se = math.hypot(res['ILN'][1], res['GBM'][1])
    p = 1 / (1 + math.exp(-d)) if abs(d) < 700 else (1.0 if d > 0 else 0.0)
    print(f'log Bayes factor ILN vs GBM = {d:.2f} (SE {se:.2f}); posterior probability ILN = {p:.4f}, GBM = {1-p:.4f}')
    print('preferred clock model:', 'independent rates (ILN)' if d > 0 else 'autocorrelated rates (GBM)')
