"""Independent SciPy checks of R forecasting arithmetic; no model refitting."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t as student
from scipy.optimize import brentq
from scipy.special import logsumexp
from scipy.integrate import quad

root = Path(__file__).resolve().parents[1]
out = root / 'results/main_sample'
checks = []
def check(name, value, tol):
    checks.append(dict(check=name, maximum_error=float(value), tolerance=tol,
                       passed=bool(np.isfinite(value) and value <= tol)))

for folder in sorted(out.iterdir()):
    if not (folder / 'forecast_bank.csv').exists():
        continue
    ticker = folder.name
    b = pd.read_csv(folder / 'forecast_bank.csv')
    a = b[b.model == 'GARCH'].sort_values('target_date')
    c = b[b.model == 'GJR'].sort_values('target_date')
    assert np.array_equal(a.target_date, c.target_date)
    y = a.ret.to_numpy()
    var = np.column_stack([a.sigma2, c.sigma2])
    nu = np.column_stack([a.df, c.df])
    scale = np.sqrt(var * (nu - 2) / nu)
    logf = student.logpdf(y[:, None] / scale, nu) - np.log(scale)
    states = np.column_stack([np.linspace(0, 1, 21), np.linspace(1, 0, 21)])
    candidates, scores = [], []
    validation = a.target_date.to_numpy() < '2022-01-01'
    for hazard in [0, .002, .01, .05]:
        prob = np.full(21, 1/21)
        weights, ls = [], []
        for row in logf:
            weights.append(prob @ states)
            mx = row.max()
            emission = np.log(states @ np.exp(row-mx)) + mx
            ls.append(logsumexp(np.log(prob) + emission))
            logpost = np.log(prob) + emission
            prob = (1-hazard)*np.exp(logpost-logsumexp(logpost)) + hazard/21
        candidates.append(np.array(weights))
        scores.append(np.mean(np.array(ls)[validation]))
    w = candidates[int(np.argmax(scores))]
    rw = pd.read_csv(folder / 'pair_forecast_weights.csv')
    rw = rw[rw.method == 'Pool_filter'].sort_values('date')
    check(ticker+'_filter_weights', np.max(np.abs(w-rw[['GARCH','GJR']].to_numpy())), 1e-10)
    chosen = pd.read_csv(folder / 'validation_selection.csv')
    chosen = chosen[chosen.family == 'Pool_filter']
    check(ticker+'_validation_scores', np.max(np.abs(chosen.validation_score-scores)), 1e-10)
    keep = ~validation
    r = pd.read_csv(root / 'results/independent_audit' / (ticker+'_Pool_filter.csv'))
    nll = -logsumexp(np.log(w[keep])+logf[keep],axis=1)
    pit = np.sum(w[keep]*student.cdf(y[keep,None]/scale[keep],nu[keep]),axis=1)
    check(ticker+'_negative_log_score',np.max(np.abs(nll-r.NLL)),1e-10)
    check(ticker+'_PIT',np.max(np.abs(pit-r.PIT)),1e-10)
    idx = np.flatnonzero(keep)
    qerr, eerr, ferr = [], [], []
    for ri in np.linspace(0,len(idx)-1,16).astype(int):
        ti = idx[ri]; ss, nn, ww = scale[ti], nu[ti], w[ti]
        cq = ss*student.ppf(.025,nn)
        func = lambda q: float(ww @ student.cdf(q/ss,nn)) - .025
        q = brentq(func,cq.min()-1e-10,cq.max()+1e-10,xtol=1e-12)
        density = lambda x: float(ww @ (student.pdf(x/ss,nn)/ss))
        es = quad(lambda x:x*density(x),-np.inf,q,epsabs=1e-10)[0]/.025
        fz = -float(y[ti]<=q)*(q-y[ti])/(.025*es)+q/es+np.log(-es)-1
        qerr.append(abs(q-r.VaR_0p025.iloc[ri]))
        eerr.append(abs(es-r.ES_0p025.iloc[ri]))
        ferr.append(abs(fz-r.FZ_0p025.iloc[ri]))
    check(ticker+'_VaR_root',max(qerr),1e-8)
    check(ticker+'_ES_quadrature',max(eerr),1e-7)
    check(ticker+'_FZ0',max(ferr),1e-8)
    ref = pd.read_csv(root / 'results/independent_audit' / (ticker+'_Static_pair.csv'))
    d = r.NLL.to_numpy()-ref.NLL.to_numpy(); n=len(d); lag=int(n**(1/3));z=d-d.mean()
    lrv=np.dot(z,z)/n
    for k in range(1,lag+1): lrv+=2*(1-k/(lag+1))*np.dot(z[k:],z[:-k])/n
    cmp=pd.read_csv(folder/'comparisons.csv')
    cmp=cmp[(cmp.candidate=='Pool_filter')&(cmp.comparator=='Static_pair')&(cmp.metric=='NLL')].iloc[0]
    check(ticker+'_HAC_difference',abs(d.mean()-cmp.difference),1e-11)
    check(ticker+'_HAC_SE',abs(np.sqrt(lrv/n)-cmp.se),1e-11)

result=pd.DataFrame(checks)
result.to_csv(root/'results/independent_audit/python_checks.csv',index=False)
print(result.to_string(index=False))
assert result.passed.all(), 'Independent numerical audit failed'
print(f'{len(result)} independent SciPy checks PASS')
