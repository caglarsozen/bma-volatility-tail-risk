"""Independent SciPy replay of the new scale experiments on saved forecast banks."""
from pathlib import Path
import hashlib

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.special import logsumexp
from scipy.stats import norm, t

root = Path(__file__).resolve().parents[1]
base = root / 'results/main_sample'
results = root / 'results/diagnostics'
checks = []

def check(name, actual, expected, tol):
    error = np.max(np.abs(np.asarray(actual) - np.asarray(expected)))
    checks.append({'check': name, 'maximum_error': error,
                   'tolerance': tol, 'pass': bool(np.isfinite(error) and error <= tol)})

def read_result(file):
    return pd.read_csv(results / file)

def loss_cells(y, v_smearing, v_normal, df, test):
    """Aggregate last likelihoods explicitly rather than using R's recursion."""
    losses = {}
    times = np.arange(len(y))
    powers = .97 ** np.arange(len(y))
    for label, var in [('Normal', v_normal), ('Smearing', v_smearing)]:
        scales = np.sqrt(var * (df-2)/df)
        logf = t.logpdf(y[:,None]/scales, df) - np.log(scales)
        # The lag starts at one, so the contemporaneous density is excluded.
        logs = np.column_stack([
            np.r_[0, np.convolve(col, powers)[:len(y)-1]]
            for col in logf.T])
        dma = np.exp(logs - logsumexp(logs, axis=1, keepdims=True))
        ws = [('DMA', dma), ('Equal', np.ones_like(logf)/logf.shape[1])]
        for kind, w in ws:
            mv = np.sum(w * var, axis=1)
            with np.errstate(divide='ignore'):
                nll = -logsumexp(np.log(w) + logf, axis=1)
            # Independent bracketed quantiles with a fixed number of bisections.
            a = .025
            qv = scales * t.ppf(a, df)
            lo = qv.min(axis=1)
            hi = qv.max(axis=1)
            for _ in range(45):
                mid = (lo+hi)/2
                low_side = np.sum(w*t.cdf(mid[:,None]/scales,df),axis=1) < a
                lo = np.where(low_side, mid, lo)
                hi = np.where(low_side, hi, mid)
            q = (lo+hi)/2
            z = q[:,None]/scales
            es = -np.sum(w*scales*(df+z**2)/(df-1)*t.pdf(z,df),axis=1)/a
            fz = -(y<=q).astype(float)*(q-y)/(a*es) + q/es + np.log(-es) - 1
            losses[(label,kind)] = np.column_stack([nll,np.log(mv)+y*y/mv,fz])[test]
    return losses

def bank_arrays(bank, models, split):
    cols = [bank[bank.model.eq(m)].sort_values('target_date') for m in models]
    ref = cols[0]
    assert all(np.array_equal(b.target_date,ref.target_date) for b in cols)
    assert all(np.array_equal(b.ret,ref.ret) for b in cols)
    v = np.column_stack([b.sigma2 for b in cols])
    old = np.column_stack([b.sigma2_lognormal.fillna(b.sigma2) for b in cols])
    df = np.column_stack([b.df for b in cols])
    return ref.ret.to_numpy(), v, old, df, ref.target_date.to_numpy()>=split

def compare_report(name,cells,rows):
    keys = [('Normal','Equal'),('Normal','DMA'),('Smearing','Equal'),('Smearing','DMA')]
    columns = ['old_equal','old_DMA','corrected_equal','corrected_DMA']
    r = rows.set_index('metric').loc[['NLL','QLIKE','FZ_0p025']]
    for key, col in zip(keys, columns):
        check(name+' '+col, cells[key].mean(axis=0), r[col], 2e-9)
    a,b,c,d = (cells[k] for k in keys)
    ix = a-b-c+d
    check(name+' interaction', ix.mean(axis=0), r.difference, 3e-9)
    n = len(ix)
    lag = int(n**(1/3))
    e = ix - ix.mean(axis=0)
    omega = np.sum(e*e,axis=0)
    for l in range(1,lag+1):
        omega += 2*(1-l/(lag+1))*np.sum(e[l:]*e[:-l],axis=0)
    check(name+' interaction HAC standard errors',np.sqrt(omega)/n,r.se,2e-9)
    return ix

models = ['GARCH','GJR','HAR','Range_HAR','HAR66','Log_AR1']
m = read_result('scale_weight/main_factorial.csv')
for ticker in m.ticker.unique():
    bank = pd.read_csv(base / ticker / 'forecast_bank.csv')
    cells = loss_cells(*bank_arrays(bank,models,'2022-01-01'))
    compare_report(ticker,cells,m[m.ticker.eq(ticker)])

e = read_result('scale_weight/external_factorial.csv')
bank = read_result('realized_benchmark/forecast_bank.csv')
external = {}
for proxy,model,filename in [('Realized kernel','HAR_RK','external_HAR_RK.csv'),
                             ('Squared open-to-close return','HAR','external_HAR_R2.csv')]:
    common = bank[bank.model.ne('HAR_RK')]
    other = read_result('scale_weight/'+filename)
    # CSV roundoff may differ in the final digit on separately written inputs.
    ref_y = common[common.model.eq('GARCH')].set_index('target_date').ret
    check(proxy+' target equality',other.ret,other.target_date.map(ref_y),1e-12)
    other['ret'] = other.target_date.map(ref_y)
    cells = loss_cells(*bank_arrays(pd.concat([common,other]),['GARCH','GJR','RealGARCH',model],'2005-01-01'))
    external[proxy] = compare_report(proxy,cells,e[e.proxy.eq(proxy)])
paired = read_result('scale_weight/external_paired_proxy_contrast.csv')
check('matched dates proxy contrast',(external['Squared open-to-close return'] -
    external['Realized kernel']).mean(axis=0),paired.difference,3e-9)

g = read_result('proxy_granularity/analytic_granularity.csv')
for r in g.itertuples():
    def laplace(s):
        if s < 1e-3:
            numerator = s**3/6 - s**4/24 + s**5/120 - s**6/720
        else:
            numerator = -np.expm1(-s)-s+s*s/2
        return np.exp(-r.d*s/2) * numerator / (s*(-np.expm1(-s)))
    check(f'd={r.d} analytical Laplace constant',quad(laplace,0,np.inf,epsabs=1e-12)[0],r.log_c,1e-10)
    kl = quad(lambda x:-norm.pdf(x)*np.log1p((2/3)*np.expm1(
        -.5*np.log(r.c)+.5*x*x*(1-1/r.c))),-12,12,epsabs=1e-11)[0]
    check(f'd={r.d} Gaussian mixture KL',kl,r.NLL_deficit,5e-10)

manifest = root/'docs/windows_run_input_manifest.csv'
if manifest.is_file():
    mf = pd.read_csv(manifest)
    for row in mf[mf.file.str.startswith('data/raw/')].itertuples():
        got = hashlib.md5((root/row.file).read_bytes()).hexdigest()
        check(row.file+' original execution MD5',int(got!=row.md5),0,0)

out = pd.DataFrame(checks)
directory = results/'verification'
directory.mkdir(exist_ok=True)
out.to_csv(directory/'independent_revision_checks.csv',index=False)
assert out['pass'].all(),out[~out['pass']].to_string(index=False)
print(f"{len(out)} independent revision checks passed")
