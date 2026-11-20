"""Numerical integration, independent of R's digamma/trigamma calculation."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm, t
from scipy.integrate import quad

root=Path(__file__).resolve().parents[1]
data=pd.read_csv(root/'results/log_proxy_mechanism/analytic_constants.csv')
ans=[]
for r in data.itertuples():
    scale=1.0 if np.isinf(r.nu) else np.sqrt((r.nu-2)/r.nu)
    def logp(x, ratio=1.0):
        s=scale*np.sqrt(ratio)
        return (norm.logpdf(x/s) if np.isinf(r.nu) else t.logpdf(x/s,r.nu))-np.log(s)
    def log_square_density(u):
        x=np.exp(u/2)
        return np.exp(logp(x)+u/2)
    # Log-square tails beyond +/-80 contribute less than the stated tolerance
    # for these log moments, even for nu=2.1. This truncation is not used for E Z².
    m=quad(lambda u:u*log_square_density(u),-80,80,epsabs=1e-10)[0]
    s2=quad(lambda u:(u-m)**2*log_square_density(u),-80,80,epsabs=1e-10)[0]
    c=np.exp(m+s2/2)
    kl=quad(lambda x:np.exp(logp(x))*(logp(x)-logp(x,c)),
            -np.inf,np.inf,epsabs=1e-10,limit=1000)[0]
    error=max(abs(m-r.mean_log_z2),abs(s2-r.var_log_z2),
              abs(c-r.lognormal_variance_ratio),abs(kl-r.expected_logscore_penalty))
    ans.append(dict(distribution=r.distribution,maximum_error=error,
                    tolerance=1e-7,passed=error<1e-7))
out=pd.DataFrame(ans)
out.to_csv(root/'results/log_proxy_mechanism/independent_checks.csv',index=False)
print(out.to_string(index=False))
assert out.passed.all()
print(f'{len(out)} independent analytic-mechanism checks PASS')
