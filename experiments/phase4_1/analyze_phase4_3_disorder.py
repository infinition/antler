from __future__ import annotations
import csv,json,glob,os,re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];files=glob.glob(str(ROOT/'results/phase4_3/disorder/*.json'));rows=[]
for f in files:
 d=json.load(open(f));m=d['metrics'];rows.append({'label':Path(f).stem,'sigma':d['sigma'],'seed':d['seed'],**m})
rows.sort(key=lambda r:(r['sigma'],r['seed']))
out=ROOT/'results/phase4_3/disorder_summary.csv'
keys=['label','sigma','seed','sigma_min','leak_worst','unitarity_frob_max','odd_phase','odd_slope','odd_offdiag_norm','favg_target']
with out.open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=keys);w.writeheader();w.writerows([{k:r[k] for k in keys} for r in rows])
summary=[]
for s in sorted(set(r['sigma'] for r in rows)):
 rr=[r for r in rows if r['sigma']==s]
 summary.append({'sigma':s,'n':len(rr),'slope_mean':float(np.mean([r['odd_slope'] for r in rr])),'slope_std':float(np.std([r['odd_slope'] for r in rr])),'leak_max':float(max(r['leak_worst'] for r in rr)),'favg_min':float(min(r['favg_target'] for r in rr))})
(ROOT/'results/phase4_3/disorder_summary.json').write_text(json.dumps({'runs':rows,'summary':summary},indent=2))
plt.figure(figsize=(7,4.5))
for seed in sorted(set(r['seed'] for r in rows)):
 rr=[r for r in rows if r['seed']==seed];plt.plot([r['sigma'] for r in rr],[r['odd_slope'] for r in rr],'o-',label=f'seed {seed}')
plt.axhline(-.9735108633262304,linestyle='--',label='clean digital')
plt.xlabel(r'On-site disorder $\sigma_\mu$');plt.ylabel(r'Odd phase slope $\Delta\phi/\theta$');plt.grid(True,alpha=.25);plt.legend();plt.tight_layout();plt.savefig(ROOT/'results/phase4_3/disorder_slope.png',dpi=180);plt.close()
print(json.dumps(summary,indent=2))
