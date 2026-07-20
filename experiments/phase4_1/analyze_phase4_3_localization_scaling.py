from __future__ import annotations
import csv,json
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
DEPTHS=[3.5,4.0,4.5,5.0,6.0]
FILES=[ROOT/f'results/phase4_3/D{str(d).replace(".","") if d%1 else str(int(d))+"0"}_seq.json' for d in DEPTHS]
# explicit names avoid formatting ambiguity
FILES=[ROOT/'results/phase4_3/D35_seq.json',ROOT/'results/phase4_3/D40_seq.json',ROOT/'results/phase4_3/D45_seq.json',ROOT/'results/phase4_3/D50_seq.json',ROOT/'results/phase4_3/D60_seq.json']
rows=[]
for D,f in zip(DEPTHS,FILES):
 d=json.loads(f.read_text());m=d['metrics'];rows.append({'depth':D,**m,'phase_error':1-abs(m['odd_slope'])})
x=np.array([r['depth'] for r in rows]);y=np.array([r['phase_error'] for r in rows])
def pure(x,a,p):return a*x**(-p)
def offset(x,c,a,p):return c+a*x**(-p)
pp,cov=curve_fit(pure,x,y,p0=[.8,2.5]);po,covo=curve_fit(offset,x,y,p0=[0,.8,2.5],maxfev=100000)
pred=pure(x,*pp);r2=1-np.sum((y-pred)**2)/np.sum((y-y.mean())**2)
summary={'pure_power':{'a':float(pp[0]),'p':float(pp[1]),'stderr':np.sqrt(np.diag(cov)).tolist(),'r2':float(r2)},'offset_power':{'c':float(po[0]),'a':float(po[1]),'p':float(po[2]),'stderr':np.sqrt(np.diag(covo)).tolist()},'rows':rows}
outdir=ROOT/'results/phase4_3';(outdir/'localization_scaling.json').write_text(json.dumps(summary,indent=2))
with (outdir/'localization_scaling.csv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
xx=np.linspace(3.4,6.2,200)
plt.figure(figsize=(7,4.5));plt.loglog(x,y,'o',label='ED gate data');plt.loglog(xx,pure(xx,*pp),label=fr'fit $a|\Delta|^{{-p}}$, $p={pp[1]:.3f}$')
plt.xlabel(r'Trap depth $|\Delta|$');plt.ylabel(r'Quantization error $1-|\Delta\phi/\theta|$');plt.grid(True,which='both',alpha=.25);plt.legend();plt.tight_layout();plt.savefig(outdir/'localization_scaling.png',dpi=180);plt.close()
print(json.dumps(summary['pure_power'],indent=2));print(json.dumps(summary['offset_power'],indent=2))
