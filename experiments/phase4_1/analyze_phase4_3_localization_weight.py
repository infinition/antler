from __future__ import annotations
import csv,json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
labels=['D35','D40','D45','D50','D60'];rows=[]
for lab in labels:
 d=json.loads((ROOT/f'results/phase4_3/{lab}_seq.json').read_text())
 m=d['metrics']; f=d['runs']['rt_+0.300']['frame'];P=.5*(f['bare_LL']+f['bare_RR'])
 rows.append({'label':lab,'depth':abs(d['config']['DEPTH']),'bare_weight':P,'delocalized_weight':1-P,'odd_slope_abs':abs(m['odd_slope']),'phase_error':1-abs(m['odd_slope'])})
x=np.array([r['delocalized_weight'] for r in rows]);y=np.array([r['phase_error'] for r in rows]);coef=np.polyfit(x,y,1);corr=float(np.corrcoef(x,y)[0,1])
summary={'pearson_r':corr,'linear_slope':float(coef[0]),'linear_intercept':float(coef[1]),'rows':rows}
out=ROOT/'results/phase4_3';(out/'localization_weight_correlation.json').write_text(json.dumps(summary,indent=2))
with (out/'localization_weight_correlation.csv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
xx=np.linspace(0,x.max()*1.05,200)
plt.figure(figsize=(7,4.5));plt.plot(x,y,'o',label='digital shuttle');plt.plot(xx,coef[0]*xx+coef[1],label=f'linear fit, r={corr:.7f}')
plt.xlabel('Delocalized logical weight $1-P_{bare}$');plt.ylabel('Phase quantization error $1-|\\Delta\\phi/\\theta|$');plt.grid(True,alpha=.25);plt.legend();plt.tight_layout();plt.savefig(out/'localization_weight_correlation.png',dpi=180);plt.close();print(json.dumps(summary,indent=2))
