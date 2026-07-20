from __future__ import annotations
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
files=sorted((ROOT/'results/phase4_2/path_variants').glob('*.json'))
rows=[]
for f in files:
 d=json.loads(f.read_text()); c=d['config']; m=d['metrics']
 rows.append({'label':f.stem,'R':c['R_LOOP'],'A':c['A_WELL'],'W':c['W_WELL'],'T':c['T_TOTAL'],'dt':d['dt'],**m})
out=ROOT/'results/phase4_2/path_robustness.csv'; out.parent.mkdir(parents=True,exist_ok=True)
keys=['label','R','A','W','T','dt','sigma_min','sigma_max','leak_worst','unitarity_frob_max','odd_phase','odd_slope','odd_offdiag_norm','favg_target']
with out.open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=keys);w.writeheader();w.writerows([{k:r[k] for k in keys} for r in rows])
print('\n'.join(str({k:r[k] for k in keys}) for r in rows));print('saved',out)
