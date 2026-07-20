from pathlib import Path
import json,csv,re
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
RDIR=ROOT/'results/phase4_3'; rows=[]
for f in sorted(RDIR.glob('*_seq*.json')):
 d=json.loads(f.read_text()); c=d['config']; m=d['metrics']
 rows.append({'label':f.stem,'kind':'digital','depth':abs(c['DEPTH']),'R':c['R_LOOP'],'T':c['T_TOTAL'],'dt':d['dt'],**m})
for f in sorted((RDIR/'disorder').glob('*.json')):
 d=json.loads(f.read_text());m=d['metrics'];rows.append({'label':f.stem,'kind':'disorder','sigma':d['sigma'],'seed':d['seed'],'depth':abs(d['depth']),**m})
keys=sorted(set().union(*(r.keys() for r in rows)))
out=RDIR/'phase4_3_summary.csv'
with out.open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=keys);w.writeheader();w.writerows(rows)
# depth plot: canonical T=20000 dt=.25 R=4 only
p=[r for r in rows if r['kind']=='digital' and r.get('R')==4 and r.get('T')==20000 and r.get('dt')==.25 and re.fullmatch(r'D\d+_seq',r['label'])]
p=sorted(p,key=lambda r:r['depth'])
plt.figure(figsize=(7,4.5));plt.plot([r['depth'] for r in p],[-r['odd_slope'] for r in p],'o-');plt.axhline(1,ls='--');plt.xlabel('|D| (trap depth)');plt.ylabel('Exchange slope magnitude');plt.title('Strong-localization convergence');plt.grid(alpha=.25);plt.tight_layout();plt.savefig(RDIR/'phase_vs_depth.png',dpi=180);plt.close()
# R plot
p=[r for r in rows if r['kind']=='digital' and r.get('depth')==4 and r.get('T')==20000 and r.get('dt')==.25 and ('R' in r['label'] or r['label']=='D40_seq')]
p=sorted(p,key=lambda r:r['R'])
plt.figure(figsize=(7,4.5));plt.plot([r['R'] for r in p],[-r['odd_slope'] for r in p],'o-');plt.axhline(1,ls='--');plt.xlabel('Shuttle distance R');plt.ylabel('Exchange slope magnitude');plt.title('Path-length invariance');plt.grid(alpha=.25);plt.tight_layout();plt.savefig(RDIR/'phase_vs_distance.png',dpi=180);plt.close()
# disorder aggregate
p=[r for r in rows if r['kind']=='disorder']; sig=sorted(set(r['sigma'] for r in p)); means=[];mins=[];maxs=[];leaks=[]
for s in sig:
 q=[r for r in p if r['sigma']==s];v=np.array([-r['odd_slope'] for r in q]);means.append(v.mean());mins.append(v.min());maxs.append(v.max());leaks.append(max(r['leak_worst'] for r in q))
plt.figure(figsize=(7,4.5));plt.errorbar(sig,means,yerr=[np.array(means)-mins,np.array(maxs)-means],fmt='o-');plt.axhline(-json.loads((RDIR/'D40_seq.json').read_text())['metrics']['odd_slope'],ls='--');plt.xlabel('Static disorder sigma / J2');plt.ylabel('Exchange slope magnitude');plt.title('Static-disorder robustness');plt.grid(alpha=.25);plt.tight_layout();plt.savefig(RDIR/'disorder_robustness.png',dpi=180);plt.close()
print('saved',out,'and plots')
