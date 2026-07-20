"""Consolidated audit of ANTLER Phase 4.1 logical-gate runs."""
from __future__ import annotations
import csv, json
from pathlib import Path
import numpy as np
from scipy.linalg import polar, sqrtm

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'
RUNS=[
 ('dt0.50_R4_th0.3',RESULTS/'t03dt05.json',0.3),
 ('dt0.25_R4_th0.3',RESULTS/'phase4_1_gate_strang_dt025.json',0.3),
 ('dt0.125_R4_th0.3',RESULTS/'t03dt0125.json',0.3),
 ('dt0.25_R4_th0.6',RESULTS/'t06dt025.json',0.6),
 ('dt0.25_R4_th0.9',RESULTS/'t09dt025.json',0.9),
 ('dt0.25_R6_th0.3',RESULTS/'t03_R6_dt025.json',0.3),
]

def mat(x): return np.asarray(x['real'])+1j*np.asarray(x['imag'])
def no_global(U): return U*np.exp(-.5j*np.angle(np.linalg.det(U)))
def favg(U,V): return float((abs(np.trace(V.conj().T@U))**2+2)/6)

def audit(path,theta):
 p=json.loads(path.read_text())
 allsv=[]; leak_avgs=[]; unitarity=[]
 Ds={}
 for th in (theta,-theta):
  for proc in ('rt','ex'):
   r=p['runs'][f'{proc}_{th:+.3f}']; S=mat(r['S']); sv=np.linalg.svd(S,compute_uv=False)
   allsv.extend(sv); leak_avgs.append(1-np.trace(S.conj().T@S).real/2)
   unitarity.append(np.linalg.norm(S.conj().T@S-np.eye(2)))
  ex=mat(p['runs'][f'ex_{th:+.3f}']['U']); rt=mat(p['runs'][f'rt_{th:+.3f}']['U'])
  Ds[th]=no_global(rt.conj().T@ex)
 Q=no_global(Ds[theta]@Ds[-theta].conj().T)
 # SU(2) has the central ambiguity Q ~ -Q at the projective-gate level.
 # Choose the representative continuously connected to identity.
 if np.linalg.norm(-Q-np.eye(2)) < np.linalg.norm(Q-np.eye(2)):
  Q=-Q
 Uodd=polar(sqrtm(Q).astype(np.complex128))[0]; Uodd=no_global(Uodd)
 rel=float(np.angle(np.exp(1j*(np.angle(Uodd[0,0])-np.angle(Uodd[1,1])))))
 target=np.diag([np.exp(-1j*theta/2),np.exp(1j*theta/2)])
 return {
  'theta':theta,'dt':float(p['dt']),'T':float(p['config']['T_TOTAL']),'R':float(p['config']['R_LOOP']),
  'sigma_min':float(min(allsv)),'sigma_max':float(max(allsv)),
  'leak_worst':float(1-min(allsv)**2),'leak_avg_max':float(max(leak_avgs)),
  'unitarity_frob_max':float(max(unitarity)),
  'odd_phase':rel,'odd_slope':rel/theta,
  'odd_offdiag_norm':float(np.linalg.norm(Uodd-np.diag(np.diag(Uodd)))),
  'favg_target':favg(Uodd,target),
  'Uodd_real':Uodd.real.tolist(),'Uodd_imag':Uodd.imag.tolist(),
 }

def main():
 rows=[]
 for label,path,theta in RUNS:
  if path.exists():
   r=audit(path,theta); r['label']=label; rows.append(r)
 keys=['label','theta','dt','T','R','sigma_min','sigma_max','leak_worst','leak_avg_max','unitarity_frob_max','odd_phase','odd_slope','odd_offdiag_norm','favg_target']
 out=RESULTS/'phase4_1_gate_metrics.csv'
 with out.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys); w.writeheader(); w.writerows([{k:r[k] for k in keys} for r in rows])
 (RESULTS/'phase4_1_gate_audit.json').write_text(json.dumps(rows,indent=2))
 print('\n'.join(str({k:r[k] for k in keys}) for r in rows))
 print('saved',out)
if __name__=='__main__': main()
