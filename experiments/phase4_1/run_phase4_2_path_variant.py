"""ANTLER Phase 4.2: full logical-gate audit for a deformed shuttle path."""
from __future__ import annotations
import argparse, json, time, sys
from pathlib import Path
import numpy as np
from scipy.linalg import polar, sqrtm

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from run_phase4_1_logical_gate import Config, build_occ
from run_phase4_1_logical_gate_strang import simulate
from antler.basis import build_basis
from antler.phase1 import hop_table


def mat(x): return np.asarray(x['real'])+1j*np.asarray(x['imag'])
def arr(A): return {'real':A.real.tolist(),'imag':A.imag.tolist()}
def no_global(U): return U*np.exp(-.5j*np.angle(np.linalg.det(U)))
def favg(U,V): return float((abs(np.trace(V.conj().T@U))**2+2)/6)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--theta',type=float,default=.3)
    ap.add_argument('--dt',type=float,default=.25)
    ap.add_argument('--T',type=float,default=20000.)
    ap.add_argument('--R',type=float,default=4.)
    ap.add_argument('--A',type=float,default=2.6)
    ap.add_argument('--W',type=float,default=1.0)
    ap.add_argument('--delta',type=float,default=-4.0)
    ap.add_argument('--jperp',type=float,default=.1)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    cfg=Config(T_TOTAL=args.T,R_LOOP=args.R,A_WELL=args.A,W_WELL=args.W,DELTA=args.delta,JPERP=args.jperp)
    M=2*cfg.L; states,index=build_basis(M,cfg.N)
    table=hop_table(cfg.L,cfg.J1,cfg.J2,cfg.JPERP,states,index); OCC=build_occ(states,M)
    t0=time.time(); runs={}
    for th in (args.theta,-args.theta):
        for name,ex in [('rt',False),('ex',True)]:
            print(f'run {name} theta={th:+.3f} R={args.R} A={args.A} W={args.W}',flush=True)
            r=simulate(th,ex,cfg,args.dt,states,index,table,OCC); runs[(name,th)]=r
            print(' leak',r['leak'],flush=True)
    Ds={}; allsv=[]; leak_avgs=[]; unit=[]
    for th in (args.theta,-args.theta):
        for proc in ('rt','ex'):
            S=runs[(proc,th)]['S']; sv=np.linalg.svd(S,compute_uv=False)
            allsv.extend(sv.tolist()); leak_avgs.append(float(1-np.trace(S.conj().T@S).real/2))
            unit.append(float(np.linalg.norm(S.conj().T@S-np.eye(2))))
        Ds[th]=no_global(runs[('rt',th)]['U'].conj().T@runs[('ex',th)]['U'])
    Q=no_global(Ds[args.theta]@Ds[-args.theta].conj().T)
    if np.linalg.norm(-Q-np.eye(2))<np.linalg.norm(Q-np.eye(2)): Q=-Q
    Uodd=polar(sqrtm(Q).astype(np.complex128))[0]; Uodd=no_global(Uodd)
    rel=float(np.angle(np.exp(1j*(np.angle(Uodd[0,0])-np.angle(Uodd[1,1])))))
    target=np.diag([np.exp(-1j*args.theta/2),np.exp(1j*args.theta/2)])
    payload={
      'config':cfg.__dict__,'theta':args.theta,'dt':args.dt,
      'metrics':{
        'sigma_min':float(min(allsv)),'sigma_max':float(max(allsv)),
        'leak_worst':float(1-min(allsv)**2),'leak_avg_max':float(max(leak_avgs)),
        'unitarity_frob_max':float(max(unit)),
        'odd_phase':rel,'odd_slope':float(rel/args.theta),
        'odd_offdiag_norm':float(np.linalg.norm(Uodd-np.diag(np.diag(Uodd)))),
        'favg_target':favg(Uodd,target),
      },
      'Uodd':arr(Uodd),'runs':{},'runtime_s':time.time()-t0,
    }
    for (name,th),r in runs.items():
      payload['runs'][f'{name}_{th:+.3f}']={'S':arr(r['S']),'U':arr(r['U']),'leak':r['leak'].tolist(),'frame':r['frame'],'nseg':r['nseg'],'dt':r['dt']}
    args.out.parent.mkdir(parents=True,exist_ok=True); args.out.write_text(json.dumps(payload,indent=2))
    print(json.dumps(payload['metrics'],indent=2)); print('saved',args.out,'runtime',payload['runtime_s'])
if __name__=='__main__': main()
