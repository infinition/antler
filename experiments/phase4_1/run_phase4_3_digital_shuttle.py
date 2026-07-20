"""ANTLER Phase 4.3: discrete site-to-site shuttle.

The Gaussian well is replaced by a compact two-site cross-fade.  At every
handoff the total trap depth is constant.  The exchange path has exactly one
JW-string crossing; the round trip has a forward/backward pair.
"""
from __future__ import annotations
import argparse,json,time,sys
from dataclasses import dataclass,asdict
from pathlib import Path
import numpy as np
from scipy.linalg import polar,sqrtm
from scipy.sparse import csr_matrix,diags
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from run_phase4_1_logical_gate import build_occ,bare_index
from run_phase4_1_logical_gate_strang import exact_logical_frame
from antler.basis import build_basis
from antler.phase1 import hop_table

@dataclass
class DConfig:
 L:int=14;N:int=2;J1:float=.4;J2:float=1.;JPERP:float=.1
 DEPTH:float=-4.;R_LOOP:int=4;T_TOTAL:float=20000.

def sin2(s): return float(np.sin(.5*np.pi*np.clip(s,0,1))**2)
def add_discrete_trap(mu,leg,x,depth,L):
 x=float(np.clip(x,0,L-1));i=int(np.floor(x));s=x-i
 if i>=L-1:mu[2*(L-1)+leg]+=depth;return
 q=sin2(s);mu[2*i+leg]+=depth*(1-q);mu[2*(i+1)+leg]+=depth*q

def mu_digital(u,exchange,cfg):
 mu=np.zeros(2*cfg.L);D=cfg.DEPTH;R=cfg.R_LOOP
 # spectator cat branch at right, unchanged
 mu[-2]=D;mu[-1]=D
 if u<.4:
  x=R*(u/.4); add_discrete_trap(mu,0,x,D,cfg.L);mu[1]+=D
 elif u<.6:
  s=(u-.4)/.2
  if exchange:
   q=sin2(s)
   mu[2*R]+=D*(1-q);mu[2*R+1]+=D*q  # moving trap crosses rung R
   mu[1]+=D*(1-q);mu[0]+=D*q        # stationary particle swaps at rung 0
  else:
   mu[2*R]+=D;mu[1]+=D
 else:
  x=R*(1-(u-.6)/.4)
  if exchange:
   add_discrete_trap(mu,1,x,D,cfg.L);mu[0]+=D
  else:
   add_discrete_trap(mu,0,x,D,cfg.L);mu[1]+=D
 return mu

def no_global(U):return U*np.exp(-.5j*np.angle(np.linalg.det(U)))
def favg(U,V):return float((abs(np.trace(V.conj().T@U))**2+2)/6)
def arr(A):return {'real':A.real.tolist(),'imag':A.imag.tolist()}

def simulate(theta,exchange,cfg,dt,states,index,table,OCC):
 d=len(states);M=2*cfg.L;rows,cols,mJ,nmid=table
 amp=mJ*np.exp(1j*theta*nmid);one=csr_matrix((amp,(rows,cols)),shape=(d,d));Hhop=one+one.conj().T
 mu0=mu_digital(0,exchange,cfg);muf=mu_digital(1,exchange,cfg)
 assert np.linalg.norm(mu0-muf)<1e-12
 H0=Hhop+diags(OCC@mu0);U0,info=exact_logical_frame(H0,index,M)
 E,V=np.linalg.eigh(Hhop.toarray());nseg=int(round(cfg.T_TOTAL/dt));dt=cfg.T_TOTAL/nseg
 Uhop=(V*np.exp(-1j*E*dt))@V.conj().T;Psi=U0.copy()
 for a in range(nseg):
  u=(a+.5)/nseg;v=OCC@mu_digital(u,exchange,cfg);ph=np.exp(-.5j*dt*v)[:,None];Psi=ph*(Uhop@(ph*Psi))
 S=U0.conj().T@Psi;return {'S':S,'U':polar(S)[0],'leak':1-np.sum(abs(S)**2,axis=0),'frame':info,'nseg':nseg,'dt':dt}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--theta',type=float,default=.3);ap.add_argument('--dt',type=float,default=.25)
 ap.add_argument('--T',type=float,default=20000.);ap.add_argument('--R',type=int,default=4);ap.add_argument('--depth',type=float,default=-4.)
 ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();cfg=DConfig(DEPTH=args.depth,R_LOOP=args.R,T_TOTAL=args.T)
 M=2*cfg.L;states,index=build_basis(M,cfg.N);table=hop_table(cfg.L,cfg.J1,cfg.J2,cfg.JPERP,states,index);OCC=build_occ(states,M)
 t0=time.time();runs={}
 for th in (args.theta,-args.theta):
  for name,ex in [('rt',False),('ex',True)]:
   print('run',name,th,'depth',args.depth,'R',args.R,flush=True);r=simulate(th,ex,cfg,args.dt,states,index,table,OCC);runs[(name,th)]=r;print('leak',r['leak'],flush=True)
 Ds={};svs=[];units=[]
 for th in (args.theta,-args.theta):
  for name in ('rt','ex'):
   S=runs[(name,th)]['S'];svs+=np.linalg.svd(S,compute_uv=False).tolist();units.append(np.linalg.norm(S.conj().T@S-np.eye(2)))
  Ds[th]=no_global(runs[('rt',th)]['U'].conj().T@runs[('ex',th)]['U'])
 Q=no_global(Ds[args.theta]@Ds[-args.theta].conj().T)
 if np.linalg.norm(-Q-np.eye(2))<np.linalg.norm(Q-np.eye(2)):Q=-Q
 Uodd=polar(sqrtm(Q).astype(complex))[0];Uodd=no_global(Uodd)
 rel=float(np.angle(np.exp(1j*(np.angle(Uodd[0,0])-np.angle(Uodd[1,1])))));target=np.diag([np.exp(-1j*args.theta/2),np.exp(1j*args.theta/2)])
 metrics={'sigma_min':float(min(svs)),'sigma_max':float(max(svs)),'leak_worst':float(1-min(svs)**2),'unitarity_frob_max':float(max(units)),'odd_phase':rel,'odd_slope':rel/args.theta,'odd_offdiag_norm':float(np.linalg.norm(Uodd-np.diag(np.diag(Uodd)))),'favg_target':favg(Uodd,target)}
 payload={'config':asdict(cfg),'theta':args.theta,'dt':args.dt,'metrics':metrics,'Uodd':arr(Uodd),'runtime_s':time.time()-t0,'runs':{}}
 for (name,th),r in runs.items():payload['runs'][f'{name}_{th:+.3f}']={'S':arr(r['S']),'U':arr(r['U']),'leak':r['leak'].tolist(),'frame':r['frame'],'nseg':r['nseg'],'dt':r['dt']}
 args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(payload,indent=2));print(json.dumps(metrics,indent=2));print('saved',args.out)
if __name__=='__main__':main()
