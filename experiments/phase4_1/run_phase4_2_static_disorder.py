"""ANTLER Phase 4.2: logical gate under a fixed on-site disorder realization.
The same disorder vector is used for exchange/round-trip and +/-theta.
"""
from __future__ import annotations
import argparse,json,time,sys
from pathlib import Path
import numpy as np
from scipy.linalg import polar,sqrtm
from scipy.sparse import csr_matrix,diags
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from run_phase4_1_logical_gate import Config,protocol_mu,build_occ
from run_phase4_1_logical_gate_strang import exact_logical_frame
from antler.basis import build_basis
from antler.phase1 import hop_table

def no_global(U):return U*np.exp(-.5j*np.angle(np.linalg.det(U)))
def favg(U,V):return float((abs(np.trace(V.conj().T@U))**2+2)/6)
def arr(A):return {'real':A.real.tolist(),'imag':A.imag.tolist()}

def simulate(theta,exchange,cfg,dt,states,index,table,OCC,eps):
 d=len(states);M=2*cfg.L; rows,cols,mJ,nmid=table
 amp=mJ*np.exp(1j*theta*nmid); one=csr_matrix((amp,(rows,cols)),shape=(d,d)); Hhop=one+one.conj().T
 mu0=protocol_mu(0,exchange,cfg)+eps; muf=protocol_mu(1,exchange,cfg)+eps
 assert np.linalg.norm(mu0-muf)<1e-12
 H0=Hhop+diags(OCC@mu0); U0,info=exact_logical_frame(H0,index,M)
 E,V=np.linalg.eigh(Hhop.toarray());nseg=int(round(cfg.T_TOTAL/dt));dt=cfg.T_TOTAL/nseg
 Uhop=(V*np.exp(-1j*E*dt))@V.conj().T;Psi=U0.copy()
 for a in range(nseg):
  u=(a+.5)/nseg;v=OCC@(protocol_mu(u,exchange,cfg)+eps);ph=np.exp(-.5j*dt*v)[:,None];Psi=ph*(Uhop@(ph*Psi))
 S=U0.conj().T@Psi; leak=1-np.sum(abs(S)**2,axis=0);U=polar(S)[0]
 return {'S':S,'U':U,'leak':leak,'frame':info}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--sigma',type=float,required=True);ap.add_argument('--seed',type=int,required=True)
 ap.add_argument('--theta',type=float,default=.3);ap.add_argument('--dt',type=float,default=.25);ap.add_argument('--T',type=float,default=20000.)
 ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();cfg=Config(T_TOTAL=args.T,R_LOOP=4.)
 rng=np.random.default_rng(args.seed);eps=rng.normal(0,args.sigma,2*cfg.L);eps-=eps.mean()
 M=2*cfg.L;states,index=build_basis(M,cfg.N);table=hop_table(cfg.L,cfg.J1,cfg.J2,cfg.JPERP,states,index);OCC=build_occ(states,M)
 t0=time.time();runs={}
 for th in (args.theta,-args.theta):
  for name,ex in [('rt',False),('ex',True)]:
   print('run',name,th,'sigma',args.sigma,'seed',args.seed,flush=True);runs[(name,th)]=simulate(th,ex,cfg,args.dt,states,index,table,OCC,eps)
 Ds={};svs=[];units=[];lav=[]
 for th in (args.theta,-args.theta):
  for name in ('rt','ex'):
   S=runs[(name,th)]['S'];s=np.linalg.svd(S,compute_uv=False);svs+=s.tolist();units.append(np.linalg.norm(S.conj().T@S-np.eye(2)));lav.append(1-np.trace(S.conj().T@S).real/2)
  Ds[th]=no_global(runs[('rt',th)]['U'].conj().T@runs[('ex',th)]['U'])
 Q=no_global(Ds[args.theta]@Ds[-args.theta].conj().T)
 if np.linalg.norm(-Q-np.eye(2))<np.linalg.norm(Q-np.eye(2)):Q=-Q
 Uodd=polar(sqrtm(Q).astype(complex))[0];Uodd=no_global(Uodd)
 rel=float(np.angle(np.exp(1j*(np.angle(Uodd[0,0])-np.angle(Uodd[1,1])))));target=np.diag([np.exp(-1j*args.theta/2),np.exp(1j*args.theta/2)])
 m={'sigma_min':float(min(svs)),'leak_worst':float(1-min(svs)**2),'unitarity_frob_max':float(max(units)),'odd_phase':rel,'odd_slope':rel/args.theta,'odd_offdiag_norm':float(np.linalg.norm(Uodd-np.diag(np.diag(Uodd)))),'favg_target':favg(Uodd,target)}
 payload={'sigma':args.sigma,'seed':args.seed,'epsilon':eps.tolist(),'config':cfg.__dict__,'metrics':m,'Uodd':arr(Uodd),'runtime_s':time.time()-t0}
 args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(payload,indent=2));print(json.dumps(m,indent=2))
if __name__=='__main__':main()
