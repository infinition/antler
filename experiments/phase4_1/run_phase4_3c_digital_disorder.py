"""ANTLER Phase 4.3c: sequential digital exchange under static on-site disorder."""
from __future__ import annotations
import argparse,json,time,sys
from pathlib import Path
import numpy as np
from scipy.linalg import polar,sqrtm
from scipy.sparse import csr_matrix,diags
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from run_phase4_3b_digital_sequential import DConfig,mu_digital
from run_phase4_1_logical_gate import build_occ
from run_phase4_1_logical_gate_strang import exact_logical_frame
from antler.basis import build_basis
from antler.phase1 import hop_table

def ng(U):return U*np.exp(-.5j*np.angle(np.linalg.det(U)))
def favg(U,V):return float((abs(np.trace(V.conj().T@U))**2+2)/6)
def simulate(th,ex,cfg,dt,states,index,table,OCC,eps):
 d=len(states);M=2*cfg.L;rows,cols,mJ,nmid=table;amp=mJ*np.exp(1j*th*nmid);one=csr_matrix((amp,(rows,cols)),shape=(d,d));Hhop=one+one.conj().T
 H0=Hhop+diags(OCC@(mu_digital(0,ex,cfg)+eps));U0,_=exact_logical_frame(H0,index,M);E,V=np.linalg.eigh(Hhop.toarray());n=int(round(cfg.T_TOTAL/dt));dt=cfg.T_TOTAL/n;Uh=(V*np.exp(-1j*E*dt))@V.conj().T;P=U0.copy()
 for a in range(n):
  u=(a+.5)/n;v=OCC@(mu_digital(u,ex,cfg)+eps);ph=np.exp(-.5j*dt*v)[:,None];P=ph*(Uh@(ph*P))
 S=U0.conj().T@P;return S,polar(S)[0]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--sigma',type=float,required=True);ap.add_argument('--seed',type=int,required=True);ap.add_argument('--depth',type=float,default=-4.);ap.add_argument('--theta',type=float,default=.3);ap.add_argument('--dt',type=float,default=.25);ap.add_argument('--T',type=float,default=20000.);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();cfg=DConfig(DEPTH=a.depth,T_TOTAL=a.T)
 rng=np.random.default_rng(a.seed);eps=rng.normal(0,a.sigma,2*cfg.L);eps-=eps.mean();M=2*cfg.L;states,index=build_basis(M,cfg.N);table=hop_table(cfg.L,cfg.J1,cfg.J2,cfg.JPERP,states,index);OCC=build_occ(states,M);runs={};t0=time.time()
 for th in (a.theta,-a.theta):
  for name,ex in [('rt',False),('ex',True)]:print('run',name,th,'sig',a.sigma,'seed',a.seed,flush=True);runs[(name,th)]=simulate(th,ex,cfg,a.dt,states,index,table,OCC,eps)
 sv=[];un=[];D={}
 for th in (a.theta,-a.theta):
  for name in ('rt','ex'):
   S,U=runs[(name,th)];sv+=np.linalg.svd(S,compute_uv=False).tolist();un.append(np.linalg.norm(S.conj().T@S-np.eye(2)))
  D[th]=ng(runs[('rt',th)][1].conj().T@runs[('ex',th)][1])
 Q=ng(D[a.theta]@D[-a.theta].conj().T)
 if np.linalg.norm(-Q-np.eye(2))<np.linalg.norm(Q-np.eye(2)):Q=-Q
 Uo=ng(polar(sqrtm(Q).astype(complex))[0]);rel=float(np.angle(np.exp(1j*(np.angle(Uo[0,0])-np.angle(Uo[1,1])))));target=np.diag([np.exp(-1j*a.theta/2),np.exp(1j*a.theta/2)])
 m={'sigma_min':float(min(sv)),'leak_worst':float(1-min(sv)**2),'unitarity_frob_max':float(max(un)),'odd_phase':rel,'odd_slope':rel/a.theta,'odd_offdiag_norm':float(np.linalg.norm(Uo-np.diag(np.diag(Uo)))),'favg_target':favg(Uo,target)}
 out={'sigma':a.sigma,'seed':a.seed,'depth':a.depth,'epsilon':eps.tolist(),'metrics':m,'runtime_s':time.time()-t0};a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2));print(json.dumps(m,indent=2))
if __name__=='__main__':main()
