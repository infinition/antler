"""ANTLER Phase 4.4: dynamic JW-string response at theta=0.

For each protocol P in {exchange, round-trip}, propagate the full logical frame
at theta=0 and integrate the Hermitian response operator dH/dtheta|_0.
For a cyclic adiabatic evolution the derivative of the dynamical phase is
  d phi_dyn / dtheta = - integral <dH/dtheta> dt.
The differential logical-Z slope predicted by this response is compared with
the finite-theta odd phase measured by Phase 4.2.
"""
from __future__ import annotations
import argparse, json, sys, time
from dataclasses import asdict
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix, diags

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from run_phase4_1_logical_gate import Config, protocol_mu, build_occ
from run_phase4_1_logical_gate_strang import exact_logical_frame
from antler.basis import build_basis
from antler.phase1 import hop_table


def simulate_response(exchange,cfg,dt,states,index,table,OCC,stride=20):
    d=len(states);M=2*cfg.L;rows,cols,mJ,nmid=table
    one=csr_matrix((mJ,(rows,cols)),shape=(d,d));Hhop=one+one.conj().T
    done=csr_matrix((1j*mJ*nmid,(rows,cols)),shape=(d,d));Htheta=(done+done.conj().T).tocsr()
    mu0=protocol_mu(0,exchange,cfg); assert np.linalg.norm(mu0-protocol_mu(1,exchange,cfg))<1e-12
    H0=Hhop+diags(OCC@mu0);U0,info=exact_logical_frame(H0,index,M)
    E,V=np.linalg.eigh(Hhop.toarray());nseg=int(round(cfg.T_TOTAL/dt));dt=cfg.T_TOTAL/nseg
    Uhop=(V*np.exp(-1j*E*dt))@V.conj().T
    Psi=U0.copy(); integral=np.zeros(2); abs_integral=np.zeros(2)
    exp_last=np.real(np.diag(Psi.conj().T@(Htheta@Psi))); t_last=0.0
    checkpoints=[]
    marks=[.1,.38,.58,.88,1.0]; imark=0
    for a in range(nseg):
        u=(a+.5)/nseg; v=OCC@protocol_mu(u,exchange,cfg); ph=np.exp(-.5j*dt*v)[:,None]
        Psi=ph*(Uhop@(ph*Psi))
        sample=((a+1)%stride==0) or (a==nseg-1)
        if sample:
            t_now=(a+1)*dt
            exp_now=np.real(np.diag(Psi.conj().T@(Htheta@Psi)))
            h=t_now-t_last
            integral += .5*h*(exp_last+exp_now)
            abs_integral += .5*h*(np.abs(exp_last)+np.abs(exp_now))
            exp_last=exp_now; t_last=t_now
            while imark<len(marks) and (a+1)/nseg>=marks[imark]:
                checkpoints.append({'u':marks[imark], 'integral':integral.tolist(),'expect':exp_now.tolist()});imark+=1
    S=U0.conj().T@Psi
    return {'integral':integral,'abs_integral':abs_integral,'final_leak':(1-np.sum(abs(S)**2,axis=0)),
            'S':S,'frame':info,'checkpoints':checkpoints,'nseg':nseg,'dt':dt,'stride':stride}

def arr(A):return {'real':A.real.tolist(),'imag':A.imag.tolist()}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dt',type=float,default=.5);ap.add_argument('--T',type=float,default=20000.)
    ap.add_argument('--R',type=float,default=4.);ap.add_argument('--A',type=float,default=2.6);ap.add_argument('--W',type=float,default=1.)
    ap.add_argument('--stride',type=int,default=20);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    cfg=Config(T_TOTAL=args.T,R_LOOP=args.R,A_WELL=args.A,W_WELL=args.W)
    M=2*cfg.L;states,index=build_basis(M,cfg.N);table=hop_table(cfg.L,cfg.J1,cfg.J2,cfg.JPERP,states,index);OCC=build_occ(states,M)
    t0=time.time(); runs={}
    for name,ex in [('rt',False),('ex',True)]:
        print('response',name,'R/A/W',args.R,args.A,args.W,flush=True)
        runs[name]=simulate_response(ex,cfg,args.dt,states,index,table,OCC,args.stride)
        print(' integral',runs[name]['integral'],'leak',runs[name]['final_leak'],flush=True)
    # logical relative phase = phi_L - phi_R. Relative ex-vs-rt derivative:
    dI=(runs['ex']['integral'][0]-runs['ex']['integral'][1])-(runs['rt']['integral'][0]-runs['rt']['integral'][1])
    dyn_slope=-float(dI)
    payload={'config':asdict(cfg),'dt':args.dt,'predicted_dynamic_odd_slope':dyn_slope,'runs':{},'runtime_s':time.time()-t0}
    for name,r in runs.items():
        payload['runs'][name]={'integral':r['integral'].tolist(),'abs_integral':r['abs_integral'].tolist(),'final_leak':r['final_leak'].tolist(),'S':arr(r['S']),'frame':r['frame'],'checkpoints':r['checkpoints'],'nseg':r['nseg'],'dt':r['dt'],'stride':r['stride']}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(payload,indent=2))
    print('predicted dynamic odd slope',dyn_slope);print('saved',args.out,'runtime',payload['runtime_s'])
if __name__=='__main__':main()
