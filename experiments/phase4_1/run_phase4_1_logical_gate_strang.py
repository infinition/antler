"""ANTLER Phase 4.1b: converged logical-subspace gate with Strang splitting.

This replaces the coarse piecewise-exponential probe.  The hopping exponential
is exact and precomputed once per theta; the diagonal potential is applied
exactly at the midpoint of every step.  The dressed cat frame is obtained by
dense Hermitian diagonalization and is orthonormal to machine precision.
"""
from __future__ import annotations
import argparse, json, time
from dataclasses import dataclass, asdict
from pathlib import Path
import numpy as np
from scipy.linalg import polar
from scipy.sparse import csr_matrix, diags
from antler.basis import build_basis
from antler.phase1 import hop_table
from run_phase4_1_logical_gate import Config, protocol_mu, build_occ, bare_index


def exact_logical_frame(H0, index, M):
    Hd=H0.toarray() if hasattr(H0,'toarray') else np.asarray(H0)
    E,V=np.linalg.eigh(Hd)
    iLL=bare_index(index,0,1); iRR=bare_index(index,M-2,M-1)
    score=np.abs(V[iLL])**2+np.abs(V[iRR])**2
    pair=np.argsort(-score)[:2]
    W,_=np.linalg.qr(V[:,pair])
    z=np.zeros(Hd.shape[0]); z[iLL]=1; z[iRR]=-1
    Z=W.conj().T@(z[:,None]*W)
    _,Q=np.linalg.eigh(Z)
    U=W@Q[:,::-1]
    for j in range(2):
        k=int(np.argmax(np.abs(U[:,j]))); U[:,j]*=np.conj(U[k,j])/abs(U[k,j])
    return U, {'orth_error':float(np.linalg.norm(U.conj().T@U-np.eye(2))),
               'bare_LL':float(abs(U[iLL,0])**2),'bare_RR':float(abs(U[iRR,1])**2),
               'energies':E[pair].tolist()}


def remove_global(U): return U*np.exp(-.5j*np.angle(np.linalg.det(U)))
def unitary_part(S): return polar(S)[0]
def gate_fidelity(U,target):
    d=2; return float((abs(np.trace(target.conj().T@U))**2+d)/(d*(d+1)))


def simulate(theta, exchange, cfg, dt, states,index,table,OCC):
    d=len(states); M=2*cfg.L
    rows,cols,mJ,nmid=table
    amp=mJ*np.exp(1j*theta*nmid)
    one=csr_matrix((amp,(rows,cols)),shape=(d,d)); Hhop=one+one.conj().T
    mu0=protocol_mu(0,exchange,cfg); muf=protocol_mu(1,exchange,cfg)
    assert np.linalg.norm(mu0-muf)<1e-12
    H0=Hhop+diags(OCC@mu0)
    U0,info=exact_logical_frame(H0,index,M)
    E,V=np.linalg.eigh(Hhop.toarray())
    Uhop=(V*np.exp(-1j*E*dt))@V.conj().T
    nseg=int(round(cfg.T_TOTAL/dt)); dt_eff=cfg.T_TOTAL/nseg
    if abs(dt_eff-dt)>1e-12:
        Uhop=(V*np.exp(-1j*E*dt_eff))@V.conj().T
    Psi=U0.copy()
    for a in range(nseg):
        u=(a+.5)/nseg
        v=OCC@protocol_mu(u,exchange,cfg)
        ph=np.exp(-.5j*dt_eff*v)[:,None]
        Psi=ph*(Uhop@(ph*Psi))
    S=U0.conj().T@Psi
    leak=1-np.sum(np.abs(S)**2,axis=0)
    return {'S':S,'U':unitary_part(S),'leak':leak,'frame':info,'nseg':nseg,'dt':dt_eff}


def pair_analysis(ex,rt,theta):
    D=remove_global(rt['U'].conj().T@ex['U'])
    rel=float(np.angle(np.exp(1j*(np.angle(D[0,0])-np.angle(D[1,1])))))
    off=float(np.linalg.norm(D-np.diag(np.diag(D))))
    target_minus=remove_global(np.diag([np.exp(-1j*theta),1]))
    target_plus=remove_global(np.diag([np.exp(1j*theta),1]))
    return D,rel,off,{'minus':gate_fidelity(D,target_minus),'plus':gate_fidelity(D,target_plus)}


def arr(A): return {'real':A.real.tolist(),'imag':A.imag.tolist()}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--theta',type=float,default=.3)
    ap.add_argument('--dt',type=float,default=.25); ap.add_argument('--T',type=float,default=20000.)
    ap.add_argument('--R',type=float,default=4.); ap.add_argument('--out',type=Path,default=Path('results/phase4_1_gate_strang.json'))
    args=ap.parse_args(); cfg=Config(T_TOTAL=args.T,R_LOOP=args.R)
    M=2*cfg.L; states,index=build_basis(M,cfg.N); table=hop_table(cfg.L,cfg.J1,cfg.J2,cfg.JPERP,states,index); OCC=build_occ(states,M)
    t0=time.time(); runs={}
    for th in (args.theta,-args.theta):
        for name,ex in [('rt',False),('ex',True)]:
            print(f'run {name} theta={th:+.3f} dt={args.dt}',flush=True)
            r=simulate(th,ex,cfg,args.dt,states,index,table,OCC); runs[(name,th)]=r
            print(' leak',r['leak'],'orth',r['frame']['orth_error'],flush=True)
    ana={}
    for th in (args.theta,-args.theta):
        D,rel,off,fids=pair_analysis(runs[('ex',th)],runs[('rt',th)],th)
        ana[str(th)]={'D':arr(D),'relative_phase':rel,'offdiag':off,'fidelities':fids}
        print(f'theta {th:+.3f}: rel={rel:+.8f} off={off:.3e} fids={fids}',flush=True)
    odd=.5*np.angle(np.exp(1j*(ana[str(args.theta)]['relative_phase']-ana[str(-args.theta)]['relative_phase'])))
    print(f'odd phase={odd:+.8f}; slope={odd/args.theta:+.8f}',flush=True)
    payload={'config':asdict(cfg),'dt':args.dt,'runs':{},'analysis':ana,'odd_phase':float(odd),'odd_slope':float(odd/args.theta),'runtime_s':time.time()-t0}
    for (name,th),r in runs.items(): payload['runs'][f'{name}_{th:+.3f}']={'S':arr(r['S']),'U':arr(r['U']),'leak':r['leak'].tolist(),'frame':r['frame'],'nseg':r['nseg'],'dt':r['dt']}
    args.out.parent.mkdir(parents=True,exist_ok=True); args.out.write_text(json.dumps(payload,indent=2)); print('saved',args.out,'runtime',payload['runtime_s'])
if __name__=='__main__': main()
