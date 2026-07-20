"""ANTLER Phase 4.3b: sequential compact exchange shuttle.

Compared with the first digital protocol, the two rung transfers are no longer
performed simultaneously. The stationary particle is first transferred at the
left rung, then the mobile particle crosses the remote rung. This removes the
near-degenerate four-state manifold created by the simultaneous swap.
"""
from __future__ import annotations
import argparse, json, sys, time
from dataclasses import dataclass, asdict
from pathlib import Path
import numpy as np
from scipy.linalg import polar, sqrtm
from scipy.sparse import csr_matrix, diags

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from run_phase4_1_logical_gate import build_occ
from run_phase4_1_logical_gate_strang import exact_logical_frame
from antler.basis import build_basis
from antler.phase1 import hop_table

@dataclass
class DConfig:
    L:int=14; N:int=2; J1:float=.4; J2:float=1.; JPERP:float=.1
    DEPTH:float=-4.; R_LOOP:int=4; T_TOTAL:float=20000.
    OUT_END:float=.35; LEFT_SWAP_END:float=.45; REMOTE_SWAP_END:float=.55; RETURN_END:float=.90

def sin2(s): return float(np.sin(.5*np.pi*np.clip(s,0,1))**2)

def add_discrete_trap(mu,leg,x,depth,L):
    x=float(np.clip(x,0,L-1)); i=int(np.floor(x)); s=x-i
    if i>=L-1:
        mu[2*(L-1)+leg]+=depth; return
    q=sin2(s)
    mu[2*i+leg]+=depth*(1-q)
    mu[2*(i+1)+leg]+=depth*q

def rung_crossfade(mu,rung,leg_from,leg_to,s,depth):
    q=sin2(s)
    mu[2*rung+leg_from]+=depth*(1-q)
    mu[2*rung+leg_to]+=depth*q

def mu_seq(u,exchange,cfg):
    mu=np.zeros(2*cfg.L); D=cfg.DEPTH; R=cfg.R_LOOP
    # spectator branch
    mu[-2]=D; mu[-1]=D
    a,b,c,d=cfg.OUT_END,cfg.LEFT_SWAP_END,cfg.REMOTE_SWAP_END,cfg.RETURN_END
    if u<a:
        x=R*(u/a); add_discrete_trap(mu,0,x,D,cfg.L); mu[1]+=D
    elif u<b:
        # mobile particle held at remote leg 0; move stationary particle 1 -> 0
        mu[2*R]+=D
        if exchange: rung_crossfade(mu,0,1,0,(u-a)/(b-a),D)
        else: mu[1]+=D
    elif u<c:
        if exchange:
            # stationary particle now fixed at (0,0), mobile particle crosses R: 0 -> 1
            mu[0]+=D; rung_crossfade(mu,R,0,1,(u-b)/(c-b),D)
        else:
            mu[2*R]+=D; mu[1]+=D
    elif u<d:
        s=(u-c)/(d-c); x=R*(1-s)
        if exchange:
            add_discrete_trap(mu,1,x,D,cfg.L); mu[0]+=D
        else:
            add_discrete_trap(mu,0,x,D,cfg.L); mu[1]+=D
    else:
        # closed final Hamiltonian
        mu[0]+=D; mu[1]+=D
    return mu

def no_global(U): return U*np.exp(-.5j*np.angle(np.linalg.det(U)))
def favg(U,V): return float((abs(np.trace(V.conj().T@U))**2+2)/6)
def arr(A): return {'real':A.real.tolist(),'imag':A.imag.tolist()}

def simulate(theta,exchange,cfg,dt,states,index,table,OCC,diagnostic=False):
    d=len(states);M=2*cfg.L;rows,cols,mJ,nmid=table
    amp=mJ*np.exp(1j*theta*nmid);one=csr_matrix((amp,(rows,cols)),shape=(d,d));Hhop=one+one.conj().T
    mu0=mu_seq(0,exchange,cfg);muf=mu_seq(1,exchange,cfg)
    assert np.linalg.norm(mu0-muf)<1e-12
    H0=Hhop+diags(OCC@mu0);U0,info=exact_logical_frame(H0,index,M)
    E,V=np.linalg.eigh(Hhop.toarray());nseg=int(round(cfg.T_TOTAL/dt));dt=cfg.T_TOTAL/nseg
    Uhop=(V*np.exp(-1j*E*dt))@V.conj().T;Psi=U0.copy(); checkpoints=[]
    bounds=[cfg.OUT_END,cfg.LEFT_SWAP_END,cfg.REMOTE_SWAP_END,cfg.RETURN_END,1.0]
    next_b=0
    for k in range(nseg):
        u=(k+.5)/nseg;v=OCC@mu_seq(u,exchange,cfg);ph=np.exp(-.5j*dt*v)[:,None];Psi=ph*(Uhop@(ph*Psi))
        while diagnostic and next_b<len(bounds) and (k+1)/nseg>=bounds[next_b]:
            S=U0.conj().T@Psi
            checkpoints.append({'u':bounds[next_b], 'fixed_code_leak':(1-np.sum(abs(S)**2,axis=0)).tolist()})
            next_b+=1
    S=U0.conj().T@Psi
    return {'S':S,'U':polar(S)[0],'leak':1-np.sum(abs(S)**2,axis=0),'frame':info,'nseg':nseg,'dt':dt,'checkpoints':checkpoints}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--theta',type=float,default=.3);ap.add_argument('--dt',type=float,default=.25)
    ap.add_argument('--T',type=float,default=20000.);ap.add_argument('--R',type=int,default=4);ap.add_argument('--depth',type=float,default=-4.)
    ap.add_argument('--left-swap',type=float,default=.10,help='fraction of total time for left rung transfer')
    ap.add_argument('--remote-swap',type=float,default=.10,help='fraction of total time for remote rung transfer')
    ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    out_end=.35; left_end=out_end+args.left_swap; remote_end=left_end+args.remote_swap
    if remote_end>=.78: raise ValueError('swap windows leave too little return time')
    cfg=DConfig(DEPTH=args.depth,R_LOOP=args.R,T_TOTAL=args.T,OUT_END=out_end,LEFT_SWAP_END=left_end,REMOTE_SWAP_END=remote_end,RETURN_END=.90)
    M=2*cfg.L;states,index=build_basis(M,cfg.N);table=hop_table(cfg.L,cfg.J1,cfg.J2,cfg.JPERP,states,index);OCC=build_occ(states,M)
    t0=time.time();runs={}
    for th in (args.theta,-args.theta):
        for name,ex in [('rt',False),('ex',True)]:
            print('run',name,th,'depth',args.depth,'R',args.R,'windows',args.left_swap,args.remote_swap,flush=True)
            r=simulate(th,ex,cfg,args.dt,states,index,table,OCC,diagnostic=True);runs[(name,th)]=r;print('leak',r['leak'],flush=True)
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
    for (name,th),r in runs.items():payload['runs'][f'{name}_{th:+.3f}']={'S':arr(r['S']),'U':arr(r['U']),'leak':r['leak'].tolist(),'frame':r['frame'],'nseg':r['nseg'],'dt':r['dt'],'checkpoints':r['checkpoints']}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(payload,indent=2));print(json.dumps(metrics,indent=2));print('saved',args.out,'runtime',payload['runtime_s'])
if __name__=='__main__':main()
