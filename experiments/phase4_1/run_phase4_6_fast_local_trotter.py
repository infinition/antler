"""ANTLER Phase 4.6: fast local-bond symmetric Trotter propagator.

Purpose: reproduce the exact-hopping Strang logical gate at much lower cost,
then enable deep-localization, composition, and disorder ensembles.

The kinetic propagator is decomposed into exact two-level rotations for each
physical hopping bond. A forward/reverse half-step sweep gives a symmetric
second-order formula. The on-site potential remains exact and diagonal.
"""
from __future__ import annotations
import argparse, json, time, sys
from dataclasses import dataclass, asdict
from pathlib import Path
import numpy as np
from scipy.linalg import polar, sqrtm
from scipy.sparse import csr_matrix, diags
from numba import njit

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE)); sys.path.insert(0,str(ROOT))
from run_phase4_3b_digital_sequential import DConfig, mu_digital
from run_phase4_1_logical_gate import build_occ
from run_phase4_1_logical_gate_strang import exact_logical_frame
from antler.basis import build_basis, between_mask, popcount
from antler.model import hop_list
from antler.phase1 import hop_table

@njit(cache=True)
def _apply_local_step(Psi, vdiag, dt, offsets, rows, cols, phase, Jbond):
    # V/2
    for i in range(Psi.shape[0]):
        z=np.exp(-0.5j*dt*vdiag[i])
        Psi[i,0]*=z; Psi[i,1]*=z
    # symmetric forward/reverse half-bond sweeps
    for direction in range(2):
        if direction==0:
            b0,b1,bs=0,len(Jbond),1
        else:
            b0,b1,bs=len(Jbond)-1,-1,-1
        for b in range(b0,b1,bs):
            ang=Jbond[b]*dt*0.5
            c=np.cos(ang); s=np.sin(ang)
            for q in range(offsets[b],offsets[b+1]):
                r=rows[q]; col=cols[q]; p=phase[q]
                x0=Psi[r,0]; y0=Psi[col,0]
                x1=Psi[r,1]; y1=Psi[col,1]
                # H[r,col]/J = -exp(i theta nmid) = phase
                Psi[r,0]=c*x0-1j*s*p*y0
                Psi[col,0]=c*y0-1j*s*np.conj(p)*x0
                Psi[r,1]=c*x1-1j*s*p*y1
                Psi[col,1]=c*y1-1j*s*np.conj(p)*x1
    # V/2
    for i in range(Psi.shape[0]):
        z=np.exp(-0.5j*dt*vdiag[i])
        Psi[i,0]*=z; Psi[i,1]*=z


@njit(cache=True)
def _sin2_numba(s):
    if s<0.0: s=0.0
    elif s>1.0: s=1.0
    return np.sin(0.5*np.pi*s)**2

@njit(cache=True)
def _fill_mu_digital(mu,u,exchange,L,D,R):
    for k in range(mu.size): mu[k]=0.0
    mu[2*L-2]=D; mu[2*L-1]=D
    if u<0.35:
        x=R*(u/0.35); leg=0
        i=int(np.floor(x)); ss=x-i
        if i>=L-1: mu[2*(L-1)+leg]+=D
        else:
            q=_sin2_numba(ss); mu[2*i+leg]+=D*(1-q); mu[2*(i+1)+leg]+=D*q
        mu[1]+=D
    elif u<0.45:
        ss=(u-0.35)/0.10; mu[2*R]+=D
        if exchange:
            q=_sin2_numba(ss); mu[1]+=D*(1-q); mu[0]+=D*q
        else: mu[1]+=D
    elif u<0.55:
        ss=(u-0.45)/0.10
        if exchange:
            q=_sin2_numba(ss); mu[2*R]+=D*(1-q); mu[2*R+1]+=D*q; mu[0]+=D
        else: mu[2*R]+=D; mu[1]+=D
    else:
        x=R*(1-(u-0.55)/0.45); leg=1 if exchange else 0
        i=int(np.floor(x)); ss=x-i
        if i>=L-1: mu[2*(L-1)+leg]+=D
        else:
            q=_sin2_numba(ss); mu[2*i+leg]+=D*(1-q); mu[2*(i+1)+leg]+=D*q
        if exchange: mu[0]+=D
        else: mu[1]+=D

@njit(cache=True)
def _propagate_all(Psi,nseg,cycles,dt,exchange,L,D,R,occ_a,occ_b,offsets,rows,cols,phase,Jbond):
    mu=np.zeros(2*L,np.float64); vd=np.empty(Psi.shape[0],np.float64)
    for cyc in range(cycles):
        for a in range(nseg):
            u=(a+0.5)/nseg
            _fill_mu_digital(mu,u,exchange,L,D,R)
            for i in range(Psi.shape[0]): vd[i]=mu[occ_a[i]]+mu[occ_b[i]]
            _apply_local_step(Psi,vd,dt,offsets,rows,cols,phase,Jbond)

def bond_pairs(L,J1,J2,Jperp,states,index,theta):
    offs=[0]; rr=[]; cc=[]; pp=[]; jj=[]
    for k,l,J in hop_list(L,J1,J2,Jperp):
        for col,s0 in enumerate(states):
            s0=int(s0)
            if ((s0>>l)&1) and not ((s0>>k)&1):
                new=s0^(1<<l)^(1<<k)
                nm=popcount(s0 & between_mask(k,l))
                rr.append(index[new]); cc.append(col)
                pp.append(-np.exp(1j*theta*nm))
        offs.append(len(rr)); jj.append(abs(J))
    return (np.asarray(offs,np.int64),np.asarray(rr,np.int64),
            np.asarray(cc,np.int64),np.asarray(pp,np.complex128),
            np.asarray(jj,np.float64))

def ng(U): return U*np.exp(-.5j*np.angle(np.linalg.det(U)))
def favg(U,V): return float((abs(np.trace(V.conj().T@U))**2+2)/6)
def arr(A): return {'real':A.real.tolist(),'imag':A.imag.tolist()}

def simulate(theta,exchange,cfg,dt,states,index,table,OCC,cycles=1,eps=None):
    d=len(states); M=2*cfg.L
    rows,cols,mJ,nmid=table
    amp=mJ*np.exp(1j*theta*nmid); one=csr_matrix((amp,(rows,cols)),shape=(d,d)); Hhop=one+one.conj().T
    mu0=mu_digital(0,exchange,cfg)+(0 if eps is None else eps)
    H0=Hhop+diags(OCC@mu0); U0,info=exact_logical_frame(H0,index,M)
    offsets,br,bc,bp,bj=bond_pairs(cfg.L,cfg.J1,cfg.J2,cfg.JPERP,states,index,theta)
    nseg=int(round(cfg.T_TOTAL/dt)); dt=cfg.T_TOTAL/nseg
    Psi=np.ascontiguousarray(U0.copy())
    occ_a=np.empty(d,np.int64); occ_b=np.empty(d,np.int64)
    for ii,st in enumerate(states):
        oo=[k for k in range(M) if (int(st)>>k)&1]
        occ_a[ii],occ_b[ii]=oo[0],oo[1]
    if eps is not None:
        raise NotImplementedError("fast kernel disorder path not yet enabled")
    _propagate_all(Psi,nseg,cycles,dt,exchange,cfg.L,cfg.DEPTH,cfg.R_LOOP,
                   occ_a,occ_b,offsets,br,bc,bp,bj)
    S=U0.conj().T@Psi
    return {'S':S,'U':polar(S)[0],'leak':1-np.sum(abs(S)**2,axis=0),'frame':info,'nseg':nseg,'dt':dt}

def audit(runs,theta,cycles):
    Ds={}; sv=[]; un=[]
    for th in (theta,-theta):
        for name in ('rt','ex'):
            S=runs[(name,th)]['S']; sv+=np.linalg.svd(S,compute_uv=False).tolist(); un.append(np.linalg.norm(S.conj().T@S-np.eye(2)))
        Ds[th]=ng(runs[('rt',th)]['U'].conj().T@runs[('ex',th)]['U'])
    Q=ng(Ds[theta]@Ds[-theta].conj().T)
    if np.linalg.norm(-Q-np.eye(2))<np.linalg.norm(Q-np.eye(2)): Q=-Q
    Uo=ng(polar(sqrtm(Q).astype(complex))[0])
    rel=float(np.angle(np.exp(1j*(np.angle(Uo[0,0])-np.angle(Uo[1,1])))))
    target=np.diag([np.exp(-1j*cycles*theta/2),np.exp(1j*cycles*theta/2)])
    return {'sigma_min':float(min(sv)),'sigma_max':float(max(sv)),
            'leak_worst':float(1-min(sv)**2),'unitarity_frob_max':float(max(un)),
            'odd_phase':rel,'odd_slope_per_cycle':rel/(cycles*theta),
            'odd_offdiag_norm':float(np.linalg.norm(Uo-np.diag(np.diag(Uo)))),
            'favg_target':favg(Uo,target)},Uo

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--theta',type=float,default=.3); ap.add_argument('--dt',type=float,default=.05)
    ap.add_argument('--T',type=float,default=20000.); ap.add_argument('--R',type=int,default=4); ap.add_argument('--depth',type=float,default=-4.)
    ap.add_argument('--cycles',type=int,default=1); ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args(); cfg=DConfig(DEPTH=a.depth,R_LOOP=a.R,T_TOTAL=a.T)
    M=2*cfg.L; states,index=build_basis(M,cfg.N); table=hop_table(cfg.L,cfg.J1,cfg.J2,cfg.JPERP,states,index); OCC=build_occ(states,M)
    t0=time.time(); runs={}
    for th in (a.theta,-a.theta):
        for name,ex in [('rt',False),('ex',True)]:
            print('run',name,th,'D',a.depth,'cycles',a.cycles,flush=True)
            runs[(name,th)]=simulate(th,ex,cfg,a.dt,states,index,table,OCC,a.cycles)
            print(' leak',runs[(name,th)]['leak'],flush=True)
    metrics,Uo=audit(runs,a.theta,a.cycles)
    payload={'config':asdict(cfg),'theta':a.theta,'dt':a.dt,'cycles':a.cycles,'metrics':metrics,'Uodd':arr(Uo),'runtime_s':time.time()-t0,'runs':{}}
    for (name,th),r in runs.items(): payload['runs'][f'{name}_{th:+.3f}']={'S':arr(r['S']),'U':arr(r['U']),'leak':r['leak'].tolist(),'frame':r['frame'],'nseg':r['nseg'],'dt':r['dt']}
    a.out.parent.mkdir(parents=True,exist_ok=True); a.out.write_text(json.dumps(payload,indent=2)); print(json.dumps(metrics,indent=2)); print('runtime',payload['runtime_s'])
if __name__=='__main__': main()
