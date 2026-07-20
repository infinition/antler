"""ANTLER Phase 4.4b: adiabatic Hellmann-Feynman string response.

Tracks the left/right logical instantaneous eigenstates at theta=0 along the
exchange and round-trip schedules. Integrates <dH/dtheta> to predict the
odd-in-theta dynamical phase slope. The residual against the finite-theta gate
slope is the geometric/statistical contribution not explained dynamically.
"""
from __future__ import annotations
import argparse, json, sys, time
from dataclasses import asdict
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import eigsh

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from run_phase4_1_logical_gate import Config, protocol_mu, build_occ
from run_phase4_1_logical_gate_strang import exact_logical_frame
from antler.basis import build_basis
from antler.phase1 import hop_table


def track_pair(H, prev, k=14):
    # Low-energy window; jointly assign two distinct eigenvectors by continuity.
    E,V=eigsh(H,k=k,which='SA',tol=1e-10)
    order=np.argsort(E); E=E[order]; V=V[:,order]
    O=np.abs(prev.conj().T@V)**2
    best=None
    for a in range(k):
        for b in range(k):
            if a==b: continue
            val=O[0,a]+O[1,b]
            if best is None or val>best[0]: best=(val,a,b)
    inds=[best[1],best[2]]
    out=V[:,inds].copy()
    ovs=np.array([O[0,inds[0]],O[1,inds[1]]],float)
    for j in range(2):
        z=np.vdot(prev[:,j],out[:,j])
        if abs(z)>0: out[:,j]*=np.conj(z)/abs(z)
    return out, E[inds].astype(float), ovs


def response(exchange,cfg,nsamp,states,index,table,OCC):
    d=len(states);M=2*cfg.L
    rows,cols,mJ,nmid=table
    one=csr_matrix((mJ,(rows,cols)),shape=(d,d))
    Hhop=(one+one.conj().T).tocsr()
    done=csr_matrix((1j*mJ*nmid,(rows,cols)),shape=(d,d))
    Htheta=(done+done.conj().T).tocsr()
    H0=Hhop+diags(OCC@protocol_mu(0.,exchange,cfg))
    U,info=exact_logical_frame(H0,index,M)
    # Exact energies of the localized degenerate frame.
    energies=np.real(np.diag(U.conj().T@(H0@U))).copy()
    us=np.linspace(0.,1.,nsamp)
    G=np.zeros((nsamp,2)); Es=np.zeros((nsamp,2)); ovs=np.ones((nsamp,2))
    G[0]=np.real(np.diag(U.conj().T@(Htheta@U))); Es[0]=energies
    frames=[U[:,0].copy(),U[:,1].copy()]
    for q,u in enumerate(us[1:],1):
        H=Hhop+diags(OCC@protocol_mu(float(u),exchange,cfg))
        prev=np.column_stack(frames)
        fr,energies,ovs[q]=track_pair(H,prev)
        frames=[fr[:,0].copy(),fr[:,1].copy()]
        Es[q]=energies
        G[q]=np.real(np.diag(fr.conj().T@(Htheta@fr)))
    # Integral in physical time t=T*u.
    I=cfg.T_TOTAL*np.trapezoid(G,us,axis=0)
    Edyn=cfg.T_TOTAL*np.trapezoid(Es,us,axis=0)
    return {
        'integral_dH':I.tolist(),
        'integral_E':Edyn.tolist(),
        'min_step_overlap':np.min(ovs[1:],axis=0).tolist(),
        'mean_step_overlap':np.mean(ovs[1:],axis=0).tolist(),
        'max_abs_response':np.max(np.abs(G),axis=0).tolist(),
        'frame':info,
        'samples':nsamp,
        'trace':{'u':us[::max(1,nsamp//100)].tolist(),
                 'g':G[::max(1,nsamp//100)].tolist(),
                 'E':Es[::max(1,nsamp//100)].tolist()}
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--samples',type=int,default=401)
    ap.add_argument('--T',type=float,default=20000.)
    ap.add_argument('--R',type=float,default=4.)
    ap.add_argument('--A',type=float,default=2.6)
    ap.add_argument('--W',type=float,default=1.)
    ap.add_argument('--measured-slope',type=float,default=-0.996426)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    cfg=Config(T_TOTAL=args.T,R_LOOP=args.R,A_WELL=args.A,W_WELL=args.W)
    M=2*cfg.L;states,index=build_basis(M,cfg.N)
    table=hop_table(cfg.L,cfg.J1,cfg.J2,cfg.JPERP,states,index)
    OCC=build_occ(states,M)
    t0=time.time(); runs={}
    for name,ex in [('rt',False),('ex',True)]:
        print('track',name,flush=True)
        runs[name]=response(ex,cfg,args.samples,states,index,table,OCC)
        print(name,'I',runs[name]['integral_dH'],'min overlap',runs[name]['min_step_overlap'],flush=True)
    Iex=np.array(runs['ex']['integral_dH']); Irt=np.array(runs['rt']['integral_dH'])
    dyn=-float((Iex[0]-Iex[1])-(Irt[0]-Irt[1]))
    residual=float(args.measured_slope-dyn)
    payload={'config':asdict(cfg),'samples':args.samples,'measured_odd_slope':args.measured_slope,
             'predicted_dynamic_odd_slope':dyn,'residual_geometric_slope':residual,
             'runs':runs,'runtime_s':time.time()-t0}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(payload,indent=2))
    print(json.dumps({k:payload[k] for k in ('measured_odd_slope','predicted_dynamic_odd_slope','residual_geometric_slope','runtime_s')},indent=2))

if __name__=='__main__':main()
