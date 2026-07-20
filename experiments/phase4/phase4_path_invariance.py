from __future__ import annotations
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('MKL_NUM_THREADS','1')
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import expm_multiply, eigsh
from antler.basis import build_basis
from antler.phase1 import hop_table
import run_phase4_shuttle_validated as base

THETA = 0.3
R_TEST = 4


def one_run(args):
    T_total, nseg, proc, theta = args
    exchange = proc == 'ex'
    base.R_LOOP = R_TEST
    L, N = base.L, base.N
    M = 2 * L
    states, index = build_basis(M, N)
    d = len(states)
    table = hop_table(L, base.J1, base.J2, base.JPERP, states, index)
    OCC = np.zeros((d, M))
    for p, state in enumerate(states):
        for k in range(M):
            if (int(state) >> k) & 1:
                OCC[p, k] = 1.0
    rows, cols, mJ, nmid = table
    amp = mJ * np.exp(1j * theta * nmid)
    Hs = csr_matrix((amp, (rows, cols)), shape=(d, d))
    Hhop = Hs + Hs.conj().T
    H0 = Hhop + diags(OCC @ base.mu_sites(0.0, exchange, M))
    E0, V0 = eigsh(H0, k=2, which='SA')
    psi = V0[:, np.argmin(E0)].astype(complex)
    j = int(np.argmax(np.abs(psi)))
    psi *= np.conj(psi[j]) / abs(psi[j])
    dt = T_total / nseg
    for a in range(nseg):
        u = (a + 0.5) / nseg
        H = Hhop + diags(OCC @ base.mu_sites(u, exchange, M))
        psi = expm_multiply(-1j * dt * H.tocsc(), psi)
    psi /= np.linalg.norm(psi)
    Hf = Hhop + diags(OCC @ base.mu_sites(1.0, exchange, M))
    Ef, Vf = eigsh(Hf, k=2, which='SA')
    fid_final = float(abs(np.vdot(Vf[:, np.argmin(Ef)], psi)))
    return {'T':T_total,'nseg':nseg,'proc':proc,'theta':theta,
            'psi_real':psi.real.tolist(),'psi_imag':psi.imag.tolist(),
            'fid_final':fid_final}


def analyze(records):
    out=[]
    for T,nseg in sorted({(r['T'],r['nseg']) for r in records}):
        rr=[r for r in records if r['T']==T and r['nseg']==nseg]
        def psi(proc,th):
            r=next(x for x in rr if x['proc']==proc and abs(x['theta']-th)<1e-12)
            return np.array(r['psi_real'])+1j*np.array(r['psi_imag'])
        vals={}
        for proc in ('rt','ex'):
            p0=psi(proc,0.0)
            for th in (THETA,-THETA):
                vals[(proc,th)]=np.vdot(p0,psi(proc,th))
        odd={}
        for proc in ('rt','ex'):
            ratio=vals[(proc,THETA)]/vals[(proc,-THETA)]
            odd[proc]=0.5*float(np.angle(ratio))
        dphi=float(np.angle(np.exp(1j*(odd['ex']-odd['rt']))))
        out.append({'T':T,'nseg':nseg,'odd_ex':odd['ex'],'odd_rt':odd['rt'],
                    'dphi':dphi,'slope':dphi/THETA,
                    'min_fid_final':min(r['fid_final'] for r in rr)})
    return out

if __name__=='__main__':
    settings=[(20000.0,5000)]
    tasks=[(T,nseg,proc,th) for T,nseg in settings for proc in ('rt','ex') for th in (0.0,THETA,-THETA)]
    records=[]
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs=[ex.submit(one_run,t) for t in tasks]
        for f in as_completed(futs):
            r=f.result(); records.append(r)
            print(f"done T={r['T']:.0f} n={r['nseg']} {r['proc']} th={r['theta']:+.1f} fid={r['fid_final']:.4f}",flush=True)
    summary=analyze(records)
    with open('results/phase4_path_R4.json','w') as f:
        json.dump({'summary':summary,'records':records},f,indent=2)
    print(json.dumps(summary,indent=2))
