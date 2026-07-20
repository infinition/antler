"""ANTLER Phase 4.5: isolated two-level handoff lemma.

A localized particle is adiabatically transferred between two sites through
H(s) = [[-D(1-q), -J exp(+i phi)], [-J exp(-i phi), -D q]].
Parallel transport of the instantaneous ground state must inherit exactly the
oriented link phase. This is the local building block of the digital braid.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np


def ramp(s,kind):
    if kind=='linear': return s
    if kind=='sin2': return np.sin(.5*np.pi*s)**2
    if kind=='smoothstep': return s*s*(3-2*s)
    raise ValueError(kind)


def transport(phi,D,J,kind,n=4001):
    prev=None
    for s in np.linspace(0,1,n):
        q=ramp(float(s),kind)
        H=np.array([[-D*(1-q),-J*np.exp(1j*phi)],
                    [-J*np.exp(-1j*phi),-D*q]],complex)
        E,V=np.linalg.eigh(H); v=V[:,0]
        if prev is None:
            # initial gauge: dominant source amplitude real positive
            v*=np.conj(v[np.argmax(abs(v))])/abs(v[np.argmax(abs(v))])
        else:
            z=np.vdot(prev,v); v*=np.conj(z)/abs(z)
        prev=v
    amp=prev[1]
    return {'phase':float(np.angle(amp)),'target_weight':float(abs(amp)**2)}


def main():
    phis=[.2,.6,1.1]; depths=[2.,4.,8.,16.]; kinds=['linear','sin2','smoothstep']
    rows=[]
    for phi in phis:
        for D in depths:
            for kind in kinds:
                r=transport(phi,D,1.,kind)
                # In this convention source 0 -> target 1 carries -phi.
                err=float(np.angle(np.exp(1j*(r['phase']+phi))))
                rows.append({'phi':phi,'D_over_J':D,'ramp':kind,**r,'phase_error_to_minus_phi':err})
    out={'max_abs_phase_error':max(abs(r['phase_error_to_minus_phi']) for r in rows),
         'min_target_weight':min(r['target_weight'] for r in rows),'rows':rows}
    p=Path('results/phase4_5/two_level_handoff_lemma.json');p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(out,indent=2));print(json.dumps({k:out[k] for k in ('max_abs_phase_error','min_target_weight')},indent=2))

if __name__=='__main__':main()
