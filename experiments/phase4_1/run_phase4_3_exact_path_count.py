"""Exact JW-string phase count for the sequential digital shuttle.

No Schrödinger evolution is used here. We enumerate the ideal localized Fock
configurations and multiply the phase exponents dictated by the rung-major
fractional Jordan-Wigner hopping convention.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path


def between(p,q):
    lo,hi=sorted((p,q)); return range(lo+1,hi)

def hop_exponent(occupied,src,dst):
    """Exponent nu in matrix element ~ exp(i theta nu).

    In the implemented convention, a decreasing-index hop l->k has +n_mid;
    its Hermitian-conjugate increasing-index hop k->l has -n_mid.
    """
    if src not in occupied or dst in occupied: raise ValueError((occupied,src,dst))
    nmid=sum(1 for m in between(src,dst) if m in occupied)
    return nmid if src>dst else -nmid

def apply(occupied,src,dst):
    nu=hop_exponent(occupied,src,dst)
    out=set(occupied);out.remove(src);out.add(dst)
    return out,nu

def exchange_path(R):
    occ={0,1}; seq=[]
    # mobile particle: leg 0, left -> remote
    for r in range(R): seq.append((2*r,2*(r+1),'out leg0'))
    # stationary particle swaps rung at left; mobile swaps rung at R
    seq.append((1,0,'left rung swap'))
    seq.append((2*R,2*R+1,'remote rung swap'))
    # mobile particle returns on leg 1
    for r in range(R,0,-1): seq.append((2*r+1,2*(r-1)+1,'return leg1'))
    return occ,seq

def roundtrip_path(R):
    occ={0,1};seq=[]
    for r in range(R):seq.append((2*r,2*(r+1),'out leg0'))
    for r in range(R,0,-1):seq.append((2*r,2*(r-1),'return leg0'))
    return occ,seq

def evaluate(kind,R):
    occ,seq=exchange_path(R) if kind=='exchange' else roundtrip_path(R)
    history=[];total=0
    for src,dst,label in seq:
        before=sorted(occ);occ,nu=apply(occ,src,dst);total+=nu
        history.append({'label':label,'src':src,'dst':dst,'before':before,'after':sorted(occ),'nu':nu})
    return {'kind':kind,'R':R,'total_exponent':total,'phase':f'exp(i theta * {total})','final':sorted(occ),'history':history}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--R',type=int,default=4);ap.add_argument('--out',type=Path,default=Path('results/phase4_3/exact_path_count.json'));args=ap.parse_args()
    ex=evaluate('exchange',args.R);rt=evaluate('roundtrip',args.R)
    payload={'exchange':ex,'roundtrip':rt,'differential_exponent':ex['total_exponent']-rt['total_exponent']}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(payload,indent=2))
    print('exchange exponent',ex['total_exponent']);print('round-trip exponent',rt['total_exponent']);print('differential',payload['differential_exponent'])
if __name__=='__main__':main()
