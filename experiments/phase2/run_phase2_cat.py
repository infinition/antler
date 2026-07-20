"""ANTLER 0.1, run_phase2_cat.py
Phase 2: correction de l'encodage chat et validation du plan de controle.

Contenu:
  TEST A: decomposition par jambe du zero intra-code (code bords opposes)
  TEST B: encodage chat corrige {|L1 L2>, |R1 R2>}:
          action premier ordre P dH/dth P, action differentielle sous
          delta_asym, courbure F(theta, delta_asym) avec scaling h^2.

Note: l'encodage |B_L A_L> du run precedent etait invalide: le produit
bosonique symetrise de deux combinaisons orthogonales de jambes
interfere destructivement sur la configuration dominante (capture 0.34).
La base chat correcte est {|L1 L2>, |R1 R2>}.

Usage: python run_phase2_cat.py
"""
import numpy as np
from scipy.linalg import logm
from antler.basis import build_basis
from antler.edge import edge_orbitals
from antler.logical import projector, two_particle_state
from antler.model import single_particle_matrix
from antler.phase1 import build_h, hop_table, instantaneous_doublet, mu_diagonal

L=14; N=2; J1,J2,JP=0.4,1.0,0.1; DREF=-0.35
M=2*L
states,index=build_basis(M,N); d=len(states)
table=hop_table(L,J1,J2,JP,states,index)
rows,cols,mJ,nmid=table
mu_ref=np.zeros(M); mu_ref[[0,1,M-2,M-1]]=DREF
h1=single_particle_matrix(L,J1,J2,JP,mu=mu_ref)
orbs,_,_=edge_orbitals(h1,L)

# ------- TEST A: le zero exact vient-il d'une compensation inter-jambes?
v0=two_particle_state((orbs["L1"]+orbs["L2"])/np.sqrt(2),(orbs["R1"]-orbs["R2"])/np.sqrt(2),states,index)
v1=two_particle_state((orbs["L1"]-orbs["L2"])/np.sqrt(2),(orbs["R1"]+orbs["R2"])/np.sqrt(2),states,index)
U_sep,_=np.linalg.qr(np.column_stack([v0,v1])); P_sep=projector(U_sep)

def dHdth_leg(theta,sigma):
    # sauts de jambe du mover sur la jambe sigma uniquement (k%2==sigma)
    Hd=np.zeros((d,d),dtype=complex)
    # nmid>0 <=> saut de jambe; la jambe du mover = parite de l'indice k du site bas
    # on doit reconstruire k depuis la table: refaisons une table annotée
    return Hd

# reconstruire une table annotee par jambe
from antler.model import hop_list
from antler.basis import between_mask, popcount
hops=hop_list(L,J1,J2,JP)
ann=[]
for col,s in enumerate(states):
    s=int(s)
    for k,l,J in hops:
        if ((s>>l)&1) and not((s>>k)&1):
            new=s^(1<<l)^(1<<k)
            nm=popcount(s&between_mask(k,l))
            ann.append((index[new],col,-J,nm,k%2 if l-k==2 else -1))
ann=np.array(ann,dtype=float)
def build_dH(theta,leg=None):
    Hd=np.zeros((d,d),dtype=complex)
    sel=ann[:,3]>0 if leg is None else (ann[:,4]==leg)
    a=ann[sel]
    np.add.at(Hd,(a[:,0].astype(int),a[:,1].astype(int)),
              a[:,2]*(1j*a[:,3])*np.exp(1j*0.742*a[:,3]))
    return Hd+Hd.conj().T
G1=U_sep.conj().T@build_dH(0.742,leg=0)@U_sep
G2=U_sep.conj().T@build_dH(0.742,leg=1)@U_sep
print("TEST A (code bords opposes):")
print(f"  ||P dH_leg1 P|| = {np.linalg.norm(G1):.3e}")
print(f"  ||P dH_leg2 P|| = {np.linalg.norm(G2):.3e}")
print(f"  ||somme||       = {np.linalg.norm(G1+G2):.3e}  (compensation si sommes nulles et termes non nuls)")

# ------- TEST B: vraie base chat {|L1 L2>, |R1 R2>}
vLL=two_particle_state(orbs["L1"],orbs["L2"],states,index)
vRR=two_particle_state(orbs["R1"],orbs["R2"],states,index)
U_cat,_=np.linalg.qr(np.column_stack([vLL,vRR])); P_cat=projector(U_cat)
mu_edge=np.zeros(M); mu_edge[[0,1,M-2,M-1]]=1.0
muL=np.zeros(M); muL[[0,1]]=1.0
muR=np.zeros(M); muR[[M-2,M-1]]=1.0
dL=mu_diagonal(states,muL); dR=mu_diagonal(states,muR)
def H_at(th,da):
    return build_h(th,table,d,(DREF+da/2)*dL+(DREF-da/2)*dR)
fr,cap,sp,Ep=instantaneous_doublet(H_at(0,0),P_cat)
print(f"\nTEST B (chat L1L2/R1R2): E={Ep.round(6)}, split={sp:.2e}, capture={cap:.5f}/2")
Gc=U_cat.conj().T@build_dH(0.742)@U_cat
print(f"  ||P dH/dth P|| (1er ordre) = {np.linalg.norm(Gc):.3e}")
print("  E_doublet(theta) suivi:")
for da in (0.0,0.10):
    Es=[]
    for th in (0.0,0.3,0.6,0.9):
        _,c2,s2_,E2=instantaneous_doublet(H_at(th,da),P_cat)
        Es.append((E2,c2))
    E0,c0=Es[0]; E3,c3=Es[-1]
    dE=E3-E0
    print(f"   da={da:.2f}: E(0)={E0.round(5)} (cap {c0:.3f}), E(0.9)={E3.round(5)} (cap {c3:.3f})")
    print(f"            shift branches={dE.round(6)}, differentiel={dE[1]-dE[0]:+.3e}")
def frame_at(th,da): return instantaneous_doublet(H_at(th,da),P_cat)[0]
def logW(th,da,h):
    fs=[frame_at(th,da),frame_at(th+h,da),frame_at(th+h,da+h),frame_at(th,da+h)]
    W=np.eye(2,dtype=complex)
    for a in range(4):
        A,B=fs[a],fs[(a+1)%4]
        u,_,vh=np.linalg.svd(A.conj().T@B); W=W@(u@vh)
    return float(np.linalg.norm(logm(W)))
print("  courbure F(theta,da) en (0.6,0.05), scaling:")
prev=None
for h in (0.08,0.04,0.02):
    lw=logW(0.6,0.05,h)
    r="" if prev is None else f", ratio={prev/lw:.2f}"
    print(f"   h={h:.2f}: |logW|={lw:.3e}, ||F||={lw/h**2:.4f}{r}")
    prev=lw
