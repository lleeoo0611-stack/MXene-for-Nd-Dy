#!/usr/bin/env python3
"""V3: PARTIAL Nd/Dy–citrate aqueous speciation (not adsorption, MD or DFT).

Constants: ACS Earth Space Chem (2024) DOI 10.1021/acsearthspacechem.4c00009
Ln3+ + Cit3- + HCit2- = Ln(Cit)(HCit)2-:
log10(beta0): Nd=13.05, Dy=13.28
Ln3+ + OH- = LnOH2+: log10(beta0): Nd=5.82, Dy=6.41
Bulk citrate pKa=3.13,4.76,6.40: ACS Omega (2019),
DOI 10.1021/acsomega.9b00270.

EXCLUDES all other citrate complexes, nitrate, carbonate, precipitation,
surface-immobilized citrate and MXene. Fixed ionic strength Davies
approximation at 25 degC; results MUST NOT be treated as full speciation.
"""
from __future__ import annotations
import argparse, csv, math
from pathlib import Path
PKA=(3.13,4.76,6.40)
BETA_MIXED={"Nd":13.05,"Dy":13.28}
BETA_HYDROXO={"Nd":5.82,"Dy":6.41}
def gamma(z,I):
    if not 0<=I<=0.5:raise ValueError("I outside screening range 0..0.5M")
    t=math.sqrt(I)
    return 10**(-0.51*z*z*(t/(1+t)-.3*I))
def acid_fractions(pH,I,pka=PKA):
    """H3Cit/H2Cit-/HCit2-/Cit3-; activity-corrected concentrations."""
    H=10**(-pH);g1,g2,g3=(gamma(z,I) for z in (1,2,3))
    k1,k2,k3=(10**(-x) for x in pka)
    w0=1.;w1=w0*k1/(H*g1);w2=w1*k2*g1/(H*g2);w3=w2*k3*g2/(H*g3)
    sm=w0+w1+w2+w3
    return tuple(x/sm for x in (w0,w1,w2,w3))
def state(pH,I,L,metal_total,shifts=None,pka=PKA):
    f=acid_fractions(pH,I,pka);hc=L*f[2];cit=L*f[3]
    g3=gamma(3,I);g2=gamma(2,I)
    aoh=10**(pH-14)
    shifts=shifts or {}
    rows={};bound=0
    for el,total in metal_total.items():
        b=10**(BETA_MIXED[el]+shifts.get(el,0))*g3*g3*hc*cit
        o=10**BETA_HYDROXO[el]*(g3/g2)*aoh
        free=total/(1+b+o);complexed=free*b;hydroxo=free*o
        bound+=2*complexed
        rows[el]={"free":free,"mixed":complexed,"hydroxo":hydroxo}
    return rows,bound
def solve(pH,I,Ltotal,Nd=1e-4,Dy=1e-4,shifts=None,pka=PKA):
    if min(Ltotal,Nd,Dy)<0:raise ValueError("negative total")
    metals={"Nd":Nd,"Dy":Dy}
    low,high=0.,Ltotal
    for _ in range(150):
        half=(low+high)/2
        _,bound=state(pH,I,half,metals,shifts,pka)
        if half+bound>Ltotal:high=half
        else:low=half
    free_lig=(low+high)/2
    states,bound=state(pH,I,free_lig,metals,shifts,pka)
    out={"pH":pH,"I_M":I,"citrate_total_M":Ltotal,
         "free_citrate_total_M":free_lig,"bound_citrate_M":bound,
         "citrate_balance_error_M":free_lig+bound-Ltotal}
    for el,total in metals.items():
        for label,val in states[el].items():
            out[el+"_"+label+"_M"]=val
            out[el+"_"+label+"_fraction"]=val/total if total else 0.
        out[el+"_metal_balance_error_M"]=sum(states[el].values())-total
    return out
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",default="results/v3_citrate_aqueous_partial.csv")
    args=p.parse_args()
    rows=[solve(ph,I,CmM*0.001) for I in (0.01,0.1)
          for ph in (4.5,5.5,6.5) for CmM in (0,0.1,0.5)]
    target=Path(args.output);target.parent.mkdir(exist_ok=True,parents=True)
    with target.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    print("Computed",len(rows),"partial-equilibrium conditions;",target)
    r=solve(5.5,0.01,0.0001)
    print("pH5.5 I0.01M citrate0.1mM:",
          "Nd",round(100*r["Nd_mixed_fraction"],2),"Dy",round(100*r["Dy_mixed_fraction"],2),
          "percent Ln(Cit)(HCit)2-. NOT full Ln–citrate speciation.")
if __name__=="__main__":main()
