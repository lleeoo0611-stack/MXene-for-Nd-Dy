#!/usr/bin/env python3
"""First CP2K Ln3+ monohydrate single-point input; geometry is NOT optimized."""
from __future__ import annotations
import argparse
from math import radians,cos,sin
from pathlib import Path
def make_input(ion,box=12,cutoff=250,eps_scf=5e-5):
    if ion not in ("Nd","Dy"):raise ValueError(ion)
    if box<10:raise ValueError("Too small nonperiodic cell")
    d={"Nd":2.50,"Dy":2.38}[ion] # initial separation guesses
    center=box/2
    dx=.9572*cos(radians(52.26))
    dy=.9572*sin(radians(52.26))
    xyz=[
       (ion,center,center,center),
       ("O",center+d,center,center),
       ("H",center+d+dx,center+dy,center),
       ("H",center+d+dx,center-dy,center),
    ]
    positions="\n".join(
        "     %-2s %.8f %.8f %.8f"%p for p in xyz
    )
    return f"""&GLOBAL
  PROJECT {ion.lower()}_w1_dft
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_LnPP2
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME POTENTIALS_COMBINED
    CHARGE 3
    MULTIPLICITY 1
    &QS
      METHOD GPW
    &END QS
    &MGRID
      CUTOFF {cutoff}
      REL_CUTOFF 40
      NGRIDS 4
    &END MGRID
    &POISSON
      PERIODIC NONE
      PSOLVER MT
    &END POISSON
    &SCF
      EPS_SCF {eps_scf:.7g}
      MAX_SCF 160
      SCF_GUESS ATOMIC
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF
        MAX_SCF 8
        EPS_SCF {eps_scf:.7g}
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC {box} {box} {box}
      PERIODIC NONE
    &END CELL
    &COORD
{positions}
    &END COORD
    &KIND {ion}
      BASIS_SET DZVP-MOLOPT-GTH-q11
      POTENTIAL GTH-PBE-q11
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
      BASIS_SET DZVP-MOLOPT-GTH-q1
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
"""
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--ion",choices=["Nd","Dy"],required=True)
    p.add_argument("--out",default="build/dft")
    p.add_argument("--cutoff",type=int,default=250)
    p.add_argument("--cell",type=float,default=12)
    a=p.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    f=out/f"{a.ion.lower()}_w1.inp"
    f.write_text(make_input(a.ion,box=a.cell,cutoff=a.cutoff),encoding="utf8")
    print(f)
if __name__=="__main__":
    main()
