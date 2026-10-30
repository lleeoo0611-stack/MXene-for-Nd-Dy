#!/usr/bin/env python3
"""Generate CP2K 2026.2 Nd3+/Dy3+ hydration cluster *starting geometries*.

NOT a computed hydration preference, bulk-water coordination, or membrane result.
Based on Lu et al., JCTC 19 (2023) 82-96, DOI 10.1021/acs.jctc.2c00922.
"""
from __future__ import annotations
import argparse
import itertools
import math
from pathlib import Path

LAN_THANIDES = {"Nd": 2.50, "Dy": 2.38}  # START guesses, Angstrom; NOT fitted equilibrium bonds
BOX = 26.0
ORIGIN = (13.0, 13.0, 13.0)

def vadd(a, b):
    return tuple(x+y for x, y in zip(a, b))
def vscale(s, a):
    return tuple(s*x for x in a)
def dot(a,b):
    return sum(x*y for x,y in zip(a,b))
def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def unit(a):
    n=math.sqrt(dot(a,a))
    if n<1e-9:raise ValueError("Zero vector")
    return vscale(1./n,a)

def directions(n):
    if n == 8:
        return [unit(p) for p in itertools.product((-1.,1.), repeat=3)]
    if n == 9:
        # Deterministic Fibonacci-shell positions; NOT claimed to be optimized
        return [((1-z*z)**.5*math.cos(i*math.pi*(3-math.sqrt(5))),
                 (1-z*z)**.5*math.sin(i*math.pi*(3-math.sqrt(5))),z)
                for i in range(n) for z in [1-2*(i+.5)/n]]
    raise ValueError("Only 8/9 water initial shells implemented")

def coords(ion, nw):
    if ion not in LAN_THANIDES:raise ValueError(ion)
    xyz=[(ion,ORIGIN)]
    r=LAN_THANIDES[ion]
    theta=math.radians(52.26)
    for u in directions(nw):
        o=vadd(ORIGIN,vscale(r,u))
        t=unit(cross(u,(0,0,1) if abs(u[2]) < .9 else (0,1,0)))
        h1=vadd(o,vscale(.9572,vadd(vscale(math.cos(theta),u),vscale(math.sin(theta),t))))
        h2=vadd(o,vscale(.9572,vadd(vscale(math.cos(theta),u),vscale(-math.sin(theta),t))))
        xyz.extend((("O",o),("H",h1),("H",h2)))
    return xyz

def make_input(ion, nw, run_type="ENERGY", cutoff=450):
    if run_type not in ("ENERGY","GEO_OPT"): raise ValueError(run_type)
    lines=[f"{el:<2s} {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}" for el,p in coords(ion,nw)]
    # 4f-in-core LnPP2 PBE q11 for both Nd and Dy; charges account for ion +3.
    # POTENTIALS_COMBINED = LnPP2_POTENTIALS + GTH_POTENTIALS (no redistribution).
    return f"""&GLOBAL
  PROJECT {ion.lower()}_w{nw}
  RUN_TYPE {run_type}
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
      REL_CUTOFF 55
      NGRIDS 4
    &END MGRID
    &POISSON
      PERIODIC NONE
      PSOLVER MT
    &END POISSON
    &SCF
      EPS_SCF 1.0E-6
      MAX_SCF 100
      SCF_GUESS ATOMIC
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF
        MAX_SCF 8
        EPS_SCF 1.0E-6
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC {BOX:.1f} {BOX:.1f} {BOX:.1f}
      PERIODIC NONE
    &END CELL
    &COORD
{chr(10).join('      '+line for line in lines)}
    &END COORD
    &KIND {ion}
      BASIS_SET DZVP-MOLOPT-GTH-q11
      POTENTIAL GTH-PBE-q11
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
{'''&MOTION
  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER 80
    MAX_DR 0.04
  &END GEO_OPT
&END MOTION
''' if run_type == 'GEO_OPT' else ''}
"""

def write_all(out,run_type,cutoff):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    for ion in LAN_THANIDES:
        for nw in (8,9):
            f=out/f"{ion.lower()}_w{nw}.inp"
            f.write_text(make_input(ion,nw,run_type,cutoff))
            print(f)
    (out/"GEOMETRIES_ONLY.txt").write_text(
        "The generated structures are geometrical guesses. The CP2K calculations must"
        " converge before any electronic-structure results are reported. Never compare"
        " total energies of Nd versus Dy directly; different Hamiltonians/atomic references."
        " No isolated hydration cluster can establish bulk CN, solvation free energy or Nd/Dy permeance.\n")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",default="simulations/cp2k_hydration/generated")
    p.add_argument("--run-type",choices=["ENERGY","GEO_OPT"],default="ENERGY")
    p.add_argument("--cutoff",type=int,default=450)
    a=p.parse_args()
    write_all(a.out,a.run_type,a.cutoff)

if __name__=="__main__":main()
