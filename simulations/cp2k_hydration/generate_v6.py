"""Create four 8/9-water CP2K GEO_OPT inputs by modifying V5 initial geometries."""
from __future__ import annotations
import argparse
from pathlib import Path
from generate import make_input

def build(ion, n, box=18., cutoff=220, eps_scf=1e-5, max_iter=55):
    if ion not in ("Nd","Dy") or n not in (8,9): raise ValueError((ion,n))
    if not (15 <= box <= 30): raise ValueError("Box must be >=15 A")
    s=make_input(ion,n,run_type="GEO_OPT",cutoff=cutoff)
    s=s.replace("PRINT_LEVEL MEDIUM","PRINT_LEVEL LOW")
    s=s.replace("ABC 26.0 26.0 26.0",f"ABC {box:.1f} {box:.1f} {box:.1f}")
    s=s.replace("REL_CUTOFF 55","REL_CUTOFF 45")
    s=s.replace("EPS_SCF 1.0E-6",f"EPS_SCF {eps_scf:.2E}")
    s=s.replace("MAX_SCF 100","MAX_SCF 120")
    s=s.replace("MAX_ITER 80",f"MAX_ITER {max_iter}")
    s=s.replace("MAX_DR 0.04","MAX_DR 3.0E-3\n    RMS_DR 1.5E-3\n    MAX_FORCE 4.5E-4\n    RMS_FORCE 3.0E-4")
    shift=(26-box)/2
    lines=s.splitlines();coord=False;new=[]
    for line in lines:
        if line.strip()=="&COORD":coord=True;new.append(line);continue
        if line.strip()=="&END COORD":coord=False;new.append(line);continue
        if coord:
            a=line.split()
            if len(a)!=4: raise RuntimeError("Unexpected coordinate count")
            new.append("      "+a[0]+" "+" ".join(f"{float(v)-shift:.8f}" for v in a[1:]))
        else:new.append(line)
    s="\n".join(new)+"\n"
    s=s.replace("&END MOTION","""  &PRINT
    &TRAJECTORY
      FORMAT XYZ
      &EACH
        GEO_OPT 1
      &END EACH
    &END TRAJECTORY
  &END PRINT
&END MOTION""")
    return s

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="build/v6_dft")
    p.add_argument("--box",type=float,default=18)
    p.add_argument("--cutoff",type=int,default=220)
    p.add_argument("--max-iter",type=int,default=55)
    p.add_argument("--eps-scf",type=float,default=1e-5)
    a=p.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    for ion in ("Nd","Dy"):
        for n in (8,9):
            f=out/f"{ion.lower()}_w{n}_geo.inp"
            f.write_text(build(ion,n,a.box,a.cutoff,a.eps_scf,a.max_iter),encoding="utf8")
            print(f)
if __name__=="__main__":main()
