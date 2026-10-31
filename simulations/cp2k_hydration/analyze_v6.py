"""Parse CP2K GEO_OPT stdout/XYZ and refuse to label failed jobs converged."""
from __future__ import annotations
import argparse,csv,json,math,re
from pathlib import Path

def analyze(log, workdir):
    logfile=Path(log)
    text=logfile.read_text(errors="replace") if logfile.exists() else ""
    match=re.search(r"(nd|dy)_w([89])",str(logfile).lower())
    ion=match.group(1) if match else "unknown"
    nw=int(match.group(2)) if match else 0
    completed=("GEOMETRY OPTIMIZATION COMPLETED" in text)
    failed_scf=bool(re.search(r"SCF run NOT converged|SCF did NOT converge",text,re.I))
    entries=re.findall(r"ENERGY\| Total FORCE_EVAL[^\n]*?([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)\s*$",text,re.M)
    energy=float(entries[-1]) if entries else None
    candidates=list(Path(workdir).glob("*-pos-1.xyz"))
    last=None
    if candidates:
        xyz=max(candidates,key=lambda p:p.stat().st_mtime).read_text().splitlines()
        i=0
        while i<len(xyz):
            try:n=int(xyz[i])
            except (ValueError,IndexError):break
            data=xyz[i+2:i+2+n]
            if len(data)<n:break
            last=[(p[0],tuple(float(v) for v in p[1:4])) for p in (line.split() for line in data)]
            i+=n+2
    distances=[]
    if last:
        metals=[p for el,p in last if el in ("Nd","Dy")]
        oxy=[p for el,p in last if el=="O"]
        if metals:
            distances=sorted(math.dist(p,metals[0]) for p in oxy)
    return {"ion":ion,"nwater":nw,"optimization_completed":completed,
            "scf_failure_reported":failed_scf,
            "energy_hartree":energy,
            "LnO_distances_A":distances,
            "mean_LnO_A":sum(distances)/len(distances) if distances else None,
            "LnO_within_3p5A":sum(d<3.5 for d in distances),
            "electronic_energy_evaluations":len(entries),
            "log":str(logfile)}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--log",required=True)
    p.add_argument("--directory",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    r=analyze(a.log,a.directory)
    Path(a.output).write_text(json.dumps(r,indent=2)+"\n")
    print(json.dumps(r,indent=2))
    if not r["optimization_completed"] or r["scf_failure_reported"]:
        raise SystemExit("CP2K geometry optimization NOT verified as converged.")

if __name__=="__main__":main()
