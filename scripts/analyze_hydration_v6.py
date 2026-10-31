#!/usr/bin/env python3
"""Aggregate three independent V6 OpenMM water benchmarks. No fabricated data.

Input root must contain real extracted GitHub Actions artifact directories
(each containing nd/dy_summary.json, *_rdf.csv, *_cn_timeseries.csv).
Scientific validation is strict; short 40+80 ps runs are PRELIMINARY only.
"""
from __future__ import annotations
import argparse,csv,json,math,statistics
from pathlib import Path

def analyze(root):
    paths=sorted(Path(root).rglob("*_summary.json"))
    results=[]
    for path in paths:
        x=json.loads(path.read_text())
        rdf=path.parent/(x["ion"].lower()+"_rdf.csv")
        if not rdf.exists():raise ValueError("No RDF histogram for "+str(path))
        with rdf.open(newline="") as fh:
            bins=list(csv.DictReader(fh))
        cutoff=x.get("RDF_first_min_nm") or .35
        pairs=[(float(r["r_nm"]),int(r["raw_count"])) for r in bins if float(r["r_nm"])<cutoff]
        count=sum(n for _,n in pairs)
        x["firstshell_mean_LnO_A_from_hist"]=(10*sum(r*n for r,n in pairs)/count if count else None)
        x["firstshell_CN_from_RDF_hist"]=(count/x["frames"] if count else None)
        results.append(x)
    by={}
    for entry in results:
        ion=entry["ion"];seed=entry.get("seed")
        if seed is None:raise ValueError("Missing seed in "+ion)
        by.setdefault(ion,{})
        if seed in by[ion]:raise ValueError("Duplicate "+ion+str(seed))
        by[ion][seed]=entry
    rows=[]
    for ion in ("Nd","Dy"):
        runs=list(by.get(ion,{}).values())
        if not runs:continue
        cn=[r["observed_mean_CN_r0p35nm"] for r in runs]
        peak=[r["observed_first_rdf_peak_nm"]*10 for r in runs]
        mean_distance=[r["firstshell_mean_LnO_A_from_hist"] for r in runs if r["firstshell_mean_LnO_A_from_hist"] is not None]
        firstmin=[r["CN_mean_RDF_first_min"] for r in runs
                  if r.get("CN_mean_RDF_first_min") is not None]
        def average(x):return statistics.mean(x) if x else None
        def sem(x):return statistics.stdev(x)/math.sqrt(len(x)) if len(x)>1 else None
        physics=dict(Nd={"EXAFS_distance_A":2.51,"main_CN_expected":9},
                     Dy={"EXAFS_distance_A":2.39,"main_CN_expected":8})[ion]
        min_eq=min(r["equilibration_ps"] for r in runs)
        min_prod=min(r["production_ps"] for r in runs)
        sufficient=(len(runs)>=3 and min_eq>=500 and min_prod>=2000
                    and len(firstmin)==len(runs))
        rows.append(dict(ion=ion,n_replicas=len(runs),seeds=";".join(map(str,sorted(by[ion]))),
            min_equilibration_ps=min_eq,min_production_ps=min_prod,
            mean_RDF_first_peak_A=average(peak),sem_between_replicas_RDF_peak_A=sem(peak),
            mean_firstshell_LnO_A_from_RDF=average(mean_distance),
            sem_firstshell_LnO_A_between_replicas=sem(mean_distance),
            mean_CN_0p35nm=average(cn),sem_between_replicas_CN_0p35nm=sem(cn),
            mean_CN_RDF_first_min=average(firstmin),first_min_available_replicas=len(firstmin),
            literature_EXAFS_first_shell_distance_A=physics["EXAFS_distance_A"],
            literature_dominant_CN=physics["main_CN_expected"],
            mean_firstshell_distance_deviation_from_EXAFS_A=average(mean_distance)-physics["EXAFS_distance_A"],
            validated_vs_literature=False,
            reason=("duration/repeats/forcefield/finite-size checks incomplete"
                    if not sufficient else "requires further independent forcefield and box convergence")))
    return rows

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--root",required=True)
    a.add_argument("--out",default="build/v6_aggregate")
    x=a.parse_args()
    rows=analyze(x.root)
    if not rows:raise SystemExit("No MD summaries found; cannot claim completed.")
    out=Path(x.out);out.mkdir(parents=True,exist_ok=True)
    with (out/"summary.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    (out/"qc.json").write_text(json.dumps(rows,indent=2)+"\n")
    for r in rows:
        print(r)
        if not r["validated_vs_literature"]:
            print("NOT VALIDATED: "+r["ion"]+" requires long trajectories and controls.")

if __name__=="__main__":main()
