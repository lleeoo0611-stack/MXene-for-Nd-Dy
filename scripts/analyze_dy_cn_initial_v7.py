#!/usr/bin/env python3
"""Compare actually sampled CN8/CN9-initialized Dy OpenMM trajectories.

Requires extracted Actions artifacts in separate directories. This script DOES
NOT compute an equilibrium probability from correlated short trajectories.
"""
import argparse,csv,json,statistics
from pathlib import Path

def aggregate(root):
    values=[]
    for summary_file in sorted(Path(root).rglob("dy_summary.json")):
        summary=json.loads(summary_file.read_text())
        initial=summary.get("initial_cn_requested")
        if initial not in (8,9):
            continue
        fn=summary_file.with_name("dy_cn_timeseries.csv")
        if not fn.exists():
            raise RuntimeError("Missing observed CN trajectory: "+str(fn))
        with fn.open(newline="") as f:
            series=list(csv.DictReader(f))
        counts=[]
        for x in series:
            try:
                counts.append(int(x["CN_at_RDF_first_min"]))
            except (TypeError,ValueError):
                counts.append(int(float(x["CN_within_0p35nm"])))
        if len(counts)!=summary["frames"]:
            raise RuntimeError("Frame count mismatch")
        # Correlated frames: observed fraction is a trajectory diagnostic only.
        row={
            "initial_CN":initial,"seed":summary["seed"],
            "water_count":summary["n_water"],
            "CN_before_minimization":summary["initial_cn_before_minimization"],
            "CN_after_minimization":summary["cn_after_minimization"],
            "CN_first_sampled":counts[0],
            "CN_last_sampled":counts[-1],
            "CN_mean":statistics.mean(counts),
            "CN8_fraction_of_correlated_frames":counts.count(8)/len(counts),
            "CN9_fraction_of_correlated_frames":counts.count(9)/len(counts),
            "changes_between_sampled_frames":sum(x!=y for x,y in zip(counts,counts[1:])),
            "equilibration_ps":summary["equilibration_ps"],
            "production_ps":summary["production_ps"],
            "frames":len(counts),
            "validated_equilibrium":False,
        }
        values.append(row)
    if not values:
        raise ValueError("No valid Dy CN8/CN9 initial-shell summaries")
    return values

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--artifacts-dir",required=True)
    p.add_argument("--output",default="build/v7_comparison.csv")
    a=p.parse_args()
    rows=aggregate(a.artifacts_dir)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    with dst.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    print("Observed trajectories (NOT independent equilibrium samples):")
    for row in rows:print(row)
    cn8={r["initial_CN"] for r in rows}
    if cn8!={8,9}:print("WARNING: Both CN8 and CN9 starts not yet available.")
    print("NO thermodynamic validation: longer sampling, transitions, box and force-field tests are still required.")

if __name__=="__main__":
    main()
