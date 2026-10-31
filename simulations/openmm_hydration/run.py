#!/usr/bin/env python3
"""Short classical MD proof-of-execution: Nd3+ / Dy3+ in water, 12-6-4.

This is NOT hydrated Ln in citrate-MXene, not an equilibrated CN estimate, and
NOT Duvail 2024 parametrization.

Parameters from Li, Song & Merz, J Phys Chem B 119 (2015) 883-895,
DOI 10.1021/jp505875v: first, TIP3P-compatible 12-6-4 column.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

METALS = {
    "Nd": dict(mass=144.242, rmin_half_A=1.724,
               epsilon_kcal=0.15486311, c4_kcal_A4=213),
    "Dy": dict(mass=162.500, rmin_half_A=1.637,
               epsilon_kcal=0.09900804, c4_kcal_A4=207),
}

def ion_parameters(ion):
    x=METALS[ion]
    return dict(mass=x["mass"], sigma_nm=.2*x["rmin_half_A"]/(2**(1/6)),
                epsilon_kJmol=4.184*x["epsilon_kcal"],
                c4_kJmol_nm4=4.184e-4*x["c4_kcal_A4"])

def run(ion, output_dir, equil_steps=2000, sample_steps=4000,
        sample_interval=25, box_nm=2.5, seed=20261031):
    import numpy as np
    import openmm
    from openmm import Vec3, unit, app
    if not box_nm>2.0:raise ValueError("Periodic box must exceed 2*0.9 nm")
    if sample_steps<sample_interval:raise ValueError("No MD frames")
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    mod=app.Modeller(app.Topology(),[])
    ff=app.ForceField("amber14/tip3p.xml")
    mod.addSolvent(ff,model="tip3p",boxSize=Vec3(box_nm,box_nm,box_nm)*unit.nanometer)
    xyz=np.array([[p.x,p.y,p.z] for p in mod.positions.value_in_unit(unit.nanometer)])
    waters=[]
    for res in mod.topology.residues():
        atoms=list(res.atoms())
        if len(atoms)!=3:raise ValueError("Unexpected solvent topology")
        oxygen=next(a.index for a in atoms if a.element.symbol=="O")
        hydrogens=[a.index for a in atoms if a.element.symbol=="H"]
        waters.append((oxygen,hydrogens[0],hydrogens[1]))
    # Replace the water closest to cell center with one Ln3+.
    idx=min(range(len(waters)),
            key=lambda k:float(np.linalg.norm(xyz[waters[k][0]]-box_nm/2)))
    metal_pos=xyz[waters[idx][0]]
    waters.pop(idx)

    system=openmm.System()
    v=Vec3
    system.setDefaultPeriodicBoxVectors(v(box_nm,0,0)*unit.nanometer,
        v(0,box_nm,0)*unit.nanometer,v(0,0,box_nm)*unit.nanometer)
    nb=openmm.NonbondedForce()
    nb.setNonbondedMethod(openmm.NonbondedForce.PME)
    nb.setCutoffDistance(.9*unit.nanometer)
    nb.setEwaldErrorTolerance(5e-4)
    pos=[]; oxygen_indices=[]
    for i,j,k in waters:
        a=system.addParticle(15.999*unit.dalton)
        b=system.addParticle(1.008*unit.dalton)
        c=system.addParticle(1.008*unit.dalton)
        oxygen_indices.append(a)
        pos.extend((xyz[i],xyz[j],xyz[k]))
        nb.addParticle(-.834,.315075*unit.nanometer,
                       .636386*unit.kilojoule_per_mole)
        nb.addParticle(+.417,.1*unit.nanometer,0*unit.kilojoule_per_mole)
        nb.addParticle(+.417,.1*unit.nanometer,0*unit.kilojoule_per_mole)
        system.addConstraint(a,b,.09572*unit.nanometer)
        system.addConstraint(a,c,.09572*unit.nanometer)
        system.addConstraint(b,c,.15139*unit.nanometer)
        for u,w in ((a,b),(a,c),(b,c)):
            nb.addException(u,w,0,.1*unit.nanometer,
                            0*unit.kilojoule_per_mole)
    p=ion_parameters(ion)
    mid=system.addParticle(p["mass"]*unit.dalton)
    pos.append(metal_pos)
    nb.addParticle(+3,p["sigma_nm"]*unit.nanometer,
                   p["epsilon_kJmol"]*unit.kilojoule_per_mole)
    system.addForce(nb)

    # -C4/r^4 between Ln and water O ONLY, 12-6 and PME supplied by nb.
    induction=openmm.CustomNonbondedForce("-C4/r^4")
    induction.addGlobalParameter("C4",p["c4_kJmol_nm4"])
    for _ in range(system.getNumParticles()):
        induction.addParticle([])
    induction.addInteractionGroup({mid},set(oxygen_indices))
    induction.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    induction.setCutoffDistance(.9*unit.nanometer)
    system.addForce(induction)

    integ=openmm.LangevinMiddleIntegrator(
        298.15*unit.kelvin,1/unit.picosecond,.002*unit.picosecond)
    integ.setRandomNumberSeed(seed)
    ctx=openmm.Context(system,integ,openmm.Platform.getPlatformByName("CPU"),
                       {"Threads":"2"})
    ctx.setPositions(np.asarray(pos)*unit.nanometer)
    openmm.LocalEnergyMinimizer.minimize(ctx,maxIterations=200)
    ctx.setVelocitiesToTemperature(298.15*unit.kelvin,seed)
    integ.step(equil_steps)

    edges=np.linspace(0,1,101)
    hist=np.zeros(100)
    coord=[]
    frames=sample_steps//sample_interval
    for n in range(frames):
        integ.step(sample_interval)
        a=ctx.getState(getPositions=True).getPositions(
            asNumpy=True).value_in_unit(unit.nanometer)
        dr=a[oxygen_indices]-a[mid]
        dr-=box_nm*np.rint(dr/box_nm)
        r=np.linalg.norm(dr,axis=1)
        hist+=np.histogram(r,bins=edges)[0]
        coord.append(float(np.sum(r<.35)))
    vol=4*np.pi/3*(edges[1:]**3-edges[:-1]**3)
    g=hist/(frames*vol*(len(waters)/box_nm**3))
    window=(edges[:-1]>=.20)&(edges[1:]<=.40)
    peak=int(np.argmax(np.where(window,g,-1)))
    with (out/(ion.lower()+"_rdf.csv")).open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["r_nm","g_LnO","raw_count"])
        for i in range(len(g)):
            w.writerow([float((edges[i]+edges[i+1])/2),
                        float(g[i]),int(hist[i])])
    summary={
      "status":"SHORT_MD_SHAKEDOWN_NOT_EQUILIBRATED",
      "ion":ion,"doi_parameters":"10.1021/jp505875v",
      "parameters":p,"water_model":"TIP3P",
      "n_water":len(waters),"box_nm":box_nm,"net_charge_e":3,
      "PME_background":"implicit uniform neutralizing plasma",
      "equilibration_ps":equil_steps*.002,
      "production_ps":frames*sample_interval*.002,
      "frames":frames,
      "observed_mean_CN_r0p35nm":float(np.mean(coord)),
      "observed_first_rdf_peak_nm":float((edges[peak]+edges[peak+1])/2),
      "openmm_version":openmm.version.version,
      "warnings":[
        "Short run; coordination number and RDF not statistically converged",
        "No MXene/citrate or counterions; not transport or hydration free energy",
        "Published TIP3P Li-Merz parameters, not Duvail 2024 OPC parameter set",
        "Fixed net +3 neutralizing PME background causes finite-size artifacts"
      ]
    }
    (out/(ion.lower()+"_summary.json")).write_text(
        json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
    del ctx
    return summary

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--ion",choices=["Nd","Dy"],required=True)
    p.add_argument("--out",default="build/md")
    p.add_argument("--equil-steps",type=int,default=2000)
    p.add_argument("--sample-steps",type=int,default=4000)
    p.add_argument("--interval",type=int,default=25)
    a=p.parse_args()
    run(a.ion,a.out,a.equil_steps,a.sample_steps,a.interval)
