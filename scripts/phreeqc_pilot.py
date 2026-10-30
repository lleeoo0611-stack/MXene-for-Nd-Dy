#!/usr/bin/env python3
"""72-case aqueous Nd/Dy speciation grid; excludes citrate, carbonate, membranes."""
import argparse, csv, subprocess, tempfile
from pathlib import Path
import pandas as pd

PHS=[4.0+i*0.5 for i in range(9)]
SALTS=[10.,100.]
METALS=[0.05,0.1,0.5,1.0]

def input_text(ph,salt,ln):
    return f"""TITLE Nd Dy NO-citrate NO-carbonate aqueous baseline
SOLUTION 1
 temp 25
 pH {ph:.2f}
 pe 4
 units mmol/kgw
 Na {salt:.6f}
 N(5) {salt+6*ln:.6f}
 Nd {ln:.6f}
 Dy {ln:.6f}
 -water 1
SELECTED_OUTPUT
 -file selected.tsv
 -reset false
 -pH true
 -ionic_strength true
 -totals Nd Dy
 -molalities Nd+3 Dy+3 NdNO3+2 DyNO3+2 NdOH+2 DyOH+2
 -saturation_indices Nd(OH)3 Dy(OH)3 Nd(OH)3(am) Dy(OH)3(am)
END
"""

def calculate(exe,db,output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    recs=[]
    with tempfile.TemporaryDirectory() as td:
        wd=Path(td)
        for salt in SALTS:
            for ln in METALS:
                for ph in PHS:
                    (wd/'selected.tsv').unlink(missing_ok=True)
                    (wd/'case.pqi').write_text(input_text(ph,salt,ln),encoding='utf-8')
                    r=subprocess.run([str(Path(exe).resolve()),'case.pqi','case.out',str(Path(db).resolve())],
                        cwd=wd,text=True,capture_output=True)
                    if r.returncode: raise RuntimeError(f'PHREEQC failed: {r.stderr}; {r.stdout}')
                    with (wd/'selected.tsv').open() as f:
                        row=next(csv.DictReader(f,delimiter='\t'))
                    v={k.strip():float(x) for k,x in row.items() if k and x and x.strip()}
                    def frac(spec,element):return v.get('m_'+spec,0.)/v[element]
                    def si(phase):
                        x=v.get('si_'+phase,float('nan'))
                        return float('nan') if x<=-999 else x
                    recs.append(dict(pH=ph,NaNO3_background_mmolkgw=salt,each_Ln_mmolkgw=ln,
                        Nd_free_fraction=frac('Nd+3','Nd'),Dy_free_fraction=frac('Dy+3','Dy'),
                        Nd_nitrate_fraction=frac('NdNO3+2','Nd'),Dy_nitrate_fraction=frac('DyNO3+2','Dy'),
                        SI_DyOH3_amorphous=si('Dy(OH)3(am)'),SI_DyOH3_crystalline=si('Dy(OH)3'),
                        SI_NdOH3_amorphous=si('Nd(OH)3(am)'),SI_NdOH3_crystalline=si('Nd(OH)3'),
                        ionic_strength_molkgw=v['mu']))
    df=pd.DataFrame(recs)
    assert len(df)==72
    df.to_csv(output/'phreeqc_pilot_grid.csv',index=False)
    print('Solved',len(df),'cases. Results:',output/'phreeqc_pilot_grid.csv')
    print(df[(df.pH==6.5)&(df.each_Ln_mmolkgw==0.1)][['NaNO3_background_mmolkgw','Nd_free_fraction','Dy_free_fraction','SI_DyOH3_amorphous']].to_string(index=False))
    return df

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--phreeqc',required=True,help='PHREEQC 3.8.6 executable')
    p.add_argument('--database',required=True,help='Matching llnl.dat')
    p.add_argument('--output',default='results')
    a=p.parse_args()
    calculate(a.phreeqc,a.database,a.output)

if __name__=='__main__':main()
