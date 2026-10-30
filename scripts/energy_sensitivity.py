#!/usr/bin/env python3
"""Hypothetical equal-prefactor activation barrier sensitivity. NOT DFT/MD/PMF."""
import argparse
import numpy as np
import pandas as pd

R=0.008314462618  # kJ mol^-1 K^-1

def rate_ratio(delta_barrier_dy_minus_nd,t=298.15):
    if t<=0:raise ValueError('temperature must be positive')
    return np.exp(np.asarray(delta_barrier_dy_minus_nd)/(R*t))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--temperature',type=float,default=298.15)
    p.add_argument('--output',default='results/energy_difference_selectivity.csv')
    a=p.parse_args()
    from pathlib import Path
    out=Path(a.output);out.parent.mkdir(exist_ok=True,parents=True)
    diffs=np.arange(-5,5.0001,0.25)
    pd.DataFrame({'DeltaGbarrier_Dy_minus_Nd_kJmol':diffs,
                  'hypothetical_relative_rate_Nd_over_Dy':rate_ratio(diffs,a.temperature)}).to_csv(out,index=False)
    print('Barrier for ratio 2:',float(R*a.temperature*np.log(2)),'kJ/mol; NOTE: synthetic sensitivity only')

if __name__=='__main__':main()
