#!/usr/bin/env python3
"""Discrete-site conditional exit-time toy model. NOT physical ion permeability."""
from pathlib import Path
import numpy as np
import pandas as pd

RT=0.008314462618*298.15

def mfpt(g,intrinsic_barrier=3.,attempt_rate=1.):
    """Continuous-time Markov chain, detailed-balance-compatible step rates.
    Site 0: channel already entered; left reflecting. Last site: absorbing exit.
    Time unit is arbitrary unless attempt_rate is calibrated.
    """
    g=np.asarray(g,dtype=float)
    if len(g)<3 or attempt_rate<=0 or intrinsic_barrier<0:
        raise ValueError('Invalid model parameters')
    n=len(g)-1
    q=np.zeros((n,n))
    for i in range(n):
        for j in (i-1,i+1):
            if j<0 or j>n:continue
            saddle=max(g[i],g[j])+intrinsic_barrier
            rate=attempt_rate*np.exp(-(saddle-g[i])/RT)
            q[i,i]-=rate
            if j<n:q[i,j]+=rate
    return float(np.linalg.solve(-q,np.ones(n))[0])

def main():
    output=Path('results/toy_transport_conditional_mfpt.csv')
    output.parent.mkdir(exist_ok=True,parents=True)
    out=[]
    for nd_well,dy_well in [(0,0),(2,2),(2,4),(4,2),(4,6),(6,4)]:
        a=np.zeros(11);b=np.zeros(11)
        a[4]=-nd_well;b[4]=-dy_well
        tn=mfpt(a);td=mfpt(b)
        out.append(dict(Nd_assumed_well_depth_kJmol=nd_well,Dy_assumed_well_depth_kJmol=dy_well,
                        Nd_conditional_MFPT_arbitrary_unit=tn,Dy_conditional_MFPT_arbitrary_unit=td,
                        Nd_over_Dy_conditional_exit_rate=td/tn))
    d=pd.DataFrame(out)
    d.to_csv(output,index=False)
    print(d.to_string(index=False))
    print('No entry partition, true water coordination, membrane flux, or citrate force field modeled.')

if __name__=='__main__':main()
