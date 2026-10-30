#!/usr/bin/env python3
"""Literature-informed screening: validated equations; surface numbers are SCENARIOS.

Do not interpret site K, site pKa, transport barriers or t50 as measured Nd/Dy
quantities. Does not implement soluble Ln-citrate complexation or true MD/DFT.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import brentq

BULK_PKA=(3.13,4.76,6.40) # free citric acid; ACS Omega 2019, 10.1021/acsomega.9b00270
RT=0.008314462618*298.15

def citrate_species(pH,pka=BULK_PKA):
    """fractions: H3Cit, H2Cit-, HCit2-, Cit3-, in absence of metals."""
    h=10.**(-np.asarray(pH,dtype=float))
    a,b,c=(10.**(-x) for x in pka)
    w=np.array([h**3,a*h*h,a*b*h,np.ones_like(h)*a*b*c])
    return w/w.sum(axis=0)

def site_deprot(pH,pka_surface=5.5):
    """Unknown grafted surface pKa=5.5 is SCENARIO, not bulk citrate pKa."""
    return 1/(1+10.**(pka_surface-np.asarray(pH,dtype=float)))

def competitive_occupancy(cNd,cDy,Knd,Kdy,f_site=1.):
    x=Knd*cNd*f_site
    y=Kdy*cDy*f_site
    return x/(1+x+y),y/(1+x+y),1/(1+x+y)

def normalized_downstream_flux(theta):
    """J(t)/Jsteady for 1D diffusion with instantaneously reversible fixed binding,
    empty membrane, constant upstream mobile concentration, downstream sink.
    theta=D0*t/(R*L**2); R=1+S*K in dilute linear binding limit.
    """
    t=np.asarray(theta,dtype=float)
    if np.any(t<=0):raise ValueError('theta must be positive')
    n=np.arange(1,101).reshape((-1,)+(1,)*t.ndim)
    val=1+2*np.sum((-1.)**n*np.exp(-n*n*np.pi*np.pi*t),axis=0)
    return np.clip(val,0,1)

def breakthrough_t50_dimensionless():
    return brentq(lambda z:float(normalized_downstream_flux(z))-.5,1e-3,10)

def permeability_ratio(delta_G_entry_Dy_minus_Nd_kJmol=0.,Dnd_over_Ddy=1.):
    """P_Nd/P_Dy=exp((Gentry_Dy-Gentry_Nd)/RT) * D_Nd/D_Dy.
    Energy and diffusion ratio are NOT known for Nd/Dy citrate MXene.
    """
    return np.exp(np.asarray(delta_G_entry_Dy_minus_Nd_kJmol)/RT)*Dnd_over_Ddy

def run(out='results/literature_v2'):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    phs=np.linspace(3,8,101)
    alpha=citrate_species(phs)
    pd.DataFrame(dict(pH=phs,H3Cit=alpha[0],H2Cit_minus=alpha[1],
          HCit_2minus=alpha[2],Cit_3minus=alpha[3])).to_csv(out/'citrate_uncomplexed.csv',index=False)
    fig,ax=plt.subplots()
    for y,name in zip(alpha,('H3Cit','H2Cit-','HCit2-','Cit3-')):
        ax.plot(phs,y,label=name)
    ax.set(xlabel='pH',ylabel='Uncomplexed citrate fraction',title='Bulk citrate acid/base speciation')
    ax.legend();fig.tight_layout();fig.savefig(out/'citrate_uncomplexed.png',dpi=160);plt.close(fig)
    # Scenario: reversible fixed sites alter transient breakthrough, not steady flux.
    time=np.geomspace(.005,10,250)
    dat={'dimensionless_t_D0_over_L2':time}
    fig,ax=plt.subplots()
    for R in [1,5,20,100]:
        f=normalized_downstream_flux(time/R)
        dat[f'J_over_Jsteady_retardation_{R}']=f
        ax.semilogx(time,f,label=f'R={R}')
    ax.set(xlabel='Dimensionless time t D0/L²',ylabel='J(t)/Jsteady',
           title='Fixed reversible binding: transient response')
    ax.legend();fig.tight_layout();fig.savefig(out/'breakthrough.png',dpi=160);plt.close(fig)
    pd.DataFrame(dat).to_csv(out/'breakthrough.csv',index=False)
    # Unknown relative site affinity and grafted site pKa; condition scan only.
    rows=[]
    unit_t50=breakthrough_t50_dimensionless()
    for ph in [4.5,5.5,6.5]:
        f=float(site_deprot(ph)) # pKa_surface ASSUMED
        for delta in np.linspace(-1,1,41):
            Rnd=1+20*f        # S*K_Nd at full site ionization ASSUMED=20
            Rdy=1+20*f*10**delta
            rows.append(dict(pH=ph,log10_Ksite_Dy_over_Nd_ASSUMED=delta,
                retardation_Nd=Rnd,retardation_Dy=Rdy,
                t50_Nd_D0_over_L2=unit_t50*Rnd,
                t50_Dy_D0_over_L2=unit_t50*Rdy,
                t50_Dy_over_Nd=Rdy/Rnd,
                steady_P_Nd_over_Dy_if_equal_mobile_K_D=1.))
    pd.DataFrame(rows).to_csv(out/'binding_sensitivity.csv',index=False)
    energies=np.linspace(-4,4,81)
    pd.DataFrame({'Gentry_Dy_minus_Nd_kJmol_ASSUMED':energies,
      **{f'Pnd_over_Pdy_if_Dnd_over_Ddy_{d:g}':permeability_ratio(energies,d)
         for d in [0.5,1.,2.]}}).to_csv(out/'permeability_sensitivity.csv',index=False)
    print(f'V2 written to {out}; t50 * D0/(R L²) = {unit_t50:.6f}')
    for p in [4.5,5.5,6.5]:
        print('pH',p,'Cit3- fraction (uncomplexed):',round(float(citrate_species(p)[3]),5))
    print('Binding-only retardation does not predict steady-state permeability selectivity.')

if __name__=='__main__':
    run()
