#!/usr/bin/env python3
"""Build genuinely distinct Dy3+ coordination-eight/nine initial water shells.

Standalone pure-Python geometry helper: no OpenMM/NumPy import needed.
The inner-shell H atoms point away from +3 central ion.
Shell creation is a deliberate out-of-equilibrium starting condition, NOT
force-field equilibration or a prediction of eight- versus nine-coordination.
"""
import math,random

def unit(v):
    n=math.sqrt(sum(t*t for t in v))
    if n<=1e-12:raise ValueError("zero vector")
    return tuple(x/n for x in v)

def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])

def directions(n,angle=0.0):
    if n==8:
        points=[unit((a,b,c)) for a in (-1,1) for b in (-1,1) for c in (-1,1)]
    elif n==9:
        ga=math.pi*(3-math.sqrt(5))
        points=[]
        for i in range(n):
            z=1-2*(i+.5)/n
            r=math.sqrt(1-z*z)
            points.append((r*math.cos(i*ga),r*math.sin(i*ga),z))
    else:
        raise ValueError("Only CN 8 or CN 9 supported")
    c,s=math.cos(angle),math.sin(angle)
    return [(c*x-s*y,s*x+c*y,z) for x,y,z in points]

def minimum_image_norm(point,center,box_nm):
    delta=[a-b for a,b in zip(point,center)]
    delta=[x-box_nm*round(x/box_nm) for x in delta]
    return math.sqrt(sum(x*x for x in delta))

def water_shell(center,n,radius_nm=.245,rotation=0.):
    if n not in (8,9):raise ValueError(n)
    theta=math.radians(52.26)
    sin,cos=math.sin(theta),math.cos(theta)
    out=[]
    for u in directions(n,rotation):
        ox=tuple(center[k]+radius_nm*u[k] for k in range(3))
        axis=unit(cross(u,(0,0,1) if abs(u[2])<.85 else (0,1,0)))
        h1=tuple(ox[k]+.09572*(cos*u[k]+sin*axis[k]) for k in range(3))
        h2=tuple(ox[k]+.09572*(cos*u[k]-sin*axis[k]) for k in range(3))
        out.append((ox,h1,h2))
    return out

def prepare_shell(coords,waters,box_nm,seed,target_cn,
                  exclusion_radius_nm=.52,shell_radius_nm=.245):
    """coords is a 3N×3 list/array, waters index triples (O,H,H).

    Generate independently translated solvent packing, delete nearby bulk waters,
    insert exactly CN water molecules with outward H orientation. All non-core
    solvent molecules remain the same chemistry; missing deleted waters must be
    accounted for in the reported density and finite-size validation.
    """
    if target_cn not in (8,9):raise ValueError("Use --initial-cn 8 or 9")
    if box_nm<2.0:raise ValueError("box too small")
    if not .4<exclusion_radius_nm<.7:raise ValueError("unphysical exclusion radius")
    rng=random.Random(seed)
    shift=[rng.uniform(0,box_nm) for _ in range(3)]
    center=[box_nm/2]*3
    keep=[]
    removed=0
    for i,j,k in waters:
        original_o=list(map(float,coords[i]))
        # Translate whole water as a rigid entity and wrap oxygen only.
        translated_o=[(original_o[a]+shift[a])%box_nm for a in range(3)]
        shift_water=[translated_o[a]-original_o[a] for a in range(3)]
        translated=tuple(tuple(float(coords[idx][a])+shift_water[a] for a in range(3))
                         for idx in (i,j,k))
        if minimum_image_norm(translated_o,center,box_nm)<exclusion_radius_nm:
            removed+=1
        else:
            keep.append(translated)
    shell=water_shell(center,target_cn,shell_radius_nm,2*math.pi*rng.random())
    composed=keep+shell
    o_dist=[minimum_image_norm(w[0],center,box_nm) for w in composed]
    initial_cn=sum(d<.35 for d in o_dist)
    if initial_cn!=target_cn:
        raise RuntimeError("Initial shell compromised by bulk water placement")
    return composed,tuple(center),dict(
        seed=seed,coordination_target=target_cn,
        initial_CN_0p35nm=initial_cn,
        waters_removed_from_solvent=removed,
        water_count_final=len(composed),
        exclusion_radius_nm=exclusion_radius_nm,
        initial_shell_radius_nm=shell_radius_nm,
        packing_shift_nm=shift,
        note="Prepared metastable shell; not an equilibrated CN prediction")
