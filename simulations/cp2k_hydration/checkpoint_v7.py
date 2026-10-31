#!/usr/bin/env python3
"""CP2K 2026.2 checkpointed GEO_OPT and safe continuation from V6 XYZ/restart.

Official CP2K refs: https://manual.cp2k.org/trunk/methods/restarting.html
https://manual.cp2k.org/trunk/methods/optimization/geometry.html

A last complete XYZ frame restarts *positions only*. It does not restore optimizer
Hessian/history or an SCF wavefunction. Prefer real .restart + .wfn if available.
Do not interpret a successfully copied XYZ as a converged geometry.
"""
from __future__ import annotations
import argparse, hashlib, json, re, shutil
from pathlib import Path

def read_last_xyz(path, atom_count, ion):
    lines=Path(path).read_text(encoding="utf8").splitlines()
    current=None
    index=0
    while index<len(lines):
        try:
            n=int(lines[index].strip())
        except (ValueError,IndexError):
            break
        if n!=atom_count or index+2+n>len(lines):
            break
        rows=[]
        for line in lines[index+2:index+2+n]:
            tokens=line.split()
            if len(tokens)<4:
                raise ValueError("Malformed XYZ coordinate record")
            rows.append((tokens[0],*[float(s) for s in tokens[1:4]]))
        expect=[ion]+["O","H","H"]*((atom_count-1)//3)
        if [r[0] for r in rows] != expect:
            raise ValueError("XYZ species order mismatch")
        current=rows
        index+=n+2
    if current is None:
        raise ValueError("No complete, valid XYZ frame found")
    return current

def replace_coords(input_text, atoms):
    m=re.search(r"(?im)^(\s*)&COORD\s*\n(.*?)^\s*&END COORD\s*$",input_text,re.S|re.M)
    if not m:
        raise ValueError("CP2K &COORD section missing")
    original=[x.split()[0] for x in m.group(2).splitlines() if x.strip()]
    if original != [row[0] for row in atoms]:
        raise ValueError("Coordinate atom count/order mismatch with CP2K input")
    pos="\n".join(f"      {r[0]:<2s} {r[1]:.8f} {r[2]:.8f} {r[3]:.8f}" for r in atoms)
    return input_text[:m.start()]+m.group(1)+"&COORD\n"+pos+"\n"+m.group(1)+"&END COORD"+input_text[m.end():]

def inject_print_restart(s):
    if "&RESTART ON" in s:
        raise ValueError("Already has restart section; do not duplicate")
    # &DFT/&SCF/&PRINT/&RESTART saves wfn at each geometry step.
    scf_tag="    &END SCF"
    if s.count(scf_tag)!=1:raise ValueError("Unique SCF block not found")
    scf_print="""      &PRINT
        &RESTART ON
          BACKUP_COPIES 1
          &EACH
            GEO_OPT 1
            QS_SCF 0
          &END EACH
        &END RESTART
      &END PRINT
"""
    s=s.replace(scf_tag,scf_print+scf_tag,1)
    # &MOTION/&PRINT/&RESTART saves geometry/settings each step.
    motion_pos=s.find("&MOTION")
    if motion_pos<0:raise ValueError("MOTION section not found")
    print_pos=s.find("  &PRINT\n",motion_pos)
    if print_pos<0:raise ValueError("MOTION PRINT section not found")
    geom_print="""    &RESTART ON
      BACKUP_COPIES 2
      &EACH
        GEO_OPT 1
      &END EACH
    &END RESTART
"""
    s=s[:print_pos+len("  &PRINT\n")]+geom_print+s[print_pos+len("  &PRINT\n"):]
    return s

def prepare(input_path, output_path, xyz_path=None, restart_file=None, wfn_file=None):
    original=Path(input_path).read_text(encoding="utf8")
    ionmatch=re.search(r"PROJECT\s+((nd|dy)_w([89]))",original,re.I)
    if not ionmatch:raise ValueError("Cannot identify 8/9-water Nd/Dy CP2K cluster")
    ion=ionmatch.group(2).capitalize();n=int(ionmatch.group(3))
    atom_count=1+3*n
    if xyz_path and restart_file:
        raise ValueError("Use one coordinate source: XYZ fallback OR CP2K .restart")
    s=original
    source="original_initial_guess"
    if xyz_path:
        s=replace_coords(s,read_last_xyz(xyz_path,atom_count,ion))
        source="last_complete_xyz_positions_only"
    if restart_file:
        if not Path(restart_file).is_file():raise FileNotFoundError(restart_file)
        # May restore optimizer counters when using .restart as full input;
        # here deliberately restart coordinates only to keep current known settings.
        block="""&EXT_RESTART
  RESTART_FILE_NAME {file}
  RESTART_DEFAULT .FALSE.
  RESTART_POS .TRUE.
  RESTART_COUNTERS .FALSE.
&END EXT_RESTART
""".format(file=Path(restart_file).name)
        if "&EXT_RESTART" in s:raise ValueError("Input already contains EXT_RESTART")
        s=s.replace("&FORCE_EVAL",block+"&FORCE_EVAL",1)
        source="cp2k_restart_coordinates_only"
    s=inject_print_restart(s)
    if wfn_file:
        if not Path(wfn_file).is_file():raise FileNotFoundError(wfn_file)
        s=s.replace("    SCF_GUESS ATOMIC","    SCF_GUESS RESTART",1)
        if "SCF_GUESS RESTART" not in s:raise ValueError("Cannot set SCF_GUESS RESTART")
        s=s.replace("  &DFT\n","  &DFT\n    WFN_RESTART_FILE_NAME "+Path(wfn_file).name+"\n",1)
    out=Path(output_path);out.parent.mkdir(parents=True,exist_ok=True)
    for p in [restart_file,wfn_file]:
        if p:
            src=Path(p);dst=out.parent/src.name
            if src.resolve()!=dst.resolve():shutil.copy2(src,dst)
    out.write_text(s,encoding="utf8")
    provenance={
        "ion":ion,"nwater":n,"input":str(input_path),
        "coordinate_source":source,
        "xyz_path":str(xyz_path) if xyz_path else None,
        "cp2k_restart":Path(restart_file).name if restart_file else None,
        "wavefunction_restart":Path(wfn_file).name if wfn_file else None,
        "input_sha256":hashlib.sha256(s.encode()).hexdigest(),
        "warning":"A restart/geometry is not convergence; require final CP2K force and SCF validation."
    }
    out.with_suffix(".provenance.json").write_text(json.dumps(provenance,indent=2)+"\n")
    return provenance

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",required=True)
    p.add_argument("--output",required=True)
    p.add_argument("--xyz")
    p.add_argument("--restart")
    p.add_argument("--wfn")
    a=p.parse_args()
    print(json.dumps(prepare(a.input,a.output,a.xyz,a.restart,a.wfn),indent=2))
if __name__=="__main__":main()
