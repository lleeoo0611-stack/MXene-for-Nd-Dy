"""V7 regression tests: restart semantics and actual distinct CN8/CN9 water packing."""
from __future__ import annotations
import math,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"simulations/cp2k_hydration"))
sys.path.insert(0,str(ROOT/"simulations/openmm_hydration"))
from generate_v6 import build
from checkpoint_v7 import read_last_xyz,prepare
from cn_initial_v7 import water_shell,prepare_shell,minimum_image_norm

class Checkpoints(unittest.TestCase):
    def make_xyz(self,ion,n):
        content=build(ion,n)
        ls=content.split("&COORD",1)[1].split("&END COORD",1)[0].strip().splitlines()
        a=[x.split() for x in ls]
        def fmt(delta):
            return str(len(a))+"\nLast complete frame\n"+"\n".join(
                "%s %.7f %.7f %.7f"%(v[0],float(v[1])+delta,float(v[2]),float(v[3]))
                for v in a)+"\n"
        return fmt(0)+fmt(.012), len(a)

    def test_complete_xyz_restart_and_checkpoint_output(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            f=root/"dy_w9_geo.inp";f.write_text(build("Dy",9))
            xyz,n=self.make_xyz("Dy",9)
            z=root/"last.xyz";z.write_text(xyz+"28\nincomplete\nDy 1.0 1.0 1.0\n")
            rows=read_last_xyz(z,n,"Dy")
            self.assertEqual(len(rows),28)
            out=root/"continued.inp"
            meta=prepare(f,out,xyz_path=z)
            a=out.read_text()
            self.assertEqual(meta["coordinate_source"],"last_complete_xyz_positions_only")
            self.assertEqual(a.count("&EXT_RESTART"),0)
            self.assertIn("GEO_OPT 1",a)
            self.assertEqual(a.count("&RESTART ON"),2)
            self.assertIn("&SCF",a)
            self.assertEqual(a.count("&COORD"),1)
            self.assertTrue(out.with_suffix(".provenance.json").exists())

    def test_restart_and_wfn_not_silently_claimed(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            f=root/"nd_w8_geo.inp";f.write_text(build("Nd",8))
            restart=root/"nd_w8_geo-1.restart";restart.write_text("&GLOBAL\n&END GLOBAL\n")
            wfn=root/"nd_w8_geo-RESTART.wfn";wfn.write_bytes(b"dummy-for-test")
            out=root/"continue.inp"
            p=prepare(f,out,restart_file=restart,wfn_file=wfn)
            text=out.read_text()
            self.assertIn("&EXT_RESTART",text)
            self.assertIn("RESTART_POS .TRUE.",text)
            self.assertIn("SCF_GUESS RESTART",text)
            self.assertIn("WFN_RESTART_FILE_NAME",text)
            self.assertEqual(p["coordinate_source"],"cp2k_restart_coordinates_only")
            with self.assertRaises(ValueError):
                prepare(f,root/"bad.inp",xyz_path=f,restart_file=restart)

    def test_wrong_xyz_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"bad.xyz";p.write_text("3\nbad\nO 1 1 1\nH 1 1 2\nH 1 1 3\n")
            with self.assertRaises(ValueError):
                read_last_xyz(p,28,"Dy")

class HydrationInitializations(unittest.TestCase):
    def test_shell_OH_geometry_and_counts(self):
        c=(1.25,1.25,1.25)
        for n in (8,9):
            waters=water_shell(c,n,rotation=.35)
            self.assertEqual(len(waters),n)
            for o,h1,h2 in waters:
                self.assertAlmostEqual(math.dist(c,o),.245,places=8)
                self.assertAlmostEqual(math.dist(o,h1),.09572,places=7)
                self.assertAlmostEqual(math.dist(o,h2),.09572,places=7)
                self.assertAlmostEqual(math.dist(h1,h2),2*.09572*math.sin(math.radians(52.26)),places=6)

    def test_independent_solvent_packing_and_target_CN(self):
        # Unrelated bulk waters beyond exclusion sphere, real prep inserts inner shell.
        coords=[];idx=[]
        for i in range(30):
            o=(.17*(i%5)+.2, .17*((i//5)%3)+.2, .17*(i//15)+.2)
            h1=(o[0]+.09572,o[1],o[2]);h2=(o[0]-.02,o[1]+.09,o[2])
            k=len(coords);coords.extend([o,h1,h2]);idx.append((k,k+1,k+2))
        result=[]
        for cn in (8,9):
            ws,c,meta=prepare_shell(coords,idx,2.5,4312+cn,cn)
            self.assertEqual(meta["initial_CN_0p35nm"],cn)
            self.assertGreaterEqual(len(ws),cn)
            self.assertAlmostEqual(sum(minimum_image_norm(w[0],c,2.5)<.35 for w in ws),cn)
            result.append((len(ws),meta["packing_shift_nm"]))
        self.assertNotEqual(result[0][1],result[1][1])

if __name__=="__main__":
    unittest.main()
