import sys, unittest
from pathlib import Path
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/"simulations/cp2k_hydration"))
from generate_v6 import build
from analyze_v6 import analyze

class V6Tests(unittest.TestCase):
    def test_four_systems_and_coordinate_count(self):
        for ion in ("Nd","Dy"):
            for n in (8,9):
                s=build(ion,n)
                self.assertIn("RUN_TYPE GEO_OPT",s)
                self.assertIn("CUTOFF 220",s)
                self.assertIn("MAX_FORCE 4.5E-4",s)
                self.assertIn("MAX_DR 3.0E-3",s)
                self.assertIn("ABC 18.0 18.0 18.0",s)
                geom=s.split("&COORD\n",1)[1].split("&END COORD",1)[0].strip().splitlines()
                self.assertEqual(len(geom),1+3*n)
                self.assertEqual(geom[0].split()[0],ion)
    def test_optimization_not_claimed_without_log(self):
        d=analyze("/missing/dy_w9_geo.out","/missing")
        self.assertFalse(d["optimization_completed"])
    def test_different_box_is_centered(self):
        s=build("Nd",8,box=20)
        x=s.split("&COORD\n",1)[1].strip().splitlines()[0].split()
        self.assertAlmostEqual(float(x[1]),10)
    def test_unphysical_small_box_rejected(self):
        with self.assertRaises(ValueError):
            build("Nd",9,box=12)
if __name__=="__main__":unittest.main()
