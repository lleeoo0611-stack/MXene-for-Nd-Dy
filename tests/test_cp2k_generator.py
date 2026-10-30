import sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"simulations"/"cp2k_hydration"))
from generate import directions,coords,make_input,LAN_THANIDES

class HydrationInputTests(unittest.TestCase):
    def test_water_count(self):
        for ion in LAN_THANIDES:
            for n in (8,9):
                c=coords(ion,n)
                self.assertEqual(len(c),1+3*n)
                self.assertEqual(sum(el=="O" for el,_ in c),n)
                self.assertEqual(sum(el=="H" for el,_ in c),2*n)
    def test_initial_ln_o_distances(self):
        for ion,r in LAN_THANIDES.items():
            for n in (8,9):
                c=coords(ion,n)
                center=c[0][1]
                for e,p in c:
                    if e=="O":
                        self.assertAlmostEqual(math.dist(center,p),r,places=6)
    def test_basis_charge_consistency(self):
        for ion in ("Nd","Dy"):
            p=make_input(ion,8)
            self.assertIn("CHARGE 3",p)
            self.assertIn("POTENTIAL GTH-PBE-q11",p)
            self.assertIn("BASIS_SET DZVP-MOLOPT-GTH-q11",p)
            self.assertIn("PSOLVER MT",p)
            self.assertIn("PERIODIC NONE",p)
    def test_invalid_input(self):
        with self.assertRaises(ValueError):make_input("La",8)
        with self.assertRaises(ValueError):make_input("Nd",7)
        with self.assertRaises(ValueError):make_input("Nd",8,"CELL_OPT")
if __name__=="__main__":unittest.main()
