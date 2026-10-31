"""Pure Python unit tests: published input parameters & correct unit conversion."""
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"simulations"/"openmm_hydration"))
sys.path.insert(0,str(ROOT/"simulations"/"cp2k_hydration"))
from run import METALS, ion_parameters
from generate_dimer import make_input
class TestHydration(unittest.TestCase):
    def test_limerz_tip3p(self):
        self.assertEqual(METALS["Nd"]["c4_kcal_A4"],213)
        self.assertEqual(METALS["Dy"]["c4_kcal_A4"],207)
        self.assertAlmostEqual(METALS["Nd"]["epsilon_kcal"],.15486311)
        self.assertAlmostEqual(METALS["Dy"]["epsilon_kcal"],.09900804)
    def test_unit_conversion(self):
        for el in ("Nd","Dy"):
            p=ion_parameters(el)
            self.assertAlmostEqual(p["c4_kJmol_nm4"],
                4.184e-4*METALS[el]["c4_kcal_A4"])
            self.assertAlmostEqual(p["sigma_nm"],
                METALS[el]["rmin_half_A"]*.2/2**(1/6))
    def test_input_charge_and_pseudos(self):
        for element in ("Nd","Dy"):
            s=make_input(element)
            self.assertIn("CHARGE 3",s)
            self.assertIn("DZVP-MOLOPT-GTH-q11",s)
            self.assertIn("PSOLVER MT",s)
            self.assertIn("RUN_TYPE ENERGY",s)
            self.assertIn("CUTOFF 250",s)
    def test_new_inputs_do_not_claim_optimized(self):
        for el in ("Nd","Dy"):
            self.assertNotIn("GEO_OPT",make_input(el))
if __name__=="__main__":
    unittest.main()
