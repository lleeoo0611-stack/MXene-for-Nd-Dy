"""Stdlib-only reproducibility and mass-balance tests for V3."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from citrate_speciation_v3 import gamma,acid_fractions,solve,BETA_MIXED
class SpeciationTests(unittest.TestCase):
    def test_source_constants(self):
        self.assertEqual((BETA_MIXED["Nd"],BETA_MIXED["Dy"]),(13.05,13.28))
    def test_fractions(self):
        for pH in [4.5,5.5,6.5]:
            self.assertAlmostEqual(sum(acid_fractions(pH,.01)),1,places=12)
    def test_gamma(self):
        self.assertAlmostEqual(gamma(3,0),1)
        self.assertLess(gamma(3,.01),1)
    def test_zero_ligand(self):
        r=solve(5.5,.01,0)
        for el in ["Nd","Dy"]:self.assertEqual(r[el+"_mixed_fraction"],0)
    def test_metal_and_ligand_balance(self):
        for I in [.01,.1]:
            for pH in [4.5,5.5,6.5]:
                for L in [0,.0001,.0005]:
                    r=solve(pH,I,L)
                    self.assertLess(abs(r["citrate_balance_error_M"]),1e-12)
                    for el in ["Nd","Dy"]:
                        self.assertLess(abs(r[el+"_metal_balance_error_M"]),1e-12)
                        fsum=sum(r[el+"_"+x+"_fraction"] for x in ["free","hydroxo","mixed"])
                        self.assertAlmostEqual(fsum,1,places=12)
    def test_expected_order(self):
        r=solve(5.5,.01,.0001)
        self.assertGreater(r["Dy_mixed_fraction"],r["Nd_mixed_fraction"])
        self.assertAlmostEqual(r["Nd_mixed_fraction"],.1988,delta=.002)
        self.assertAlmostEqual(r["Dy_mixed_fraction"],.2957,delta=.002)
    def test_more_citrate(self):
        a=solve(5.5,.01,.0001);b=solve(5.5,.01,.0005)
        self.assertGreater(b["Nd_mixed_fraction"],a["Nd_mixed_fraction"])
        self.assertGreater(b["Dy_mixed_fraction"],a["Dy_mixed_fraction"])
    def test_log_beta_sensitivity(self):
        a=solve(5.5,.01,.0001);b=solve(5.5,.01,.0001,shifts={"Dy":-0.3})
        self.assertLess(b["Dy_mixed_fraction"],a["Dy_mixed_fraction"])
if __name__=="__main__":unittest.main()
