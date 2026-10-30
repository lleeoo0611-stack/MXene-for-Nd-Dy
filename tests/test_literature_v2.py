import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from literature_pilot_v2 import (citrate_species,site_deprot,competitive_occupancy,
    normalized_downstream_flux,breakthrough_t50_dimensionless,permeability_ratio)
class LiteratureV2Tests(unittest.TestCase):
    def test_proton_mass_balance(self):
        q=citrate_species(np.linspace(3,8,20))
        np.testing.assert_allclose(q.sum(axis=0),1.)
    def test_citrate_deprotonation_increases(self):
        self.assertGreater(float(citrate_species(6.5)[3]),float(citrate_species(4.5)[3]))
    def test_surface_apparent_pka(self):
        self.assertAlmostEqual(float(site_deprot(5.5)),0.5)
    def test_competitive_site(self):
        n,d,v=competitive_occupancy(1e-4,1e-4,1e4,2e4,.5)
        self.assertAlmostEqual(n+d+v,1.)
        self.assertAlmostEqual(d/n,2.)
    def test_breakthrough_math(self):
        self.assertLess(float(normalized_downstream_flux(.01)),.01)
        self.assertGreater(float(normalized_downstream_flux(10)),.999)
        self.assertTrue(.1<breakthrough_t50_dimensionless()<.2)
    def test_independent_steady_perm(self):
        self.assertAlmostEqual(float(permeability_ratio(0,1)),1.)
        self.assertAlmostEqual(float(permeability_ratio(0,2)),2.)
if __name__=="__main__":unittest.main()
