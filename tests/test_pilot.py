import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from phreeqc_pilot import input_text,PHS,SALTS,METALS
from energy_sensitivity import rate_ratio,R
from toy_transport import mfpt

class PilotTests(unittest.TestCase):
    def test_grid_count_and_nitrate(self):
        self.assertEqual(len(PHS)*len(SALTS)*len(METALS),72)
        src=input_text(6.5,10,0.1)
        self.assertIn('N(5) 10.600000',src)
        self.assertIn('Nd 0.100000',src)
    def test_activation(self):
        self.assertAlmostEqual(float(rate_ratio(0)),1)
        self.assertAlmostEqual(float(rate_ratio(R*298.15*np.log(2))),2)
    def test_trap_slows_exit(self):
        a=np.zeros(11); b=a.copy();b[4]=-5
        self.assertGreater(mfpt(b),mfpt(a))
    def test_deep_equal_traps_yield_ratio_one(self):
        a=np.zeros(11);a[4]=-6
        self.assertAlmostEqual(mfpt(a)/mfpt(a),1)

if __name__=='__main__':unittest.main()
