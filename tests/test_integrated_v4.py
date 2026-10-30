"""Regression tests for V4 hybrid conditioning, no DFT/MD claimed."""
import sys, math, unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from integrated_speciation_v4 import (config, baseline_table, solve_condition,
                                        calculate, assert_matching_csv, REPO)

class V4Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = config()
        cls.b = baseline_table()

    def test_basegrid_size(self):
        self.assertEqual(len(self.b), 72)

    def test_v4_counts(self):
        x,y,z = calculate(self.cfg,self.b)
        self.assertEqual((len(x),len(y),len(z)),(36,324,324))

    def test_zero_citrate_recovers_phreeqc_reference(self):
        for pH in self.cfg["pH"]:
            for s in self.cfg["NaNO3_background_mmolkgw"]:
                baseline=self.b[(pH,float(s),.1)]
                a=solve_condition(pH,s,0,baseline,self.cfg)
                for element in ("Nd","Dy"):
                    self.assertAlmostEqual(a[element+"_free_fraction"],baseline[element],places=11)
                    self.assertEqual(a[element+"_mixed_fraction"],0)

    def test_mass_balance(self):
        _, rows, _ = calculate(self.cfg,self.b)
        for r in rows:
            self.assertLess(abs(r["citrate_balance_error_M"]),1e-13)
            self.assertAlmostEqual(r["Nd_free_fraction"]+r["Nd_mixed_fraction"]+
                                   r["Nd_other_aq_fraction"],1.,places=11)
            self.assertAlmostEqual(r["Dy_free_fraction"]+r["Dy_mixed_fraction"]+
                                   r["Dy_other_aq_fraction"],1.,places=11)

    def test_citrate_beta_scenario_direction(self):
        b=self.b[(5.5,10.,.1)]
        a=solve_condition(5.5,10,.1,b,self.cfg,0,-.5)
        z=solve_condition(5.5,10,.1,b,self.cfg,0,+.5)
        self.assertGreater(z["Dy_mixed_fraction"],a["Dy_mixed_fraction"])

    def test_committed_numeric_archive(self):
        a,b,c=calculate(self.cfg,self.b)
        for filename,rows in [
            ("v4_citrate_baseline_36.csv",a),
            ("v4_beta_scenarios_324.csv",b),
            ("v4_conditional_transport_324.csv",c)]:
            assert_matching_csv(REPO/"results"/filename,rows)

    def test_transport_no_assumed_barrier(self):
        base,_,out=calculate(self.cfg,self.b)
        row=next(x for x in out if x["pH"]==5.5
                 and x["NaNO3_background_mmolkgw"]==10
                 and x["citrate_total_mM"]==.1
                 and x["assumed_entry_GDy_minus_GNd_kJmol"]==0
                 and x["assumed_mobile_DNd_over_DDy"]==1)
        self.assertAlmostEqual(row["free_ion_feed_Nd_over_Dy"],
            row["apparent_JNd_over_JDy_if_only_free_Ln_transports"])

    def test_status_is_not_experimental_selectivity(self):
        base, _, _ = calculate(self.cfg,self.b)
        self.assertEqual(base[0]["model"],"conditional_partial_v4")

if __name__=="__main__":
    unittest.main()
