"""Archived deterministic CP2K input snapshots contain *unoptimized* clusters."""
import sys, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"simulations"/"cp2k_hydration"))
from generate import coords, make_input
ROOT=Path(__file__).resolve().parents[1]
class SnapshotTests(unittest.TestCase):
    def test_snapshot_covers_four_expected_geometries(self):
        for element in ("Nd","Dy"):
            for n in (8,9):
                file=ROOT/"simulations"/"cp2k_hydration"/"snapshots"/f"{element.lower()}_w{n}.inp"
                text=file.read_text(encoding="utf8")
                self.assertEqual(text, make_input(element,n,run_type="ENERGY",cutoff=450))
                self.assertIn("CHARGE 3",text)
                self.assertIn("POTENTIAL GTH-PBE-q11",text)
    def test_xyz_atom_counts(self):
        for element in ("Nd","Dy"):
            for n in (8,9):
                file=ROOT/"simulations"/"cp2k_hydration"/"snapshots"/f"{element.lower()}_w{n}.xyz"
                lines=file.read_text().splitlines()
                self.assertEqual(int(lines[0]),1+3*n)
                self.assertEqual(len(lines)-2,1+3*n)
if __name__=="__main__":unittest.main()
