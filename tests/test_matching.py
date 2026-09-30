import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from line_redshift.core import Pair, Transition, read_lines, match_pair, C_KMS
ROOT=Path(__file__).resolve().parents[1]

class MatchingTests(unittest.TestCase):
    def setUp(self):
        self.lines=[Transition('CO21',230.538),Transition('CO32',345.7959899)]
        self.pair=Pair('synthetic',230.538/3,345.7959899/3,'LSRK','LSRK')
    def test_recovers_known_redshift(self):
        rows=match_pair(self.pair,self.lines,tolerance_kms=1)
        self.assertEqual(len(rows),1);self.assertAlmostEqual(rows[0].z,2,places=10)
        self.assertEqual((rows[0].line_a,rows[0].line_b),('CO21','CO32'))
    def test_swap_observations_swaps_labels(self):
        p=Pair('reverse',self.pair.frequency_b_ghz,self.pair.frequency_a_ghz,'LSRK','LSRK')
        row=match_pair(p,self.lines,1)[0]
        self.assertEqual((row.line_a,row.line_b),('CO32','CO21'));self.assertAlmostEqual(row.z,2)
    def test_boundary_is_pairwise_velocity_not_twice_tolerance(self):
        delta=100*3/(C_KMS-50)
        p=Pair('edge',230.538/3,345.7959899/(3+delta),'LSRK','LSRK')
        self.assertEqual(len(match_pair(p,self.lines,100)),1)
        self.assertEqual(len(match_pair(p,self.lines,99.9)),0)
    def test_mixed_frames_require_explicit_override(self):
        p=Pair('mixed',self.pair.frequency_a_ghz,self.pair.frequency_b_ghz,'LSRK','BARYCENT')
        with self.assertRaisesRegex(ValueError,'mixed frames'):match_pair(p,self.lines)
        self.assertEqual(match_pair(p,self.lines,allow_mixed=True)[0].frame_status,'unverified_mixed_frames')
    def test_unknown_frame_is_flagged(self):
        p=Pair('unknown',self.pair.frequency_a_ghz,self.pair.frequency_b_ghz,'','')
        self.assertEqual(match_pair(p,self.lines)[0].frame_status,'unverified_unknown_frame')
    def test_same_transition_alias_is_not_two_lines(self):
        self.assertEqual(match_pair(Pair('duplicate',100,100,'LSRK','LSRK'),[Transition('x',300),Transition('alias',300)]),[])
    def test_ambiguous_ladders_are_all_retained(self):
        rows=match_pair(Pair('ambiguous',100,150,'LSRK','LSRK'),[Transition('a',200),Transition('b',300),Transition('c',400),Transition('d',600)],1)
        self.assertEqual({round(r.z) for r in rows},{1,3})
    def test_invalid_inputs_and_z_bounds(self):
        for v in [float('nan'),float('inf'),0,-1]:
            with self.assertRaises(ValueError):match_pair(Pair('bad',v,100,'LSRK','LSRK'),self.lines)
        self.assertEqual(match_pair(self.pair,self.lines,z_max=1),[])
        with self.assertRaises(ValueError):match_pair(self.pair,self.lines,z_min=-1)
    def test_legacy_units_comments_and_duplicate_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'lines.freq';p.write_text('# MHz\n230538 CO21 # comment\n345795.9899 CO32\n')
            for a,b in zip(read_lines(p),self.lines):self.assertAlmostEqual(a.rest_ghz,b.rest_ghz)
            self.assertAlmostEqual(read_lines(p,unit='GHz')[0].rest_ghz,230538)
            p.write_text('230538 CO21\n230538 CO21\n')
            with self.assertRaises(ValueError):read_lines(p)
    def test_cli_outputs_all_matches_and_no_match_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'results'
            result=subprocess.run([sys.executable,str(ROOT/'all_spectra_plot.py'),'--lines',str(ROOT/'lines.freq'),'--pairs',str(ROOT/'examples/pairs.csv'),'--out',str(out),'--no-plots'],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            payload=json.loads((out/'hypotheses.json').read_text());self.assertEqual(len(payload['summary']),2)
            self.assertTrue(any(abs(r['z']-2)<1e-8 and r['line_a']=='CO21' and r['line_b']=='CO32' for r in payload['hypotheses']))
            self.assertEqual(payload['summary'][1]['hypotheses'],0)
            self.assertIn('No compatible assignment',(out/'report.md').read_text())
    def test_cli_validation_does_not_write_partial_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'pairs.csv';p.write_text('pair_id,frequency_a_ghz,frequency_b_ghz,frame_a,frame_b\nx,100,150,LSRK,BARYCENT\n')
            out=Path(tmp)/'out'
            result=subprocess.run([sys.executable,str(ROOT/'all_spectra_plot.py'),'--lines',str(ROOT/'lines.freq'),'--pairs',str(p),'--out',str(out),'--no-plots'],capture_output=True)
            self.assertNotEqual(result.returncode,0);self.assertFalse(out.exists())

if __name__=='__main__':unittest.main()
