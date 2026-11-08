"""Scientific invariants, not merely mirrors of implementation."""
import unittest
import numpy as np
from run import cuts, fragments, properties, top_weights, score

class CoreTests(unittest.TestCase):
    def test_trypsin_proline_and_exceptions(self):
        self.assertEqual(cuts('AKPAAKAA',['trypsin']),[0,6,8])
        self.assertEqual(cuts('AWKPAA',['trypsin']),[0,3,6])
        self.assertEqual(cuts('ACKDAA',['trypsin']),[0,6])
    def test_thermolysin_n_side_and_block(self):
        self.assertEqual(cuts('GAGLQ',['thermolysin']),[0,1,3,5])
        self.assertEqual(cuts('DALP',['thermolysin']),[0,4])
    def test_conservation_and_coordinates(self):
        seq='MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE'
        products=list(fragments(seq,['trypsin'],0,False))
        self.assertEqual(''.join(p[0] for p in products),seq)
        for pep,a,b,mc in products:self.assertEqual(pep,seq[a:b])
    def test_combination_is_order_independent(self):
        seq='AKWFGLPLKRIYAA'
        self.assertEqual(cuts(seq,['trypsin','thermolysin']),cuts(seq,['thermolysin','trypsin']))
    def test_one_missed_is_adjacent_join(self):
        self.assertEqual({p[0] for p in fragments('AAKGGKTT',['trypsin'],1)}, {'AAK','GGK','TT','AAKGGK','GGKTT'})
    def test_mass_and_ties(self):
        self.assertAlmostEqual(properties('GG')['mw_Da'],132.1179,places=3)
        # Independent elemental mass for unmodified FIGR: C23H37N7O5.
        elemental_mw=23*12.0107+37*1.00794+7*14.0067+5*15.9994
        self.assertAlmostEqual(properties('FIGR')['mw_Da'],elemental_mw,places=3)
        self.assertAlmostEqual(.60*1e6/properties('FIGR')['mw_Da'],1220.545,places=2)
        self.assertTrue(np.allclose(top_weights([4,4,2],1),[.5,.5,0]))
        self.assertEqual(score('KLF'),4)

if __name__=='__main__':unittest.main()
