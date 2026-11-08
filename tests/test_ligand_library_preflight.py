"""Meaningful chemistry/cache gates for the exact public twenty-compound panel."""
import copy,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from rdkit import Chem
from ligand_library import preflight as p
from marine_agent.vina_tools import validate_config,VinaInputError

class LibraryPreflightTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library=json.loads((p.OUT/'library.json').read_text())
        cls.records=cls.library['records']
    def test_twenty_source_compounds_and_distinct_stereo_variants(self):
        primary=[r for r in self.records if r['role']=='primary_reference']
        self.assertEqual(len({r['cid'] for r in primary}),20)
        ibu=[r for r in primary if r['cid']==3672]
        self.assertEqual(len(ibu),2);self.assertEqual(len({r['smiles'] for r in ibu}),2)
        self.assertEqual(len({r['id'] for r in self.records}),len(self.records))
    def test_all_prepared_files_and_round_trip_identities(self):
        self.assertEqual(p.audit(),33)
    def test_unresolved_stereocenter_is_not_silently_embedded(self):
        with tempfile.TemporaryDirectory() as t:
            with self.assertRaisesRegex(ValueError,'stereochemistry'):
                p.prepare_one({'id':'unknown','smiles':'FC(Cl)Br'},Path(t))
    def test_disconnected_components_require_review(self):
        with tempfile.TemporaryDirectory() as t:
            with self.assertRaisesRegex(ValueError,'disconnected'):
                p.prepare_one({'id':'salt','smiles':'CC(=O)[O-].[Na+]'},Path(t))
    def test_stereo_inversion_is_detected(self):
        r=next(r for r in self.records if r['cid']==65064)
        mol=Chem.SDMolSupplier(str(p.ROOT/r['sdf']),removeHs=True)[0];expected=Chem.MolFromSmiles(r['smiles'])
        chiral=next(a for a in mol.GetAtoms() if a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED)
        chiral.InvertChirality()
        with self.assertRaisesRegex(ValueError,'stereochemistry'):p.exact_identity(mol,expected)
    def test_charge_change_is_detected(self):
        r=next(r for r in self.records if r['id']=='L002-primary-1');mol=Chem.SDMolSupplier(str(p.ROOT/r['sdf']),removeHs=True)[0]
        atom=next(a for a in mol.GetAtoms() if a.GetFormalCharge()==-1);atom.SetFormalCharge(0);atom.SetNumExplicitHs(1)
        with self.assertRaisesRegex(ValueError,'charge'):p.exact_identity(mol,Chem.MolFromSmiles(r['smiles']))
    def test_corrupted_cache_cannot_be_reused(self):
        with tempfile.TemporaryDirectory() as t:
            r=copy.deepcopy(self.records[0]);r['pdbqt']=str(Path(t)/'tampered.pdbqt');Path(r['pdbqt']).write_text('tampered')
            Path(t,'library.json').write_text(json.dumps({'records':[r]}))
            with patch.object(p,'OUT',Path(t)):
                with self.assertRaisesRegex(ValueError,'checksum'):p.audit()
    def test_all_and_subset_payloads_preserve_selected_ids(self):
        primary=[r for r in self.records if r['role']=='primary_reference']
        for selected in (primary,primary[:3]):
            config=validate_config(p.payload(selected))
            self.assertEqual([r['id'] for r in config['ligands']],[r['id'] for r in selected])
    def test_duplicate_selection_is_rejected(self):
        with self.assertRaises(VinaInputError):validate_config(p.payload([self.records[0],self.records[0]]))
    def test_missing_ligand_content_is_rejected(self):
        payload=p.payload(self.records[:1]);payload['ligands'][0]['content_base64']=''
        with self.assertRaises(VinaInputError):validate_config(payload)
    def test_resume_cannot_change_submitted_parameters(self):
        with tempfile.TemporaryDirectory() as t:
            folder=Path(t)/'api';folder.mkdir();(folder/'existing-submission.json').write_text(json.dumps({'signature':'wrong','run_id':'existing'}))
            with patch.object(p,'OUT',Path(t)):
                with self.assertRaisesRegex(ValueError,'differs'):p.submit('existing',{'changed':True})

if __name__=='__main__':unittest.main(verbosity=2)
