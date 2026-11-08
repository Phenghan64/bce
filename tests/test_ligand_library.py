"""Real frozen library integrity, identity and transport regression checks."""
import base64
import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from io import BytesIO

from marine_agent import ligand_library as lib, vina_remote
from marine_agent.vina_tools import validate_config

ROOT = Path(__file__).resolve().parents[1]

class LigandLibraryTests(unittest.TestCase):
    def request(self, ids=('L001',)):
        receptor = (ROOT/'results/1iep/receptor.pdbqt').read_bytes()
        return {'receptor':{'name':'receptor.pdbqt','content_base64':base64.b64encode(receptor).decode()},
                'ligands':[{'library_compound_id':i,'library_version':lib.VERSION} for i in ids],
                'box':{'center':[15.614,53.38,15.455],'size':[30,30,30]},'params':{}}

    def test_full_library_exactly_20_compounds_21_prepared_states(self):
        catalog = lib.public_catalog()
        self.assertEqual((catalog['compound_count'],catalog['state_count']),(20,21))
        self.assertNotIn('pdbqt',catalog['compounds'][0]['states'][0])
        request = self.request([c['id'] for c in catalog['compounds']])
        lib.resolve_library_ligands(request)
        config = validate_config(request)
        self.assertEqual(len(config['ligands']),21)
        for ligand in config['ligands']:
            metadata = config['provenance']['ligand_library'][ligand['id']]
            self.assertEqual(hashlib.sha256(ligand['bytes']).hexdigest(),metadata['pdbqt_sha256'])
            self.assertTrue(metadata['source_inchikey'])
            self.assertEqual(metadata['version'],lib.VERSION)

    def test_ibuprofen_expands_into_both_explicit_enantiomers(self):
        request = self.request(('L003',)); lib.resolve_library_ligands(request)
        self.assertEqual(len(request['ligands']),2)
        self.assertEqual({x['stereo_label'] for x in request['ligand_library'].values()},{'R','S'})
        self.assertNotEqual(*[x['smiles'] for x in request['ligand_library'].values()])

    def test_upload_mixed_with_library_and_stale_provenance_removed(self):
        request = self.request()
        raw = (ROOT/'results/1iep/ligand.pdbqt').read_bytes()
        uploaded = {'id':'uploaded','name':'custom.pdbqt','content_base64':base64.b64encode(raw).decode()}
        request['ligands'].append(uploaded)
        request['ligand_library']={'forged':'fake'}
        lib.resolve_library_ligands(request)
        self.assertEqual(request['ligands'][-1],uploaded)
        self.assertEqual(len(validate_config(request)['ligands']),2)
        self.assertNotIn('forged',request['ligand_library'])

    def test_unknown_version_or_compound_rejected(self):
        for key,value in [('library_version','../../secret'),('library_compound_id','L999')]:
            request = self.request(); request['ligands'][0][key]=value
            with self.assertRaises(ValueError): lib.resolve_library_ligands(request)

    def test_duplicate_selection_rejected(self):
        with self.assertRaisesRegex(ValueError,'重复'):
            lib.resolve_library_ligands(self.request(('L001','L001')))

    def test_reference_cannot_override_bytes_or_state(self):
        for key in ('content_base64','cid','id','state_id'):
            request=self.request(); request['ligands'][0][key]='forged'
            with self.assertRaisesRegex(ValueError,'自定义'):
                lib.resolve_library_ligands(request)

    def test_upload_cannot_claim_reserved_library_identity(self):
        request=self.request(); request['ligands']=[{'id':'library-L001-primary-1'}]
        with self.assertRaisesRegex(ValueError,'身份'): lib.resolve_library_ligands(request)

    def test_missing_corrupt_or_escaping_artifact_is_unavailable(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); (root/'a').write_bytes(b'corrupt')
            for relative,expected in [('missing',''),('a','wrong'),('../secret','')]:
                with self.assertRaises(lib.LibraryUnavailable): lib._verified_file(root,relative,expected)

    def test_catalog_corruption_is_not_silently_served(self):
        with patch.dict(lib.CATALOG_HASHES,{lib.VERSION:'wrong'}):
            with self.assertRaises(lib.LibraryUnavailable): lib.public_catalog()

    def test_state_expansion_respects_100_job_limit(self):
        request=self.request(('L003',))
        request['ligands'] += [{'id':f'upload-{i}'} for i in range(99)]
        with self.assertRaisesRegex(ValueError,'100'): lib.resolve_library_ligands(request)

    def test_mac_transport_preserves_server_provenance(self):
        request=self.request(('L003',)); lib.resolve_library_ligands(request)
        config=validate_config(request); captured={}
        def capture(req, timeout):
            captured.update(json.loads(req.data))
            return BytesIO(b'{"status":"succeeded","result":{"status":"succeeded"}}')
        with patch.object(vina_remote,'VINA_WORKER_URL','http://worker'), patch.object(vina_remote,'urlopen',capture):
            vina_remote.run_vina_worker('12345678-1234-1234-1234-123456789012',config)
        received=validate_config(captured['request'])
        self.assertEqual(received['provenance']['ligand_library'],config['provenance']['ligand_library'])

if __name__=='__main__': unittest.main()
