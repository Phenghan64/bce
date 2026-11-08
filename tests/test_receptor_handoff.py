"""Scientific handoff gates using a real OpenFold3 7L39 prediction."""
import base64
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from io import BytesIO
import tarfile

import gemmi
from marine_agent.receptor_api import source_payload, resolve_prepared_receptor
from marine_agent import receptor_api, vina_remote
from marine_agent.receptor_tools import inspect_source, prepare_receptor
from marine_agent.vina_tools import validate_config

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'results/closed-loop-7l39-validated/openfold/prediction/pdb_7L39/seed_42/pdb_7L39_seed_42_sample_1_model.cif'


class ReceptorHandoffTests(unittest.TestCase):
    def setUp(self):
        self.raw = SOURCE.read_bytes()
        self.payload = {'source_run_id':'test-openfold-7l39', 'structure_base64':base64.b64encode(self.raw).decode(),
                        'options':{'confirmed':True,'chains':['A'],'ph':7.4,'padding':5,
                                   'pocket_residues':'A:99,A:102,A:118,A:121,A:153',
                                   'pocket_evidence':'PDB 7L39 L99A cavity residue mapping; workflow verification'}}
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)

    def modified(self, modifier):
        structure = gemmi.make_structure_from_block(gemmi.cif.read_string(self.raw.decode()).sole_block())
        modifier(structure)
        data = copy.deepcopy(self.payload)
        data['structure_base64'] = base64.b64encode(structure.make_mmcif_document().as_string().encode()).decode()
        return data

    def test_real_prediction_prepares_without_mutating_source(self):
        inspection = inspect_source(self.payload)
        self.assertEqual(inspection['chains'][0]['residues'],164)
        self.assertTrue(inspection['nonprotein'])
        result = prepare_receptor(self.payload,Path(self.temp.name)/'prepared')
        self.assertEqual(result['residues'],164)
        self.assertEqual(result['source_sha256'],hashlib.sha256(self.raw).hexdigest())
        self.assertLess(result['protonation']['max_existing_heavy_atom_shift_angstrom'],.002)
        self.assertGreater(result['protonation']['hydrogens_after'],0)
        self.assertEqual(SOURCE.read_bytes(),self.raw)
        self.assertGreaterEqual(result['pocket_quality']['mean_plddt'],70)

    def test_no_confirmation_cannot_prepare(self):
        self.payload['options']['confirmed']=False
        with self.assertRaisesRegex(ValueError,'确认'):
            prepare_receptor(self.payload,Path(self.temp.name))

    def test_low_confidence_pocket_is_rejected(self):
        def low(s):
            for atom in s[0]['A'][98]: atom.b_iso=20
        with self.assertRaisesRegex(ValueError,'质量门'):
            prepare_receptor(self.modified(low),Path(self.temp.name))

    def test_unknown_residue_is_not_silently_ignored(self):
        self.payload['options']['pocket_residues']='A:999'
        with self.assertRaisesRegex(ValueError,'不存在'):
            prepare_receptor(self.payload,Path(self.temp.name))

    def test_missing_side_chain_is_not_automatically_deleted(self):
        def remove(s):
            residue=s[0]['A'][98]
            for index in range(len(residue)-1,-1,-1):
                if residue[index].name=='CB': del residue[index]
        with self.assertRaisesRegex(ValueError,'Missing heavy atoms'):
            prepare_receptor(self.modified(remove),Path(self.temp.name))

    def test_chain_break_is_not_joined(self):
        def shift(s):
            for atom in s[0]['A'][98]: atom.pos.x += 20
        with self.assertRaisesRegex(ValueError,'主链断裂'):
            prepare_receptor(self.modified(shift),Path(self.temp.name))

    def test_nonstandard_residue_requires_manual_preparation(self):
        def modify(s): s[0]['A'][0].name='MSE'
        with self.assertRaisesRegex(ValueError,'非标准氨基酸'):
            prepare_receptor(self.modified(modify),Path(self.temp.name))

    def test_multichain_selection_retains_only_selected_chain(self):
        def duplicate(s):
            chain=s[0]['A'].clone();chain.name='B';s[0].add_chain(chain)
        payload=self.modified(duplicate)
        self.assertEqual(len(inspect_source(payload)['chains']),2)
        result=prepare_receptor(payload,Path(self.temp.name))
        self.assertEqual(result['residues'],164)
        self.assertEqual({item['source'].split(':')[0] for item in result['chain_mapping']},{'A'})

    def test_source_manifest_hash_mismatch_is_rejected(self):
        directory=Path(self.temp.name)/'results'; directory.mkdir()
        (directory/'bad.json').write_text(json.dumps({'structure':{'format':'mmcif','url':'model.cif','sha256':'0'*64},'provenance':{'model':'OpenFold3'}}))
        (directory/'model.cif').write_bytes(self.raw)
        with patch.object(receptor_api,'ROOT',Path(self.temp.name)):
            with self.assertRaisesRegex(ValueError,'校验和'):
                source_payload({'id':'test','status':'succeeded','result':{'agent':{}}},lambda _: 'results/bad.json')

    def test_prepared_receptor_is_resolved_from_ledger_and_hash_checked(self):
        directory=Path(self.temp.name)/'results/vina-receptors/prep-id'; directory.mkdir(parents=True)
        raw=b'ledger receptor'; (directory/'receptor.pdbqt').write_bytes(raw)
        metadata={'source_run_id':'source-id','source_sha256':'source-hash','receptor_sha256':hashlib.sha256(raw).hexdigest(),
                  'options':{},'pocket_quality':{},'chain_mapping':[],'box':{},'protonation':{},'tools':{},'excluded_nonprotein':[],'scientific_limit':'test'}
        class Store:
            def get(self, _): return {'id':'prep-id','mode':'vina_receptor','status':'succeeded','result':metadata}
        payload={'receptor':{'preparation_id':'prep-id','content_base64':'forged'}}
        with patch.object(receptor_api,'ROOT',Path(self.temp.name)):
            evidence=resolve_prepared_receptor(payload,Store())
            self.assertEqual(base64.b64decode(payload['receptor']['content_base64']),raw)
            self.assertEqual(evidence['source_run_id'],'source-id')
            (directory/'receptor.pdbqt').write_bytes(b'tampered')
            with self.assertRaisesRegex(ValueError,'校验失败'):
                resolve_prepared_receptor({'receptor':{'preparation_id':'prep-id'}},Store())

    def test_worker_bundle_cannot_write_outside_its_job(self):
        buffer=BytesIO()
        with tarfile.open(fileobj=buffer,mode='w:gz') as archive:
            item=tarfile.TarInfo('results/other-job/overwrite.txt'); item.size=4
            archive.addfile(item,BytesIO(b'evil'))
        with patch.object(vina_remote,'ROOT',Path(self.temp.name)):
            with self.assertRaisesRegex(RuntimeError,'unsafe artifact path'):
                vina_remote._safe_extract_bundle(base64.b64encode(buffer.getvalue()).decode(),'results/vina-receptors/prep-id')


if __name__=='__main__':
    unittest.main(verbosity=2)
