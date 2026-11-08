"""Publish the already validated PubChem reference states without recalculation."""
import hashlib
import json
import shutil
from pathlib import Path
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
VERSION = 'pubchem20-20261030-v1'
SOURCE = ROOT/'results/ligand-library20-preflight-20261030'
DEST = ROOT/'results/ligand-libraries'/VERSION

def main():
    library = json.loads((SOURCE/'library.json').read_text())
    preview = json.loads((SOURCE/'catalog-preview.json').read_text())
    protocol = json.loads((SOURCE/'protocol.json').read_text())
    states = {x['id']: x for x in library['records']}
    DEST.mkdir(parents=True, exist_ok=True)
    for file in ('properties.json', 'structures-2d.sdf', 'retrieval.json'):
        (DEST/'source').mkdir(exist_ok=True)
        shutil.copy2(SOURCE/'source'/file, DEST/'source'/file)
    shutil.copy2(SOURCE/'protocol.json', DEST/'protocol.json')
    compounds = []
    for item in preview['compounds']:
        compound = {k: item[k] for k in ('id','cid','name_zh','name_en','pubchem_url','source_inchikey','chemistry_review_required')}
        compound['states'] = []
        for state_id in item['primary_state_ids']:
            state = dict(states[state_id])
            if state['role'] != 'primary_reference' or state['status'] != 'prepared':
                raise ValueError(state_id)
            stereo = [label for _, label in Chem.FindMolChiralCenters(Chem.MolFromSmiles(state['smiles']), includeUnassigned=True)]
            state['stereo_label'] = '/'.join(stereo) if stereo else ''
            state['display_name'] = item['name_zh'] + (f" ({state['stereo_label']})" if len(item['primary_state_ids']) > 1 else '')
            for field in ('pdbqt','sdf'):
                raw = (ROOT/state[field]).read_bytes()
                if hashlib.sha256(raw).hexdigest() != state[field+'_sha256']:
                    raise ValueError(f'Hash mismatch: {state_id}/{field}')
                relative = f"prepared/{state_id}/{'ligand.pdbqt' if field == 'pdbqt' else 'input.sdf'}"
                target = DEST/relative
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(raw)
                state[field] = relative
            compound['states'].append(state)
        compounds.append(compound)
    catalog = {'version':VERSION, 'title':'PubChem 参考配体库', 'compound_count':len(compounds),
               'state_count':sum(len(c['states']) for c in compounds), 'compounds':compounds,
               'state_policy':'羧酸采用羧酸根，酚基保持来源状态，布洛芬分别计算 R/S。属于明确的参考状态，不代表实验 pH 下的主要形态；多酚及姜黄素需结合实验条件复核。',
               'preparation':protocol['preparation'], 'versions':library['versions'],
               'sources':protocol['sources'], 'protocol_sha256':library['protocol_sha256']}
    raw = (json.dumps(catalog,ensure_ascii=False,indent=2)+'\n').encode()
    (DEST/'catalog.json').write_bytes(raw)
    print(json.dumps({'version':VERSION,'compounds':len(compounds),'states':catalog['state_count'],'catalog_sha256':hashlib.sha256(raw).hexdigest()}))

if __name__ == '__main__': main()
