"""Score observed chain-A substitutions and cross-check the upstream example."""
import csv, json, time, sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import torch
from omegaconf import OmegaConf
base = Path(__file__).resolve().parent
repo = base/'ThermoMPNN-main'
sys.path.insert(0,str(repo))
from transfer_model import TransferModel, ALPHABET
from protein_mpnn_utils import alt_parse_PDB
out=base.parents[1]/'results/enzyme-stability'
cfg=OmegaConf.create({'platform':{'thermompnn_dir':str(repo)},'model':{'hidden_dims':[64,32],'subtract_mut':True,'num_final_layers':2,'freeze_weights':True,'load_pretrained':True,'lightattn':True}})
torch.set_num_threads(4)
model=TransferModel(cfg).eval()
state=torch.load(repo/'models/thermoMPNN_default.pt',map_location='cpu')['state_dict']
model.load_state_dict({k.removeprefix('model.'):v for k,v in state.items() if k.startswith('model.')},strict=True)
def score(pdb, muts):
    with torch.inference_mode():
        predictions,_=model(pdb,muts)
    return np.array([p['ddG'].item() for p in predictions])
def mutation(i,w,m):return SimpleNamespace(position=i,wildtype=w,mutation=m)
with (repo/'examples/ThermoMPNN_inference_2OCJ.csv').open() as file:
    reference=list(csv.DictReader(file))[::17]
example=alt_parse_PDB(str(repo/'examples/2OCJ.pdb'),'A')
actual=score(example,[mutation(int(x['position']),x['wildtype'],x['mutation']) for x in reference])
err=np.max(abs(actual-np.array([float(x['ddG_pred']) for x in reference])))
start=time.monotonic()
p=alt_parse_PDB(str(out/'1BLI.pdb'),'A')
coords=p[0]['coords_chain_A']
seq=p[0]['seq']
positions=[i for i,w in enumerate(seq) if w in ALPHABET[:-1] and all(np.isfinite(v[i]).all() for v in coords.values())]
muts=[mutation(i,seq[i],a) for i in positions for a in ALPHABET[:-1]]
y=score(p,muts)
assert np.isfinite(y).all() and all(v==0 for m,v in zip(muts,y) if m.wildtype==m.mutation)
with (out/'thermo.csv').open('w') as f:
    writer=csv.writer(f); writer.writerow(['chain','pdb_position','parent','substitution','ddG_pred_kcal_mol'])
    for m,v in zip(muts,y):
        if m.wildtype!=m.mutation:writer.writerow(['A',p[0]['resn_list'][m.position],m.wildtype,m.mutation,v])
stats=dict(residues=len(positions),substitutions=len(positions)*19,elapsed_s=time.monotonic()-start,reference_rows=len(reference),reference_max_absolute_error=float(err))
print(stats);(out/'thermo-check.json').write_text(json.dumps(stats,indent=2))
import pandas as pd
frame=pd.read_csv(repo/'data_all/testing/fireprot_test.csv').query("pdb_id_corrected == '1E0W'").copy()
xyl=alt_parse_PDB(str(out/'1E0W.pdb'),'A')
muts=[mutation(int(x.pdb_position),x.wild_type,x.mutation) for x in frame.itertuples()]
assert all(xyl[0]['seq'][m.position]==m.wildtype for m in muts)
frame['predicted']=score(xyl,muts)
frame['pdb_author_number']=[xyl[0]['resn_list'][m.position] for m in muts]
frame[['pdb_author_number','wild_type','mutation','ddG','predicted','publication_doi']].to_csv(out/'xylanase-validation.csv',index=False)
print('Spearman',frame.ddG.rank().corr(frame.predicted.rank()),'sign correct',sum(np.sign(frame.ddG)==np.sign(frame.predicted)), 'pred_negative',sum(frame.predicted<0),'actual_negative',sum(frame.ddG<0))
