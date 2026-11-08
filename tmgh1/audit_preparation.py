"""Check receptor heavy-atom preservation, pocket repairs and catalytic protonation."""
import json
from pathlib import Path
import numpy as np
from rdkit import Chem
out=Path(__file__).resolve().parents[1]/'results/tmgh1-round1'
def xyz(p,qt=False):
 return {(l[21],l[22:27].strip(),l[12:16].strip()):np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]) for l in p.read_text().splitlines() if l.startswith(('ATOM  ','HETATM')) and (l.split()[-1] not in ['H','HD'] if qt else l[76:78].strip() not in ['H','D'])}
reports=[]
for d in sorted((out/'receptors').iterdir()):
 if not d.is_dir():continue
 p=json.loads((d/'preparation.json').read_text());a=xyz(Path(p['source']));b=xyz(d/'receptor-h.pdb');c=xyz(d/'receptor.pdbqt',True);assert set(b)==set(c) and set(a)<=set(c);normalized=dict(c);swaps=[]
 for k,v in b.items():
  if np.linalg.norm(v-c[k])>.002:
   assert k[2] in ['O','OXT'];other=(*k[:2],'OXT' if k[2]=='O' else 'O');assert other in c and np.linalg.norm(v-c[other])<.002;normalized[k]=c[other]
   if k[2]=='O':swaps.append(':'.join(k[:2]))
 shift=max(float(np.linalg.norm(normalized[k]-v)) for k,v in a.items());assert shift<.002
 residues=set((k[0],k[1]) for k in c);assert len(residues)==p['residues'];near=[]
 for lig in ['NOJ','GI2']:
  l=Chem.SDMolSupplier(str(out/'reference'/f'{lig}.sdf'),removeHs=True)[0].GetConformer().GetPositions()
  for key,names in p['added_missing_heavy_atoms'].items():
   chain,idx=key.split(':');rr=np.array([v for k,v in a.items() if k[:2]==(chain,idx)]);distance=float(np.linalg.norm(rr[:,None]-l[None,:],axis=2).min())
   if distance<6:near.append({'ligand':lig,'residue':key,'modeled_atoms':names,'nearest_observed_distance_A':distance})
 h={i:[line[12:16].strip() for line in (d/'receptor-h.pdb').read_text().splitlines() if line.startswith('ATOM') and int(line[22:26])==i and line[76:78].strip()=='H'] for i in [166,351]};assert ('HE2' in h[166])==(p['state']=='E166H');assert 'HE2' not in h[351]
 report={'receptor':d.name,'all_original_heavy_atoms_preserved':True,'heavy_atoms':len(c),'residues':len(residues),'max_original_shift_A_after_O_OXT_mapping':shift,'terminal_O_OXT_label_swaps':swaps,'pocket_repairs_within_6A':near,'catalytic_hydrogens':h};(d/'audit.json').write_text(json.dumps(report,indent=2)+'\n');reports.append(report)
assert len(reports)==10; (out/'preparation-audit.json').write_text(json.dumps(reports,indent=2)+'\n');print('Audited',len(reports),'receptors; pocket repairs:',[(r['receptor'],r['pocket_repairs_within_6A']) for r in reports if r['pocket_repairs_within_6A']])
