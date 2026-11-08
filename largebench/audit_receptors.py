"""Verify all observed heavy coordinates survive receptor preparation, and audit pocket repairs."""
from pathlib import Path
import json,numpy as np
from rdkit import Chem
out=Path(__file__).resolve().parents[1]/'results/large-receptors'
def coordinates(p,qt=False):
 return {(l[21],l[22:27].strip(),l[12:16].strip()):np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]) for l in p.read_text().splitlines() if l.startswith(('ATOM  ','HETATM')) and (l.split()[-1] not in ['H','HD'] if qt else l[76:78].strip() not in ['H','D'])}
for name in ['albumin','mtor']:
 d=out/name/'receptor';a=coordinates(d/'observed.pdb');b=coordinates(d/'receptor-h.pdb');c=coordinates(d/'receptor.pdbqt',True);assert set(b)==set(c);assert set(a)<=set(c)
 swaps=[];normalized=dict(c)
 for k,v in b.items():
  if np.linalg.norm(v-c[k])>.002:
   assert k[2] in ['O','OXT'],'Non-terminal atom moved in Meeko'
   other=(*k[:2],'OXT' if k[2]=='O' else 'O');assert other in c and np.linalg.norm(v-c[other])<.002
   normalized[k]=c[other]
   if k[2]=='O':swaps.append(':'.join(k[:2]))
 shift=max(float(np.linalg.norm(normalized[k]-v)) for k,v in a.items());assert shift<.002
 audit={'observed_heavy_atoms':len(a),'prepared_heavy_atoms':len(b),'pdbqt_heavy_atoms':len(c),'max_observed_to_final_pdbqt_shift_A_after_terminal_O_OXT_name_mapping':shift,'terminal_carboxylate_O_OXT_label_swaps':swaps,'all_observed_atoms_preserved':True};(d/'coordinate-audit.json').write_text(json.dumps(audit,indent=2)+'\n')
 ligand=Chem.SDMolSupplier(str(d/'crystal-control.sdf'),removeHs=True)[0].GetConformer().GetPositions();prep=json.loads((d/'preparation.json').read_text());near=[]
 for key,names in prep['added_missing_heavy_atoms'].items():
  chain,idx=key.split(':');r=np.array([xyz for (ch,res,at),xyz in a.items() if ch==chain and res==idx]);dist=float(np.linalg.norm(r[:,None]-ligand[None,:],axis=2).min())
  if dist<8:near.append({'residue':key,'added_atoms':names,'nearest_observed_atom_to_ligand_A':dist})
 (d/'pocket-repairs.json').write_text(json.dumps(near,indent=2)+'\n');print(name,audit,'pocket_repairs',near)
