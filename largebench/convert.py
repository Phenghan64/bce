import sys,json,time,subprocess
from pathlib import Path
root=Path(__file__).resolve().parents[1];d=Path(sys.argv[1]).resolve();m=json.loads((d/'metadata.json').read_text());t=time.perf_counter()
cmd=[str(root/'.env/bin/mk_prepare_receptor.py'),'--read_pdb',str(d/'receptor-h.pdb'),'-o',str(d/'receptor'),'-p','-j','-v','--box_size','30','30','30','--box_center',*map(str,m['box_center'])]
r=subprocess.run(cmd,capture_output=True,text=True);(d/'meeko.log').write_text(r.stdout+r.stderr);(d/'conversion.json').write_text(json.dumps({'seconds':time.perf_counter()-t,'returncode':r.returncode,'command':cmd},indent=2)+'\n');print(r.stdout[-1500:]+r.stderr[-1500:]);r.check_returncode()
text=(d/'receptor.pdbqt').read_text();lines=[x for x in text.splitlines() if x.startswith(('ATOM  ','HETATM'))];res=set((x[21],x[22:26],x[26]) for x in lines);assert len(res)==m['observed_residues'],(len(res),m['observed_residues'])
(d/'pdbqt-counts.json').write_text(json.dumps({'residues':len(res),'atoms':len(lines)},indent=2)+'\n')
