"""Convert the repaired receptor to Meeko PDBQT and record elapsed time."""
import json,subprocess,time
from pathlib import Path
root=Path(__file__).resolve().parents[1]
p=root/'results/abl-real-case/receptor';meta=json.loads((p/'metadata.json').read_text());start=time.perf_counter()
command=[str(root/'.env/bin/mk_prepare_receptor.py'),'--read_pdb',str(p/'receptor-h.pdb'),'-o',str(p/'receptor'),
         '-p','-j','-v','--box_size','30','30','30','--box_center',*map(str,meta['box_center'])]
r=subprocess.run(command,capture_output=True,text=True)
(p/'meeko.log').write_text(r.stdout+r.stderr)
(p/'meeko-timing.json').write_text(json.dumps({'seconds':time.perf_counter()-start,'returncode':r.returncode},indent=2)+'\n')
r.check_returncode()
