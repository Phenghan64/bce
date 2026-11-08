"""Verify completed-stage resume preserves outputs and creates no new Vina calls."""
import hashlib,json,subprocess,time
from pathlib import Path

root=Path(__file__).resolve().parents[1];out=root/'results/vina-throughput-1000'
assert not (out/'runner.lock').exists(), 'Wait for the active runner to finish'
paths=sorted((out/'jobs/screen').glob('*/docked.pdbqt'))
def snapshot():
    return {str(p):{'mtime_ns':p.stat().st_mtime_ns,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths}
before=snapshot();stage_before=(out/'screen.json').read_bytes();start=time.perf_counter()
command=[str(root/'run.sh'),'python',str(root/'throughput/run.py'),'--output-dir',str(out),'--stage','screen']
r=subprocess.run(command,capture_output=True,text=True);(out/'resume.log').write_text(r.stdout+r.stderr);r.check_returncode()
after=snapshot();check=json.loads((out/'screen-resume-check.json').read_text())
assert before==after and (out/'screen.json').read_bytes()==stage_before
assert check['cache_hits']==len(check['results'])
report={'wall_seconds_including_startup':time.perf_counter()-start,'cache_hits':check['cache_hits'],
    'result_count':len(check['results']),'qc_pass':check['completed'],'failed_preserved':check['failed'],
    'pdbqt_hashes_and_mtimes_unchanged':True,'original_stage_summary_unchanged':True,
    'new_vina_jobs':sum(not r['cache_hit'] for r in check['results']),
    'scope':'Completed-stage resume, not an interruption or power-loss simulation.'}
(out/'resume-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
