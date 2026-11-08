#!/usr/bin/env python3
"""Generic, guarded backend publication for one verified release directory.

Run on lebb-node-1. --dry-run performs all read-only preflight checks without
writing files, restarting services or changing the database.
"""
from __future__ import annotations
import argparse,hashlib,json,shutil,sqlite3,subprocess,time,urllib.request
from pathlib import Path

APP=Path('/opt/marine-agent/app')
CURRENT_CONFIG=Path('/etc/systemd/system/marine-agent.service.d/99-zzzzzzzzzzzz-current-release.conf')

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def active_image():return subprocess.check_output(['docker','inspect','marine-agent-api','--format','{{.Config.Image}}'],text=True).strip()
def select_base_config(previous_image,config_root=Path('/etc/systemd/system/marine-agent.service.d')):
    candidates=sorted(config_root.glob('*.conf'))
    matches=[p for p in candidates if previous_image in p.read_text()]
    if not matches:raise RuntimeError('No systemd override contains the active image')
    return matches[-1]
def idle():
    with sqlite3.connect(APP/'results/marine-agent-api/state.sqlite3') as db:
        active=db.execute("SELECT count(*) FROM runs WHERE status IN ('queued','running','cancelling')").fetchone()[0]
    if active:raise RuntimeError(f'{active} active runs block release')

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--release-dir',type=Path,required=True)
    parser.add_argument('--image',required=True)
    parser.add_argument('--previous-image',required=True)
    parser.add_argument('--evidence',type=Path,required=True)
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args();release=args.release_dir.resolve()
    manifest=json.loads((release/'files.json').read_text())
    if not manifest:raise RuntimeError('Empty release manifest')
    for name,expected in manifest.items():
        path=(release/'incoming'/name).resolve();path.relative_to((release/'incoming').resolve())
        if sha(path)!=expected:raise RuntimeError('Fingerprint mismatch: '+name)
    evidence=json.loads(args.evidence.read_text())
    if evidence.get('status') not in {'succeeded','passed'}:raise RuntimeError('Verification evidence has not passed')
    if active_image()!=args.previous_image:raise RuntimeError('Active image differs from the reviewed base')
    if not shutil.which('docker'):raise RuntimeError('Docker unavailable')
    subprocess.check_output(['docker','image','inspect',args.image,'--format','{{.Id}}'])
    idle();base=select_base_config(args.previous_image)
    preflight={'status':'preflight_passed','release_dir':str(release),'image':args.image,'previous_image':args.previous_image,
               'files':manifest,'base_config':str(base),'current_config':str(CURRENT_CONFIG),'dry_run':args.dry_run}
    print(json.dumps(preflight,ensure_ascii=False,indent=2))
    if args.dry_run:return 0
    idle()
    with sqlite3.connect(APP/'results/marine-agent-api/state.sqlite3') as db,sqlite3.connect(release/'state-before.sqlite3') as backup:
        counts={table:db.execute('SELECT count(*) FROM '+table).fetchone()[0] for table in ['runs','conversations','research_projects','protein_candidates','experimental_measurements']}
        db.backup(backup)
    (release/'state-before.sqlite3').chmod(0o600)
    for name in manifest:
        source=APP/name
        if source.exists():
            target=release/'host-before'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
    previous_current=CURRENT_CONFIG.read_text() if CURRENT_CONFIG.exists() else None
    try:
        CURRENT_CONFIG.write_text(base.read_text().replace(args.previous_image,args.image))
        for name in manifest:
            target=APP/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(release/'incoming'/name,target)
        subprocess.run(['systemctl','daemon-reload'],check=True);subprocess.run(['systemctl','restart','marine-agent'],check=True)
        client=urllib.request.build_opener(urllib.request.ProxyHandler({}))
        for attempt in range(40):
            try:
                with client.open('http://100.64.0.1:8098/health',timeout=3) as response:
                    if json.load(response)['status']!='ok':raise RuntimeError('Health response is not ok')
                break
            except Exception:
                if attempt==39:raise
                time.sleep(.5)
        for name,expected in manifest.items():
            actual=hashlib.sha256(subprocess.check_output(['docker','exec','marine-agent-api','cat','/app/'+name])).hexdigest()
            if actual!=expected:raise RuntimeError('Running image hash mismatch: '+name)
        with sqlite3.connect(APP/'results/marine-agent-api/state.sqlite3') as db:
            for table,count in counts.items():
                if db.execute('SELECT count(*) FROM '+table).fetchone()[0]<count:raise RuntimeError('Database count regressed: '+table)
    except Exception:
        if previous_current is None:CURRENT_CONFIG.unlink(missing_ok=True)
        else:CURRENT_CONFIG.write_text(previous_current)
        for name in manifest:
            old=release/'host-before'/name
            if old.exists():shutil.copy2(old,APP/name)
            else:(APP/name).unlink(missing_ok=True)
        subprocess.run(['systemctl','daemon-reload'],check=True);subprocess.run(['systemctl','restart','marine-agent'],check=True)
        raise
    receipt={'status':'published','image':args.image,'previous_image':args.previous_image,'preserved_counts':counts,
             'database_backup':str(release/'state-before.sqlite3'),'health':'ok','files':manifest}
    (release/'published.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(receipt,ensure_ascii=False,indent=2))
    return 0
if __name__=='__main__':raise SystemExit(main())
