#!/usr/bin/env python3
"""Guarded atomic publication of a marine-agent frontend release.

Run on lebb-node-13. --dry-run is read-only and never changes the homepage,
current symlink, target files or Nginx configuration.
"""
from __future__ import annotations
import argparse,datetime,hashlib,json,os,shutil,subprocess
from pathlib import Path

SITE=Path('/var/www/bce.dohoo.chat/current');TARGET=SITE/'marine-agent'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--release-dir',type=Path,required=True)
    parser.add_argument('--previous-receipt',type=Path,required=True)
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args();release=args.release_dir.resolve()
    manifest=json.loads((release/'files.json').read_text())
    if not manifest:raise RuntimeError('Empty release manifest')
    for name,expected in manifest.items():
        if not (name=='index.html' or name.startswith('assets/') or name in {'vina-workbench.html','vina-workbench.js','vina-workbench.css'}):
            raise RuntimeError('Unexpected frontend path: '+name)
        path=(release/'incoming'/name).resolve();path.relative_to((release/'incoming').resolve())
        if sha(path)!=expected:raise RuntimeError('Fingerprint mismatch: '+name)
    previous=json.loads(args.previous_receipt.read_text())
    if sha(TARGET/'index.html')!=previous['files']['index.html']:raise RuntimeError('Current frontend differs from reviewed receipt')
    home={name:sha(SITE/name) for name in ['index.html','portal.js']};resolved=str(SITE.resolve())
    preflight={'status':'preflight_passed','release_dir':str(release),'files':manifest,'homepage_sha256':home,'current':resolved,'dry_run':args.dry_run}
    print(json.dumps(preflight,ensure_ascii=False,indent=2))
    if args.dry_run:return 0
    subprocess.run(['sudo','-n','nginx','-t'],check=True)
    before=release/'before';before.mkdir(exist_ok=True)
    existed={name:(TARGET/name).exists() for name in manifest}
    for name in manifest:
        if existed[name]:
            target=before/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(TARGET/name,target)
    try:
        for name in sorted(manifest,key=lambda name:name=='index.html'):
            destination=TARGET/name;destination.parent.mkdir(parents=True,exist_ok=True)
            temporary=destination.with_name(destination.name+'.next-release')
            shutil.copyfile(release/'incoming'/name,temporary);temporary.chmod(0o644);os.replace(temporary,destination)
        if any(sha(TARGET/name)!=expected for name,expected in manifest.items()):raise RuntimeError('Published frontend hash mismatch')
        if {name:sha(SITE/name) for name in home}!=home or str(SITE.resolve())!=resolved:raise RuntimeError('Homepage or current link changed')
        subprocess.run(['sudo','-n','nginx','-t'],check=True)
    except Exception:
        for name in manifest:
            if existed[name]:shutil.copy2(before/name,TARGET/name)
            elif (TARGET/name).exists():(TARGET/name).unlink()
        raise
    receipt={'status':'published','release':release.name,'published_at':datetime.datetime.now().astimezone().isoformat(),
             'target':str(TARGET),'homepage_sha256':home,'files':manifest,'backup':str(before)}
    (release/'published.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(receipt,ensure_ascii=False,indent=2))
    return 0
if __name__=='__main__':raise SystemExit(main())
