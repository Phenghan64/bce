"""Strict, dependency-free contract for InterProScan 6 execution artifacts."""
from __future__ import annotations
import hashlib
import json
import re

PROTOCOL = 'interproscan6-worker-v1'
ANALYSIS_ID = re.compile(r'^[0-9a-f]{64}$')
JOB_ID = re.compile(r'^(?:[0-9a-f]{32}|[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12})$')

def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False).encode()).hexdigest()

def active_library_versions(release):
    normalize = lambda value: re.sub(r'[^a-z0-9]', '', str(value).lower())
    active = {normalize(a['name']) for a in release.get('analyses', [])}
    return {k:v for k,v in release.get('library_versions', {}).items() if normalize(k) in active}

def validate_records(records):
    if not isinstance(records, list) or not 1 <= len(records) <= 200000:
        raise ValueError('expected 1–200000 protein records')
    seen, total = set(), 0
    for record in records:
        if not isinstance(record, dict): raise ValueError('invalid protein record')
        identifier, sequence = record.get('id'), record.get('sequence')
        if not isinstance(identifier, str) or not 1 <= len(identifier) <= 1000 or re.search(r'[\s>\x00-\x1f]', identifier) or identifier in seen:
            raise ValueError('invalid or duplicate protein identifier')
        if not isinstance(sequence, str) or not re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYBXZJUO]{1,1000000}', sequence):
            raise ValueError('invalid protein sequence')
        seen.add(identifier); total += len(sequence) + len(identifier)
    if total > 32 * 1024 * 1024: raise ValueError('FASTA exceeds 32 MiB transport limit')
    return records

def validate_scan(raw, records, release):
    """Require every exact sequence and every xref, including zero-hit proteins."""
    validate_records(records)
    if not isinstance(raw, dict) or raw.get('interpro-version') != release['interpro_version'] or raw.get('interproscan-version') != release['interproscan_version']:
        raise ValueError('InterProScan output version mismatch')
    expected = {}
    for record in records:
        md5 = hashlib.md5(record['sequence'].encode()).hexdigest().upper()
        entry = expected.setdefault(md5, {'sequence': record['sequence'], 'ids': set()})
        if entry['sequence'] != record['sequence']: raise ValueError('sequence MD5 collision')
        entry['ids'].add(record['id'])
    if not isinstance(raw.get('results'), list): raise ValueError('missing scan results')
    found = {}
    for item in raw['results']:
        md5 = str(item.get('md5', '')).upper()
        if md5 not in expected or md5 in found or item.get('sequence') != expected[md5]['sequence']:
            raise ValueError('unexpected, duplicate or altered sequence in scan output')
        xrefs = item.get('xref')
        if not isinstance(xrefs, list) or len(xrefs) != len(expected[md5]['ids']) or {x.get('id') for x in xrefs} != expected[md5]['ids']:
            raise ValueError('scan output identifiers do not match submitted proteins')
        if not isinstance(item.get('matches'), list): raise ValueError('missing matches array')
        normalize = lambda value: re.sub(r'[^a-z0-9]', '', str(value).lower())
        versions = {normalize(k): str(v) for k,v in release.get('library_versions', {}).items()}
        available = {normalize(a['name']) for a in release.get('analyses', [])}
        for match in item['matches']:
            library = (match.get('signature') or {}).get('signatureLibraryRelease') or {}
            name = normalize(library.get('library', ''))
            if name not in available or str(library.get('version')) != versions.get(name):
                raise ValueError('scan signature library is outside the pinned execution release')
            if match.get('source') and normalize(match['source']) != name:
                raise ValueError('unexpected analysis source in default full scan')
        found[md5] = item
    if set(found) != set(expected): raise ValueError('scan omitted one or more submitted sequences')
    return found
