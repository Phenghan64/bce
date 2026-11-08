"""Freeze inputs and the retrospective, label-free screening protocol."""
import ast, collections, datetime, hashlib, json, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/marine-screen'
SRC = OUT / 'source'
AA = set('ACDEFGHIKLMNPQRSTVWY')

def save(name, obj):
    (OUT / name).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n')

def protein(x):
    pd = x['proteinDescription']
    name = pd.get('recommendedName', pd.get('submissionNames', [{}])[0]).get('fullName', {}).get('value', '')
    return dict(accession=x['primaryAccession'], name=name, species=x['organism']['scientificName'],
                taxon=x['organism']['taxonId'], sequence=x['sequence']['value'],
                fragment=pd.get('flag') == 'Fragment', reviewed=x['entryType'].startswith('UniProtKB reviewed'))

def main():
    raw = json.loads((SRC / 'marine-panel.json').read_text())['results']
    pattern = r'^(actin(?:,|$| alpha| beta)|beta.actin|alpha.actin|tropomyosin|troponin|parvalbumin|beta.enolase|alpha.enolase|enolase|fructose-bisphosphate aldolase|glyceraldehyde-3-phosphate dehydrogenase|myosin heavy chain, fast skeletal|myosin light chain|collagen alpha-[12]\(I\)|collagen, type I,)'
    grouped, excluded = {}, []
    for x in raw:
        p = protein(x)
        reason = ('unrelated_name' if not re.search(pattern, p['name'], re.I) else
                  'fragment' if p['fragment'] else 'noncanonical_sequence' if set(p['sequence']) - AA else None)
        if reason:
            excluded.append(dict(accession=p['accession'], reason=reason)); continue
        grouped.setdefault(p['sequence'], []).append(p)
    proteins = []
    for seq, entries in sorted(grouped.items(), key=lambda t: min(x['accession'] for x in t[1])):
        proteins.append(dict(id=min(x['accession'] for x in entries), sequence=seq, entries=entries))
    ids = ['B5DGQ7','Q5XJ10','O93484','B5DGM7','Q90339','P13104','P49055','P82159','P84335']
    benchmarks = [protein(json.loads((SRC / ('uniprot-' + i + '.json')).read_text())) for i in ids]
    # Read a literal rule dictionary only. Never execute downloaded source code.
    tree = ast.parse((SRC / 'cleavage-rules.py').read_text())
    rules = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
                 and any(isinstance(t, ast.Name) and t.id == 'expasy_rules' for t in n.targets))
    keys = ['trypsin','trypsin_exception','chymotrypsin high specificity','pepsin ph1.3','pepsin ph2.0','thermolysin','proteinase k']
    rules = {k: rules[k] for k in keys}
    protocol = dict(
        version='1.0', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        design='Retrospective protocol frozen after reading published experimental outcomes, before computing rankings; not blinded.',
        target='ACE inhibition; no claim of antihypertensive efficacy in humans',
        library='Complete UniProt stream query; named muscle/metabolic/type-I-collagen families; remove fragments and noncanonical sequences; exact sequence dedup, retain all accession/species provenance.',
        protein_filter_regex=pattern, min_length=2, max_length=20,
        cleavage='Pyteomics ExPASy regex literal snapshot; trypsin exception enabled; cuts are regex match ends.',
        rules=rules,
        schemes={'trypsin':['trypsin'], 'chymotrypsin':['chymotrypsin high specificity'],
                 'pepsin_pH1.3':['pepsin ph1.3'], 'pepsin_pH2':['pepsin ph2.0'],
                 'thermolysin':['thermolysin'], 'proteinase_K':['proteinase k'],
                 'trypsin+chymotrypsin':['trypsin','chymotrypsin high specificity'],
                 'pepsin_pH2+trypsin':['pepsin ph2.0','trypsin'],
                 'thermolysin+trypsin':['thermolysin','trypsin']},
        combinations='Union of cleavage sites on original sequence; idealized maximal cleavage scenario. Not a simultaneous compatible-pH recipe or model of temporal order.',
        missed_cleavage='MC0 primary complete cleavage; MC1 sensitivity is union of zero and one internal missed site. Do not sum alternative fragments as mass yield.',
        score='2*(C-terminal residue in FWYP) + 1*(penultimate in AVILMFYW) + 1*(length <=5); range 0..4; fixed heuristic, no trained model, not potency or probability.',
        baselines=['shorter sequence (-length)','GRAVY hydrophobicity','uniform random selection'],
        experimental_reference='ACEip_DB: canonical 2..20aa, exact single positive uM or mM value, DOI present; exclude ranges, censored, malformed, zero, missing units; heterogeneous literature reference, not all primary measurements reverified.',
        validation='No fitting/tuning. Spearman against -log10(IC50/uM). Top ceil(10% N), strong <=10uM, fractional boundary ties; EF vs random within measured reference only. Bootstrap 1000, permutation 2000, seed 20261022.',
        sardine='Six numeric single-peptide IC50s in Foods2023 12:2216, convert mg/mL to uM = mg/mL*1e6/MW; six are positive selected measurements, no inactive labels; 37 LCMS peptides not full peptidome.',
        process_benchmark='Foods2022 11:1785 eight homolog accessions from TableS1; P13104 primary tropomyosin, replace by P84335 sensitivity. Equal-molecule sum of high-score (>=3) fragments per1000 total residues primary; equal-protein normalized sensitivity. Compare pepsin, thermolysin, proteinaseK only (3/6 tested enzymes); pepsin pH2.0 rule vs experimental pH2 mismatch documented.',
        process_ground_truth='Pepsin strongest ACE inhibition at0.1mg/mL (86.14% body); other five not significantly different. Winner retrieval only, no invented exact rank. Papain/bromelain/Alcalase product rules unavailable: not simulated.',
        zenodo='Auxiliary only: author-deposited 2026 dataset; some curves explicitly back-calculated, raw instrument provenance and process specifications unresolved; not independent primary validation of this pipeline.',
        sensitivity='MC0 vs MC1; equal molecule vs normalized per-protein benchmark; longest protein per species+name reduces isoform sampling bias in marine process summary.',
        limitations=['No abundances, cleavage kinetics, PTMs, structure accessibility, denaturation or peptide recovery model.',
                     'Database annotation match is not a new experimental validation. Untested is not inactive.',
                     'No toxicity/allergenicity or oral bioavailability clearance; no production-yield prediction.',
                     'No OpenFold3/Vina jobs; sequence screening does not require 3D docking.'], seed=20261022)
    if (OUT / 'protocol.json').exists():
        raise RuntimeError('Protocol already frozen; changes require explicit versioned amendment.')
    save('proteins.json', proteins); save('benchmark-proteins.json', benchmarks)
    save('protein-exclusions.json', excluded); save('protocol.json', protocol)
    save('input-manifest.json', [dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,
        sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(SRC.rglob('*')) if p.is_file()])
    save('input-audit.json', dict(queried=len(raw),accepted_entries=sum(len(p['entries']) for p in proteins),
         unique_sequences=len(proteins),total_residues=sum(len(p['sequence']) for p in proteins),
         species_entry_counts=dict(collections.Counter(e['species'] for p in proteins for e in p['entries'])),
         reviewed_entries=sum(e['reviewed'] for p in proteins for e in p['entries']),
         exclusion_reasons=dict(collections.Counter(x['reason'] for x in excluded))))
    print(json.dumps(json.loads((OUT/'input-audit.json').read_text()),indent=2))

if __name__ == '__main__': main()
