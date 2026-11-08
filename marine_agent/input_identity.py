"""Resolve explicit input labels against UniProt and compare actual sequences."""
import hashlib
import json
import re
from .openfold_tools import ACCESSION_RE, get_uniprot_sequence


def inspect_input_identity(specs, sequence_index=None):
    entries = []
    seen_specs = set()
    for spec in specs:
        inputs = spec.get('input', {})
        signature = json.dumps(inputs, sort_keys=True)
        if signature in seen_specs:
            continue
        seen_specs.add(signature)
        fasta = inputs.get('query_fasta', '')
        for block in fasta.split('>'):
            lines = block.strip().splitlines()
            if len(lines) > 1:
                entries.append((lines[0].split()[0], ''.join(''.join(lines[1:]).split())))
        if inputs.get('sequence'):
            entries.append(('输入序列', inputs['sequence']))
    if sequence_index:
        entries = entries[sequence_index-1:sequence_index]
    results = []
    for label, sequence in entries[:10]:
        ids = [part for part in re.split(r'[_|\s]', label) if ACCESSION_RE.fullmatch(part)]
        row = {'label':label, 'length':len(sequence), 'sequence_sha256':hashlib.sha256(sequence.encode()).hexdigest(), 'identity_status':'unresolved'}
        if len(ids) == 1:
            reference = json.loads(get_uniprot_sequence.invoke({'accession':ids[0]}))
            if reference.get('status') == 'succeeded':
                parent = reference['sequence']
                row.update(reference_accession=reference['accession'], reference_metadata=reference['metadata'], reference_length=len(parent))
                if sequence == parent:
                    row['identity_status'] = 'exact_reference_sequence'
                elif len(sequence) == len(parent):
                    row['identity_status'] = 'differs_from_labeled_reference'
                    row['substitutions'] = [f'{a}{i}{b}' for i,(a,b) in enumerate(zip(parent,sequence),1) if a!=b]
                else:
                    row['identity_status'] = 'length_mismatch_with_labeled_reference'
        results.append(row)
    text = []
    for row in results:
        text.append(f"{row['label']}：保存输入为 {row['length']} 个残基。")
        if row['identity_status'] == 'unresolved':
            text.append('尚未获得可核验的蛋白身份，不能仅凭长度或序列标签认定功能。')
        else:
            meta = row['reference_metadata']
            text.append(f"标签中的 UniProt {row['reference_accession']} 对应 {meta.get('protein_name','名称未知')}；来源：{meta.get('source_url') or 'UniProt'}。")
            if row['identity_status'] == 'exact_reference_sequence':
                text.append('保存序列与该参考记录完全一致。')
            elif row.get('substitutions'):
                text.append('保存序列与参考的等长逐位差异：'+', '.join(row['substitutions'][:20])+'。这是输入差异核验，不是活性或功能验证。')
            else:
                text.append('保存序列与参考长度不一致，不能把参考名称直接当作输入身份。')
    return {'records':results,'final_text':'\n\n'.join(text) or '没有找到可核验的明确输入序列。','status':'succeeded' if results and all(x['identity_status']!='unresolved' for x in results) else 'partial'}
