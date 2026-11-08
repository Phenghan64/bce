"""Append-only, owner-scoped candidate and experimental evidence ledger."""
from __future__ import annotations
import csv, hashlib, io, json, re, uuid
from datetime import datetime, timezone
from fastapi import HTTPException
from fastapi.responses import Response, JSONResponse
from pydantic import BaseModel, ConfigDict, Field, field_validator
from .experimental_data import COLUMNS, read_table, normalize_row, row_hash, comparisons


def stamp():return datetime.now(timezone.utc).isoformat()
def dump(value):return json.dumps(value,ensure_ascii=False,sort_keys=True)

class CandidateInput(BaseModel):
    model_config=ConfigDict(extra='forbid',str_strip_whitespace=True)
    key:str=Field(min_length=1,max_length=80,pattern=r'^[A-Za-z0-9_.-]+$')
    name:str=Field(min_length=1,max_length=160)
    sequence:str=Field(min_length=1,max_length=100000)
    source:str=Field(min_length=1,max_length=2000)
    construct_note:str=Field(default='',max_length=2000,alias='construct')
    @field_validator('sequence')
    @classmethod
    def sequence_value(cls,value):
        value=''.join(value.split()).upper()
        if not value or not re.fullmatch('[ACDEFGHIKLMNPQRSTVWYBXZJUO]+',value):raise ValueError('请输入单条蛋白序列，不包含 FASTA 标题或终止符')
        return value
class CandidateBatch(BaseModel):
    candidates:list[CandidateInput]=Field(min_length=1,max_length=1000)
class ImportInput(BaseModel):
    filename:str=Field(min_length=1,max_length=200)
    content_base64:str=Field(min_length=1,max_length=5400000)
class EvidenceLink(BaseModel):
    run_id:str=Field(min_length=1,max_length=64)
class CommitInput(BaseModel):
    token:str=Field(min_length=1,max_length=64)
class CandidateFromRun(BaseModel):
    model_config=ConfigDict(extra='forbid',str_strip_whitespace=True)
    run_id:str=Field(min_length=1,max_length=64)
    key:str=Field(min_length=1,max_length=80,pattern=r'^[A-Za-z0-9_.-]+$')
    name:str=Field(min_length=1,max_length=160)
    source:str=Field(min_length=1,max_length=2000)
    construct_note:str=Field(default='',max_length=2000,alias='construct')

class CandidateStore:
    def initialize_candidates(self):
        with self._connect() as db:
            db.executescript('''
            CREATE TABLE IF NOT EXISTS candidate_revisions(project_id TEXT PRIMARY KEY, revision INTEGER NOT NULL);
            CREATE TABLE IF NOT EXISTS protein_candidates(id TEXT PRIMARY KEY,project_id TEXT NOT NULL,candidate_key TEXT NOT NULL,
                revision INTEGER NOT NULL,payload TEXT NOT NULL,UNIQUE(project_id,candidate_key));
            CREATE TABLE IF NOT EXISTS candidate_evidence(project_id TEXT NOT NULL,candidate_id TEXT NOT NULL,run_id TEXT NOT NULL,
                revision INTEGER NOT NULL,payload TEXT NOT NULL,PRIMARY KEY(candidate_id,run_id));
            CREATE TABLE IF NOT EXISTS experiment_imports(id TEXT PRIMARY KEY,project_id TEXT NOT NULL,file_sha256 TEXT NOT NULL,
                revision INTEGER NOT NULL,filename TEXT NOT NULL,raw BLOB NOT NULL,created_at TEXT NOT NULL,UNIQUE(project_id,file_sha256));
            CREATE TABLE IF NOT EXISTS experimental_measurements(project_id TEXT NOT NULL,import_id TEXT NOT NULL,revision INTEGER NOT NULL,
                row_sha256 TEXT NOT NULL,payload TEXT NOT NULL,PRIMARY KEY(project_id,row_sha256));
            CREATE TABLE IF NOT EXISTS experiment_previews(id TEXT PRIMARY KEY,project_id TEXT NOT NULL,revision INTEGER NOT NULL,
                filename TEXT NOT NULL,file_sha256 TEXT NOT NULL,raw BLOB NOT NULL,payload TEXT NOT NULL,created_at TEXT NOT NULL,committed_id TEXT);
            CREATE INDEX IF NOT EXISTS idx_candidate_project ON protein_candidates(project_id,revision);
            CREATE INDEX IF NOT EXISTS idx_measurement_project ON experimental_measurements(project_id,revision);
            ''')
    def _ledger_revision(self,db,pid):
        row=db.execute('SELECT revision FROM candidate_revisions WHERE project_id=?',(pid,)).fetchone()
        return row[0] if row else 0
    def _next_ledger_revision(self,db,pid):
        db.execute('INSERT INTO candidate_revisions VALUES (?,1) ON CONFLICT(project_id) DO UPDATE SET revision=revision+1',(pid,))
        return self._ledger_revision(db,pid)
    def candidate_snapshot(self,db,pid):
        return {'revision':self._ledger_revision(db,pid),
                'candidate_count':db.execute('SELECT count(*) FROM protein_candidates WHERE project_id=?',(pid,)).fetchone()[0],
                'measurement_count':db.execute('SELECT count(*) FROM experimental_measurements WHERE project_id=?',(pid,)).fetchone()[0],
                'read_tool':'inspect_candidate_evidence','provenance':'user_registered_candidates_and_imported_measurements'}
    def register_candidates(self,pid,batch):
        with self._lock,self._connect() as db:
            db.execute('BEGIN IMMEDIATE');self._owned_project(db,pid)
            seen=set();pending=[];existing=[]
            for item in batch.candidates:
                if item.key in seen:raise HTTPException(422,'批次内候选编号重复：'+item.key)
                seen.add(item.key);data=item.model_dump(by_alias=True);data['sequence_sha256']=hashlib.sha256(item.sequence.encode()).hexdigest()
                old=db.execute('SELECT payload FROM protein_candidates WHERE project_id=? AND candidate_key=?',(pid,item.key)).fetchone()
                if old:
                    saved=json.loads(old[0])
                    if any(saved[k]!=data[k] for k in data):raise HTTPException(409,'编号已存在且内容不同，请使用新版本编号：'+item.key)
                    existing.append(saved);continue
                pending.append({**data,'id':str(uuid.uuid4()),'created_at':stamp(),'length':len(item.sequence),
                                'ambiguous_sequence':bool(re.search('[BXZJUO]',item.sequence)),'provenance':'user_registered_sequence'})
            if pending:
                rev=self._next_ledger_revision(db,pid)
                for data in pending:db.execute('INSERT INTO protein_candidates VALUES (?,?,?,?,?)',(data['id'],pid,data['key'],rev,dump(data)))
            return {'candidates':existing+pending,'created':len(pending),'revision':self._ledger_revision(db,pid)}
    def candidate_ledger(self,pid,*,revision=None,offset=0,limit=50,candidate_key=None):
        with self._connect() as db:
            self._owned_project(db,pid)
            revision=self._ledger_revision(db,pid) if revision is None else revision
            query='SELECT payload FROM protein_candidates WHERE project_id=? AND revision<=?'
            args=[pid,revision]
            if candidate_key:query+=' AND candidate_key=?';args.append(candidate_key)
            candidates=[json.loads(r[0]) for r in db.execute(query+' ORDER BY candidate_key',args)]
            total=len(candidates);all_candidates=candidates;candidates=candidates[max(0,offset):max(0,offset)+min(100,max(1,limit))]
            ids={c['id'] for c in candidates}
            measurements=[{**json.loads(r['payload']),'record_id':r['row_sha256']} for r in db.execute('SELECT row_sha256,payload FROM experimental_measurements WHERE project_id=? AND revision<=? ORDER BY row_sha256',(pid,revision))]
            # Comparison groups are computed over the full snapshot, never an accidentally paginated subset.
            groups=comparisons(measurements)
            links=[json.loads(r[0]) for r in db.execute('SELECT payload FROM candidate_evidence WHERE project_id=? AND revision<=?',(pid,revision))]
            for c in candidates:self._decorate_candidate(c,links,measurements)
            queue={}
            for c in all_candidates:
                self._decorate_candidate(c,links,measurements)
                queue[c['review_state']]=queue.get(c['review_state'],0)+1
            review_queue={'scope':'full_snapshot','total':len(all_candidates),'counts':queue}
            return {'revision':revision,'total':total,'offset':offset,'candidates':candidates,
                    'review_queue':review_queue,
                    'measurements':[m for m in measurements if m['candidate_id'] in ids],
                    'comparisons':groups if not candidate_key else [g for g in groups if any(v['candidate_key']==candidate_key for v in g['values'])],
                    'imports':[dict(r) for r in db.execute('SELECT id,filename,file_sha256,created_at FROM experiment_imports WHERE project_id=? AND revision<=?',(pid,revision))]}
    @staticmethod
    def _decorate_candidate(candidate,links,measurements):
        candidate['evidence']=[link for link in links if link['candidate_id']==candidate['id']]
        candidate['evidence'].sort(key=lambda item:item.get('created_at',''),reverse=True)
        reviews=[link['docking_review'] for link in candidate['evidence'] if link.get('docking_review')]
        candidate['docking_review_count']=len(reviews)
        candidate['docking_review_status']=reviews[0].get('review_status') if reviews else None
        candidate['measurement_count']=sum(row['candidate_id']==candidate['id'] for row in measurements)
        kinds=set()
        for link in candidate['evidence']:
            workflow=link.get('workflow') or ''
            if workflow.startswith('openfold'):kinds.add('structure')
            if workflow=='vina':kinds.add('docking')
            if workflow=='protein_screen':kinds.add('homology')
            if workflow=='protein_domains':kinds.add('domains')
            for operation in link.get('operations',[]):
                tool=operation.get('tool') or ''
                if tool.startswith('run_openfold'):kinds.add('structure')
                elif tool=='annotate_protein_domains':kinds.add('domains')
                elif tool=='run_protein_homology_screen':kinds.add('homology')
        candidate['evidence_kinds']=sorted(kinds)
        status=candidate['docking_review_status']
        missing=[]
        if not candidate['evidence']:missing.append('尚未关联计算证据')
        if 'domains' not in kinds and 'homology' not in kinds:missing.append('可补充结构域或同源证据')
        if 'structure' not in kinds:missing.append('尚无关联结构证据')
        if not reviews:missing.append('尚无对接复核')
        if not candidate['measurement_count']:missing.append('尚无实验记录')
        if status=='divergent_review_required':missing.append('存在跨种子姿势分歧')
        candidate['missing_evidence']=missing
        if status=='divergent_review_required':state='needs_pose_review'
        elif status=='analysis_unavailable':state='needs_analysis_review'
        elif status=='single_seed_review_required':state='needs_multiseed_review'
        elif status=='consistent_review_required':state='experimental_records_present' if candidate['measurement_count'] else 'ready_for_experiment'
        elif not candidate['evidence']:state='needs_computation'
        else:state='needs_evidence'
        candidate['review_state']=state
        candidate['recommended_next']='查看中心姿势、接触和实验条件；计算证据不等于实验验证。' if reviews else '先关联可追溯的计算证据，再进入实验复核。'
        return candidate

    def _run_sequence_records(self,run,depth=0):
        """Return exact saved sequences, not only their hashes."""
        from pathlib import Path
        specs=[run.get('execution_spec') or {}]+[op.get('spec') or {} for op in run.get('operations',[])]
        records=[]
        def add(sequence):
            sequence=''.join(str(sequence or '').split()).upper()
            if sequence:records.append({'sequence':sequence,'sequence_sha256':hashlib.sha256(sequence.encode()).hexdigest()})
        for spec in specs:
            value=spec.get('input') or {}
            if value.get('sequence'):add(value['sequence'])
            try:
                if value.get('query_fasta'):
                    from .protein_screening import parse_fasta
                    for row in parse_fasta(value['query_fasta']):add(row['sequence'])
                elif spec.get('workflow')=='protein_domains' and value.get('analysis_id'):
                    from .domain_evidence import resolve_inputs
                    _,rows,_=resolve_inputs(value)
                    for row in rows:add(row['sequence'])
            except (ValueError,OSError,HTTPException):pass
        if records:return records,'server_execution_sequence_sha256'
        from .api import _find_viewer_manifest
        from .result_contract import trusted_structure_manifest, verified_structure
        root=Path(__file__).resolve().parents[1]
        agent=(run.get('result') or {}).get('agent')
        manifest=trusted_structure_manifest(agent,_find_viewer_manifest)
        if manifest and verified_structure(run,root,_find_viewer_manifest):
            try:
                path=(root/manifest).resolve().parent.parent/'input.json'
                path.resolve().relative_to((root/'results').resolve())
                data=json.loads(path.read_text());sequence=data.get('sequence') or ''
                if sequence and hashlib.sha256(sequence.encode()).hexdigest()==(data.get('metadata') or {}).get('sequence_sha256'):
                    return [{'sequence':sequence,'sequence_sha256':hashlib.sha256(sequence.encode()).hexdigest()}],'verified_structure_and_saved_input_sha256'
            except (OSError,ValueError,TypeError):pass
        prep=(run.get('manual_input') or {}).get('receptor_preparation') or {}
        if depth==0 and run.get('mode')=='vina' and prep.get('preparation_id'):
            source=self.get(str(prep.get('source_run_id') or ''));preparation=self.get(prep['preparation_id'])
            metadata=(preparation or {}).get('result') or {}
            if source and preparation and preparation['mode']=='vina_receptor' and preparation['status']=='succeeded' and all(prep.get(k) and prep[k]==metadata.get(k) for k in ['source_run_id','source_sha256','receptor_sha256']):
                records,_=self._run_sequence_records(source,depth=1)
                return records,'server_receptor_preparation_and_source_sequence'
        return [],'unverified'
    def _run_sequence_identities(self,run,depth=0):
        records,source=self._run_sequence_records(run,depth)
        return list(dict.fromkeys(row['sequence_sha256'] for row in records)),source
    @staticmethod
    def _evidence_payload(run,candidate,identity_source):
        payload={'candidate_id':candidate['id'],'run_id':run['id'],'sequence_sha256':candidate['sequence_sha256'],
                 'status':run['status'],'verified_by':identity_source,'created_at':stamp(),
                 'workflow':(run.get('execution_spec') or {}).get('workflow') or run.get('mode'),
                 'reported_answer':str((run.get('result') or {}).get('final_text') or '')[:6000],
                 'operations':[{'tool':op['tool'],'status':op['status'],'facts':op['facts']} for op in run.get('operations',[])],
                 'provenance':'computational_result_not_experimental_validation'}
        review=CandidateStore._docking_review(run)
        if review:payload['docking_review']=review
        return payload

    @staticmethod
    def _docking_review(run):
        if run.get('mode') != 'vina':
            return None
        vina=(run.get('result') or {}).get('vina') or {}
        groups=[]
        for group in vina.get('top_pose_consistency') or []:
            groups.append({key:group.get(key) for key in [
                'ligand_id','name','seed_count','status','max_top_pose_rmsd_angstrom',
                'top_score_range_kcal_mol','central_top_pose_seed','recommended_review_task_id',
                'recommended_review_seed','review_status','review_action','scientific_limit']})
        statuses=[group.get('review_status') for group in groups]
        if 'divergent_review_required' in statuses:
            review_status='divergent_review_required'
        elif 'analysis_unavailable' in statuses:
            review_status='analysis_unavailable'
        elif statuses and all(status=='consistent_review_required' for status in statuses):
            review_status='consistent_review_required'
        else:
            review_status='single_seed_review_required'
        return {'review_status':review_status,'review_scope':'computational_pose_review_not_experimental_validation',
                'ligand_count':vina.get('ligand_count'),'task_count':vina.get('task_count',len(vina.get('results',[]))),
                'configuration':vina.get('configuration'),'groups':groups}

    def link_candidate_evidence(self,pid,key,run_id):
        run=self.get(run_id)
        if not run:raise HTTPException(404,'任务不存在')
        context=self.run_project_context(run_id)
        current=self.conversation_project(run.get('conversation_id')) if run.get('conversation_id') else None
        if not ((context and context['project']['id']==pid) or (current and current['id']==pid)):
            raise HTTPException(422,'请先将该任务所在对话关联到本项目')
        if run['status'] not in {'succeeded','partial'}:raise HTTPException(422,'任务尚无可用结果')
        identities, identity_source = self._run_sequence_identities(run)
        with self._lock,self._connect() as db:
            db.execute('BEGIN IMMEDIATE');self._owned_project(db,pid)
            row=db.execute('SELECT payload FROM protein_candidates WHERE project_id=? AND candidate_key=?',(pid,key)).fetchone()
            if not row:raise HTTPException(404,'候选不存在')
            candidate=json.loads(row[0])
            if candidate['sequence_sha256'] not in identities:raise HTTPException(422,'无法从任务持久化输入验证相同序列，未建立关联')
            prior=db.execute('SELECT payload FROM candidate_evidence WHERE candidate_id=? AND run_id=?',(candidate['id'],run_id)).fetchone()
            if prior:return json.loads(prior[0])
            payload=self._evidence_payload(run,candidate,identity_source)
            rev=self._next_ledger_revision(db,pid)
            db.execute('INSERT INTO candidate_evidence VALUES (?,?,?,?,?)',(pid,candidate['id'],run_id,rev,dump(payload)))
            return payload
    def register_candidate_from_run(self,pid,request):
        run=self.get(request.run_id)
        if not run:raise HTTPException(404,'任务不存在')
        context=self.run_project_context(request.run_id)
        current=self.conversation_project(run.get('conversation_id')) if run.get('conversation_id') else None
        if not ((context and context['project']['id']==pid) or (current and current['id']==pid)):
            raise HTTPException(422,'请先将该任务所在对话关联到本项目')
        if run['status'] not in {'succeeded','partial'}:raise HTTPException(422,'任务尚无可用结果')
        records,identity_source=self._run_sequence_records(run)
        unique={row['sequence_sha256']:row for row in records}
        if len(unique)!=1:raise HTTPException(422,'任务必须追溯到唯一明确的蛋白序列，才能直接登记候选')
        record=next(iter(unique.values()))
        with self._lock,self._connect() as db:
            db.execute('BEGIN IMMEDIATE');self._owned_project(db,pid)
            if db.execute('SELECT 1 FROM protein_candidates WHERE project_id=? AND candidate_key=?',(pid,request.key)).fetchone():
                raise HTTPException(409,'候选编号已存在；请在候选下方关联计算记录')
            data=request.model_dump(by_alias=True);data.pop('run_id',None)
            candidate={**data,'id':str(uuid.uuid4()),'created_at':stamp(),'sequence':record['sequence'],
                       'sequence_sha256':record['sequence_sha256'],'length':len(record['sequence']),
                       'ambiguous_sequence':bool(re.search('[BXZJUO]',record['sequence'])),
                       'provenance':'server_extracted_sequence_and_user_registered_source'}
            rev=self._next_ledger_revision(db,pid)
            db.execute('INSERT INTO protein_candidates VALUES (?,?,?,?,?)',(candidate['id'],pid,candidate['key'],rev,dump(candidate)))
            evidence=self._evidence_payload(run,candidate,identity_source)
            db.execute('INSERT INTO candidate_evidence VALUES (?,?,?,?,?)',(pid,candidate['id'],run['id'],rev,dump(evidence)))
            return {'candidate':candidate,'evidence':evidence,'revision':rev}
    def _report_ledger(self,pid):
        first=self.candidate_ledger(pid,offset=0,limit=100)
        candidates=list(first['candidates'])
        while len(candidates)<first['total']:
            page=self.candidate_ledger(pid,offset=len(candidates),limit=100)
            candidates.extend(page['candidates'])
        return {**first,'candidates':candidates}
    def project_report(self,pid):
        ledger=self._report_ledger(pid)
        with self._connect() as db:project=self._owned_project(db,pid)
        candidates=[]
        for candidate in ledger['candidates']:
            evidence=[]
            for link in candidate.get('evidence',[]):
                evidence.append({key:link.get(key) for key in ['run_id','workflow','status','verified_by','created_at','docking_review'] if link.get(key) is not None})
            candidates.append({key:candidate.get(key) for key in ['id','key','name','length','sequence_sha256','source','construct',
                'ambiguous_sequence','measurement_count','docking_review_count','docking_review_status','review_state',
                'evidence_kinds','missing_evidence','recommended_next']})
            candidates[-1]['evidence']=evidence
        payload={'report_version':'project-review-report-v1','generated_at':stamp(),'project':project,
                 'ledger_revision':ledger['revision'],'review_queue':ledger.get('review_queue'),
                 'candidates':candidates,'measurements':ledger['measurements'],'comparisons':ledger['comparisons'],
                 'imports':ledger['imports'],
                 'limitations':['计算结果、对接一致性和用户导入的实验记录均按其来源标记，不等于独立实验验证。',
                                '报告默认不包含完整蛋白序列；序列身份由长度和 SHA-256 核对。',
                                '实验比较仅在条件完整时进行；缺失值、不同底物、批次或构建体不会被静默合并。']}
        payload['snapshot_sha256']=hashlib.sha256(dump(payload).encode()).hexdigest()
        return payload
    @staticmethod
    def project_report_markdown(report):
        project=report['project'];queue=report.get('review_queue') or {}
        lines=['# '+project['name']+'：项目复核报告','',
               f"- 项目 ID：`{project['id']}`",f"- 项目版本：{project['revision']}",
               f"- 台账版本：{report['ledger_revision']}",f"- 生成时间：{report['generated_at']}",
               f"- 快照 SHA-256：`{report['snapshot_sha256']}`",'',
               '## 研究背景','',project.get('objective') or '未填写','',
               '**条件与约束**：'+ (project.get('constraints') or '未填写'),'',
               '**待验证假设**：'+ (project.get('hypotheses') or '未填写'),'',
               '## 复核队列','',f"完整项目快照共 {queue.get('total',len(report['candidates']))} 个候选。",'']
        for state,count in (queue.get('counts') or {}).items():lines.append(f'- {state}：{count}')
        lines.extend(['','## 候选与证据',''])
        for candidate in report['candidates']:
            lines.extend([f"### {candidate['name']}（{candidate['key']}）",'',
                f"- 序列：{candidate['length']} aa；SHA-256 `{candidate['sequence_sha256']}`",
                f"- 来源：{candidate['source']}",f"- 复核状态：{candidate.get('review_state') or '--'}",
                f"- 证据类型：{', '.join(candidate.get('evidence_kinds') or []) or '无'}",
                f"- 待补证据：{'；'.join(candidate.get('missing_evidence') or []) or '无'}",
                f"- 建议：{candidate.get('recommended_next') or '--'}",''])
            for link in candidate.get('evidence',[]):
                lines.append(f"- 计算记录 `{link['run_id']}` · {link.get('workflow') or '--'} · {link.get('verified_by') or '--'}")
                review=link.get('docking_review') or {}
                for group in review.get('groups',[]):
                    lines.append(f"  - {group.get('name')}：中心种子 {group.get('recommended_review_seed','--')}；最大首位 RMSD {group.get('max_top_pose_rmsd_angstrom','--')} Å；分数跨度 {group.get('top_score_range_kcal_mol','--')} kcal/mol；{group.get('review_status','--')}")
            lines.append('')
        lines.extend(['## 实验记录',''])
        if not report['measurements']:lines.append('尚无实验记录。')
        for row in report['measurements']:
            lines.append(f"- {row.get('candidate_key')} · {row.get('metric')}：{row.get('value') if row.get('value') is not None else '未测定'} {row.get('unit') or ''}；条件 `{dump(row.get('conditions') or {})}`")
        lines.extend(['','## 条件完整性与比较',''])
        for group in report['comparisons']:
            lines.append(f"- {group.get('metric')}：{'条件完整' if group.get('comparable') else '需核对条件'}；问题：{'；'.join(group.get('issues') or []) or '无'}")
        lines.extend(['','## 原始文件',''])
        for item in report['imports']:lines.append(f"- {item['filename']} · SHA-256 `{item['file_sha256']}`")
        lines.extend(['','## 适用边界','']+['- '+item for item in report['limitations']])
        return '\n'.join(lines)+'\n'

    def preview_experiments(self,pid,upload):
        with self._lock,self._connect() as db:
            db.execute('BEGIN IMMEDIATE');self._owned_project(db,pid)
            try:raw,rows=read_table(upload.filename,upload.content_base64)
            except ValueError as e:raise HTTPException(422,str(e))
            candidates={c['key']:c for c in [json.loads(r[0]) for r in db.execute('SELECT payload FROM protein_candidates WHERE project_id=?',(pid,))]}
            parsed=[];errors=[];seen=set()
            prior={r[0] for r in db.execute('SELECT row_sha256 FROM experimental_measurements WHERE project_id=?',(pid,))}
            for saved in db.execute('SELECT payload FROM experimental_measurements WHERE project_id=?',(pid,)):
                data=json.loads(saved[0]);seen.add(dump({k:data[k] for k in ['candidate_id','metric','data_kind','replicate_id','conditions']}))
            for line,row in enumerate(rows,2):
                try:
                    data=normalize_row(row,candidates);digest=row_hash(data)
                    # Ignore value for identity, so the same replicate/summary cannot be imported twice with a different result.
                    identity=dump({k:data[k] for k in ['candidate_id','metric','data_kind','replicate_id','conditions']})
                    if identity in seen or digest in prior:raise ValueError('重复测量记录，请勿重复导入')
                    seen.add(identity);parsed.append(data)
                except (ValueError,TypeError) as e:errors.append({'row':line,'message':str(e)})
            digest=hashlib.sha256(raw).hexdigest()
            existing=db.execute('SELECT id FROM experiment_imports WHERE project_id=? AND file_sha256=?',(pid,digest)).fetchone()
            if existing:raise HTTPException(409,'此文件已导入')
            token=str(uuid.uuid4());rev=self._ledger_revision(db,pid)
            report={'token':token,'revision':rev,'filename':upload.filename,'file_sha256':digest,'row_count':len(rows),
                    'valid_count':len(parsed),'errors':errors[:100],'error_count':len(errors),'can_commit':not errors,
                    'preview':parsed[:10],'comparisons':comparisons(parsed),'provenance':'uploaded_data_not_independently_verified'}
            # Keep only a bounded number of uncommitted raw previews per project.
            db.execute('DELETE FROM experiment_previews WHERE project_id=? AND committed_id IS NULL AND id NOT IN (SELECT id FROM experiment_previews WHERE project_id=? AND committed_id IS NULL ORDER BY created_at DESC LIMIT 4)',(pid,pid))
            db.execute('INSERT INTO experiment_previews VALUES (?,?,?,?,?,?,?,?,NULL)',(token,pid,rev,upload.filename,digest,raw,dump({'report':report,'rows':parsed}),stamp()))
            return report
    def commit_experiments(self,pid,token):
        with self._lock,self._connect() as db:
            db.execute('BEGIN IMMEDIATE');self._owned_project(db,pid)
            preview=db.execute('SELECT * FROM experiment_previews WHERE id=? AND project_id=?',(token,pid)).fetchone()
            if not preview:raise HTTPException(404,'预览不存在，请重新上传')
            if preview['committed_id']:return {'import_id':preview['committed_id'],'already_committed':True}
            if preview['revision']!=self._ledger_revision(db,pid):raise HTTPException(409,'候选或实验数据已变化，请重新预览')
            payload=json.loads(preview['payload'])
            if not payload['report']['can_commit']:raise HTTPException(422,'请修正所有错误后重新上传')
            rev=self._next_ledger_revision(db,pid);iid=str(uuid.uuid4())
            db.execute('INSERT INTO experiment_imports VALUES (?,?,?,?,?,?,?)',(iid,pid,preview['file_sha256'],rev,preview['filename'],preview['raw'],stamp()))
            for row in payload['rows']:
                row['import_id']=iid
                db.execute('INSERT INTO experimental_measurements VALUES (?,?,?,?,?)',(pid,iid,rev,row_hash({k:v for k,v in row.items() if k!='import_id'}),dump(row)))
            db.execute('UPDATE experiment_previews SET committed_id=? WHERE id=?',(iid,token))
            return {'import_id':iid,'row_count':len(payload['rows']),'revision':rev,'already_committed':False}


def candidate_tools(store,snapshot):
    from langchain_core.tools import tool
    ledger=snapshot.get('candidate_ledger')
    if not ledger:return []
    @tool
    def inspect_candidate_evidence(candidate_key: str = '', offset: int = 0) -> str:
        """Read immutable project candidate identities and imported experiments, 20 candidates per page.
        Reports missing conditions and comparability. Imported measurements are not independently verified.
        Use this for candidate comparison and experimental feedback; it launches no scientific computation.
        """
        data=store.candidate_ledger(snapshot['project']['id'],revision=ledger['revision'],candidate_key=candidate_key or None,offset=max(0,offset),limit=20)
        for c in data['candidates']:
            c.pop('sequence',None)
            c['older_evidence_omitted']=max(0,len(c['evidence'])-10)
            c['evidence']=c['evidence'][-10:]
        keys={c['key'] for c in data['candidates']}
        groups=[g for g in data['comparisons'] if any(v['candidate_key'] in keys for v in g['values'])]
        data['comparison_groups_omitted']=len(groups)>40
        data['comparisons']=groups[:40]
        for group in data['comparisons']:
            all_values=group['values'];group['values']=[v for v in all_values if v['candidate_key'] in keys]
            group['exact_values_not_displayed_due_to_filter_or_pagination']=len(all_values)-len(group['values'])
            group['additional_issues_omitted']=len(group['issues'])>30
            group['issues']=group['issues'][:30]
        data['additional_measurements_omitted']=len(data['measurements'])>200
        data['measurements']=data['measurements'][:200]
        data['imports']=data['imports'][-20:]
        data['next_offset']=offset+20 if offset+20<data['total'] else None
        return dump(data)
    return [inspect_candidate_evidence]


def register_candidate_routes(app,get_store):
    @app.get('/projects/{pid}/candidates')
    def get_candidates(pid:str,offset:int=0,candidate_key:str=''):
        return get_store().candidate_ledger(pid,offset=max(0,offset),candidate_key=candidate_key or None)
    @app.post('/projects/{pid}/candidates',status_code=201)
    def add_candidates(pid:str,batch:CandidateBatch):return get_store().register_candidates(pid,batch)
    @app.post('/projects/{pid}/candidates/from-run',status_code=201)
    def add_candidate_from_run(pid:str,request:CandidateFromRun):return get_store().register_candidate_from_run(pid,request)
    @app.get('/projects/{pid}/report')
    def report(pid:str,format:str='markdown'):
        report=get_store().project_report(pid)
        if format=='json':return JSONResponse(content=report,headers={'Content-Disposition':f'attachment; filename="project-report-{pid}.json"','X-Report-SHA256':report['snapshot_sha256']})
        if format!='markdown':raise HTTPException(422,'报告格式仅支持 markdown 或 json')
        return Response(content=get_store().project_report_markdown(report),media_type='text/markdown; charset=utf-8',
                        headers={'Content-Disposition':f'attachment; filename="project-report-{pid}.md"',
                                 'X-Report-SHA256':report['snapshot_sha256']})
    @app.post('/projects/{pid}/candidates/{key}/evidence')
    def link_evidence(pid:str,key:str,link:EvidenceLink):return get_store().link_candidate_evidence(pid,key,link.run_id)
    @app.get('/projects/{pid}/experiments/template')
    def template(pid:str):
        with get_store()._connect() as db:get_store()._owned_project(db,pid)
        out=io.StringIO();csv.writer(out).writerow(COLUMNS)
        return Response(out.getvalue(),media_type='text/csv',headers={'Content-Disposition':'attachment; filename="experiments.csv"'})
    @app.post('/projects/{pid}/experiments/preview')
    def preview(pid:str,upload:ImportInput):return get_store().preview_experiments(pid,upload)
    @app.post('/projects/{pid}/experiments/commit')
    def commit(pid:str,request:CommitInput):return get_store().commit_experiments(pid,request.token)
    @app.get('/projects/{pid}/experiments/{import_id}/file')
    def original(pid:str,import_id:str):
        with get_store()._connect() as db:
            get_store()._owned_project(db,pid)
            row=db.execute('SELECT raw,filename FROM experiment_imports WHERE id=? AND project_id=?',(import_id,pid)).fetchone()
            if not row:raise HTTPException(404,'文件不存在')
            extension='xlsx' if row['filename'].lower().endswith('.xlsx') else 'csv'
            return Response(row['raw'],media_type='application/octet-stream',headers={'Content-Disposition':f'attachment; filename="experiment-{import_id}.{extension}"'})
