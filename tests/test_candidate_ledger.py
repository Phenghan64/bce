import base64,csv,io,json,tempfile,unittest
from pathlib import Path
from fastapi import HTTPException
from tests.test_vina_conversation import api
from marine_agent.access_control import account_scope
from marine_agent.research_projects import ProjectBrief,project_tools
from marine_agent.candidate_ledger import CandidateInput,CandidateBatch,ImportInput,CandidateFromRun
from marine_agent.experimental_data import COLUMNS,comparisons,read_table

class CandidateTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.store=api.RunStore(Path(self.tmp.name)/'state.sqlite3')
        self.scope=account_scope('alice');self.scope.__enter__();self.addCleanup(self.scope.__exit__,None,None,None)
        self.pid=self.store.create_project(ProjectBrief(name='候选验收'))['id']
        self.cid=self.store.new_project_conversation(self.pid)['conversation_id']
        self.items=[CandidateInput(key='A',name='同名蛋白',sequence='ACDEFGHIKLMNPQRSTVWY',source='paper'),CandidateInput(key='B',name='同名蛋白',sequence='ACDEFGHIKLMNPQRSTVWA',source='paper')]
        self.store.register_candidates(self.pid,CandidateBatch(candidates=self.items))
    def row(self,**changes):
        r=dict.fromkeys(COLUMNS,'');r.update(candidate_key='A',metric='kcat',value='60',unit='min^-1',data_kind='summary',sd='6',n='3',substrate='pNPG',temperature_c='25',ph='7.5',protocol='protocol-1',control='WT',batch='batch-1',construct='C-His',source='doi:test');r.update(changes);return r
    def upload(self,rows):
        out=io.StringIO();w=csv.DictWriter(out,fieldnames=COLUMNS);w.writeheader();w.writerows(rows)
        return ImportInput(filename='experiment.csv',content_base64=base64.b64encode(out.getvalue().encode()).decode())
    def preview(self,rows):return self.store.preview_experiments(self.pid,self.upload(rows))
    def commit(self,rows):
        preview=self.preview(rows);self.assertTrue(preview['can_commit'],preview['errors'])
        return self.store.commit_experiments(self.pid,preview['token'])
    def test_sequence_immutable_names_not_identity_and_retry(self):
        self.assertEqual(self.store.register_candidates(self.pid,CandidateBatch(candidates=self.items))['created'],0)
        with self.assertRaises(HTTPException):self.store.register_candidates(self.pid,CandidateBatch(candidates=[CandidateInput(key='A',name='同名蛋白',sequence='AAAA',source='paper')]))
        self.assertEqual(self.store.candidate_ledger(self.pid)['total'],2)
    def test_batch_atomic_on_collision(self):
        with self.assertRaises(HTTPException):self.store.register_candidates(self.pid,CandidateBatch(candidates=[CandidateInput(key='C',name='c',sequence='AAAA',source='paper'),CandidateInput(key='A',name='a',sequence='AAAA',source='paper')]))
        self.assertEqual(self.store.candidate_ledger(self.pid)['total'],2)
    def test_unit_normalization_preserves_sd_and_n(self):
        self.commit([self.row(),self.row(candidate_key='B',value='2',unit='s^-1',sd='.2')])
        data=self.store.candidate_ledger(self.pid);a=next(r for r in data['measurements'] if r['candidate_key']=='A')
        self.assertEqual(a['value'],1);self.assertEqual(a['sd'],.1);self.assertEqual(a['n'],3);self.assertEqual(a['raw_unit'],'min^-1')
        self.assertTrue(data['comparisons'][0]['comparable'])
    def test_unknown_temperature_blocks_comparison(self):
        p=self.preview([self.row(temperature_c=''),self.row(candidate_key='B',temperature_c='')]);self.assertTrue(p['can_commit']);self.assertFalse(p['comparisons'][0]['comparable'])
    def test_different_substrate_batch_construct_never_pooled(self):
        p=self.preview([self.row(),self.row(candidate_key='B',substrate='collagen'),self.row(candidate_key='B',batch='b2'),self.row(candidate_key='B',construct='untagged')]);self.assertEqual(len(p['comparisons']),4)
    def test_missing_and_censor_are_not_zero(self):
        self.commit([self.row(value='NA',sd='',n=''),self.row(candidate_key='B',qualifier='<',value='10',sd='',n='')]);data=self.store.candidate_ledger(self.pid)
        self.assertIsNone(next(r for r in data['measurements'] if r['candidate_key']=='A')['value']);self.assertFalse(data['comparisons'][0]['comparable']);self.assertEqual(data['comparisons'][0]['values'],[])
    def test_summary_does_not_generate_replicates(self):
        p=self.preview([self.row(),self.row(candidate_key='B')]);self.assertEqual(p['valid_count'],2);self.assertEqual(len(p['comparisons'][0]['values']),2)
    def test_raw_replicates_aggregate_separately(self):
        p=self.preview([self.row(value=str(v),sd='',n='',data_kind='replicate',replicate_id=str(i)) for i,v in enumerate([60,120,180])]);v=p['comparisons'][0]['values'][0];self.assertEqual(v['n'],3);self.assertEqual(v['value'],2);self.assertEqual(v['sd'],1)
    def test_validation_errors_block_all_commit(self):
        for changes in [{'candidate_key':'unknown'},{'unit':'U/mg'},{'value':'nan'},{'value':'inf'},{'ph':'15'},{'sd':'-1'},{'n':'1.5'},{'data_kind':'replicate'},{'qualifier':'missing'},{'source':''}]:
            with self.subTest(changes=changes):
                p=self.preview([self.row(),self.row(candidate_key='B',**changes)] if 'candidate_key' not in changes else [self.row(**changes)])
                self.assertFalse(p['can_commit'])
                with self.assertRaises(HTTPException):self.store.commit_experiments(self.pid,p['token'])
        self.assertEqual(self.store.candidate_ledger(self.pid)['measurements'],[])
    def test_repeated_rows_values_and_files_rejected(self):
        p=self.preview([self.row(),self.row(value='70')]);self.assertFalse(p['can_commit'])
        self.commit([self.row()])
        with self.assertRaises(HTTPException):self.preview([self.row()])
        p=self.preview([self.row(value='80')]);self.assertFalse(p['can_commit'])
    def test_commit_idempotence_and_restart(self):
        p=self.preview([self.row()]);first=self.store.commit_experiments(self.pid,p['token']);second=self.store.commit_experiments(self.pid,p['token']);self.assertEqual(first['import_id'],second['import_id'])
        reopened=api.RunStore(self.store.path);self.assertEqual(len(reopened.candidate_ledger(self.pid)['measurements']),1)
    def test_preview_stale_after_another_import(self):
        old=self.preview([self.row()]);self.commit([self.row(candidate_key='B')])
        with self.assertRaises(HTTPException) as e:self.store.commit_experiments(self.pid,old['token'])
        self.assertEqual(e.exception.status_code,409)
    def test_owner_isolation_all_paths(self):
        p=self.preview([self.row()])
        with account_scope('bob'):
            for fn in [lambda:self.store.candidate_ledger(self.pid),lambda:self.store.register_candidates(self.pid,CandidateBatch(candidates=self.items)),lambda:self.store.preview_experiments(self.pid,self.upload([self.row()])),lambda:self.store.commit_experiments(self.pid,p['token'])]:
                with self.assertRaises(HTTPException) as e:fn()
                self.assertEqual(e.exception.status_code,404)
    def test_immutable_agent_snapshot_not_new_computation(self):
        before=self.store.create('比较候选','general',conversation_id=self.cid)
        self.commit([self.row()]);after=self.store.create('读取新实验','general',conversation_id=self.cid)
        a=json.loads(project_tools(self.store,before['run_id'])[1].invoke({}));b=json.loads(project_tools(self.store,after['run_id'])[1].invoke({}))
        self.assertEqual(a['measurements'],[]);self.assertEqual(len(b['measurements']),1)
        self.assertNotIn('sequence',b['candidates'][0])
        from marine_agent.result_contract import execution_receipt
        self.assertEqual(execution_receipt({'id':'x','status':'succeeded','mode':'general','operations':[],'result':{'agent':{'messages':[{'type':'tool','name':'inspect_candidate_evidence','content':json.dumps(b)}]}}})['kind'],'dialogue')
    def test_link_requires_project_and_exact_persisted_sequence(self):
        r=self.store.create('预测','general',conversation_id=self.cid)
        self.store.update(r['run_id'],status='succeeded',execution_json=json.dumps({'workflow':'openfold_sequence','input':{'sequence':self.items[0].sequence}}))
        link=self.store.link_candidate_evidence(self.pid,'A',r['run_id']);self.assertIn('sha256',link['verified_by'])
        self.assertEqual(self.store.link_candidate_evidence(self.pid,'A',r['run_id']),link)
        with self.assertRaises(HTTPException):self.store.link_candidate_evidence(self.pid,'B',r['run_id'])
        outside=self.store.create('预测','general');self.store.update(outside['run_id'],status='succeeded')
        with self.assertRaises(HTTPException):self.store.link_candidate_evidence(self.pid,'A',outside['run_id'])
    def test_vina_review_is_append_only_candidate_evidence_not_experimental_validation(self):
        run=self.store.create('多种子对接','vina',conversation_id=self.cid)
        group={'ligand_id':'ligand','name':'配体','seed_count':3,'status':'available',
               'max_top_pose_rmsd_angstrom':0.137,'top_score_range_kcal_mol':0.041,
               'central_top_pose_seed':42,'recommended_review_task_id':'ligand-seed-42',
               'recommended_review_seed':42,'review_status':'consistent_review_required',
               'review_action':'复核中心姿势','scientific_limit':'一致性不等于实验正确性'}
        result={'vina':{'status':'succeeded','ligand_count':1,'task_count':3,
                        'configuration':{'params':{'seeds':[42,2026,20261021]}},
                        'top_pose_consistency':[group]},'final_text':'对接完成'}
        self.store.update(run['run_id'],status='succeeded',result_json=json.dumps(result),
            execution_json=json.dumps({'workflow':'vina','input':{'sequence':self.items[0].sequence}}))
        link=self.store.link_candidate_evidence(self.pid,'A',run['run_id'])
        self.assertEqual(link['docking_review']['review_status'],'consistent_review_required')
        self.assertEqual(link['docking_review']['groups'][0]['recommended_review_seed'],42)
        self.assertEqual(link['provenance'],'computational_result_not_experimental_validation')
        candidate=self.store.candidate_ledger(self.pid)['candidates'][0]
        self.assertEqual(candidate['docking_review_count'],1)
        self.assertEqual(candidate['docking_review_status'],'consistent_review_required')
        self.assertEqual(self.store.link_candidate_evidence(self.pid,'A',run['run_id']),link)

    def test_vina_identity_follows_verified_receptor_preparation_back_to_source_sequence(self):
        source=self.store.create('结构预测','general',conversation_id=self.cid)
        self.store.update(source['run_id'],status='succeeded',execution_json=json.dumps({'workflow':'openfold_sequence','input':{'sequence':self.items[0].sequence}}))
        prep=self.store.create('受体准备','vina_receptor')
        metadata={'source_run_id':source['run_id'],'source_sha256':'source-hash','receptor_sha256':'receptor-hash'}
        self.store.update(prep['run_id'],status='succeeded',result_json=json.dumps(metadata))
        vina=self.store.create('对接','vina',conversation_id=self.cid,manual={'receptor_preparation':{'preparation_id':prep['run_id'],**metadata}})
        self.store.update(vina['run_id'],status='succeeded',result_json=json.dumps({'vina':{'results':[]}}))
        link=self.store.link_candidate_evidence(self.pid,'A',vina['run_id'])
        self.assertEqual(link['verified_by'],'server_receptor_preparation_and_source_sequence')
        self.assertEqual(link['docking_review']['review_status'],'single_seed_review_required')

    def test_annotation_and_blast_only_match_query_not_reference(self):
        for workflow in ['protein_domains','protein_screen']:
            r=self.store.create('功能证据','general',conversation_id=self.cid)
            self.store.update(r['run_id'],status='succeeded',execution_json=json.dumps({'workflow':workflow,'input':{'query_fasta':'>query\n'+self.items[0].sequence,'reference_fasta':'>reference\n'+self.items[1].sequence}}))
            self.store.link_candidate_evidence(self.pid,'A',r['run_id'])
            with self.assertRaises(HTTPException):self.store.link_candidate_evidence(self.pid,'B',r['run_id'])
    def test_direct_candidate_registration_uses_unique_saved_sequence_and_one_transaction(self):
        source=self.store.create('结构预测','general',conversation_id=self.cid)
        self.store.update(source['run_id'],status='succeeded',execution_json=json.dumps({'workflow':'openfold_sequence','input':{'sequence':self.items[0].sequence}}))
        prep=self.store.create('受体准备','vina_receptor');metadata={'source_run_id':source['run_id'],'source_sha256':'source-hash','receptor_sha256':'receptor-hash'}
        self.store.update(prep['run_id'],status='succeeded',result_json=json.dumps(metadata))
        vina=self.store.create('对接','vina',conversation_id=self.cid,manual={'receptor_preparation':{'preparation_id':prep['run_id'],**metadata}})
        self.store.update(vina['run_id'],status='succeeded',result_json=json.dumps({'vina':{'results':[]}}))
        request=CandidateFromRun(run_id=vina['run_id'],key='C',name='从任务登记',source='verified-run',construct='WT')
        result=self.store.register_candidate_from_run(self.pid,request)
        self.assertEqual(result['candidate']['sequence'],self.items[0].sequence)
        self.assertEqual(result['evidence']['verified_by'],'server_receptor_preparation_and_source_sequence')
        self.assertEqual(result['evidence']['docking_review']['review_status'],'single_seed_review_required')
        self.assertEqual(self.store.candidate_ledger(self.pid)['revision'],2)
        with self.assertRaises(HTTPException):self.store.register_candidate_from_run(self.pid,request)
        standalone=self.store.create('独立上传','vina',conversation_id=self.cid,manual={})
        self.store.update(standalone['run_id'],status='succeeded',result_json=json.dumps({'vina':{'results':[]}}))
        with self.assertRaises(HTTPException):self.store.register_candidate_from_run(self.pid,CandidateFromRun(run_id=standalone['run_id'],key='D',name='无来源',source='x'))
        ledger=self.store.candidate_ledger(self.pid)
        self.assertEqual(ledger['review_queue']['scope'],'full_snapshot')
        self.assertEqual(ledger['review_queue']['total'],3)
        self.assertEqual(ledger['review_queue']['counts']['needs_computation'],2)
        self.assertEqual(ledger['review_queue']['counts']['needs_multiseed_review'],1)
        candidate=next(c for c in ledger['candidates'] if c['key']=='C')
        self.assertEqual(candidate['review_state'],'needs_multiseed_review')
        self.assertIn('尚无实验记录',candidate['missing_evidence'])

    def test_project_report_is_auditable_and_preserves_scientific_boundaries(self):
        self.commit([self.row()])
        run=self.store.create('结构预测','general',conversation_id=self.cid)
        self.store.update(run['run_id'],status='succeeded',execution_json=json.dumps({'workflow':'openfold_sequence','input':{'sequence':self.items[0].sequence}}),result_json=json.dumps({'final_text':'结构预测完成'}))
        self.store.link_candidate_evidence(self.pid,'A',run['run_id'])
        report=self.store.project_report(self.pid)
        self.assertEqual(report['report_version'],'project-review-report-v1')
        self.assertEqual(report['review_queue']['total'],2)
        self.assertEqual(report['review_queue']['scope'],'full_snapshot')
        self.assertEqual(len(report['candidates']),2)
        self.assertEqual(report['candidates'][0]['evidence'][0]['run_id'],run['run_id'])
        self.assertNotIn('sequence',report['candidates'][0])
        expected=__import__('hashlib').sha256(__import__('marine_agent.candidate_ledger',fromlist=['dump']).dump({k:v for k,v in report.items() if k!='snapshot_sha256'}).encode()).hexdigest()
        self.assertEqual(report['snapshot_sha256'],expected)
        text=self.store.project_report_markdown(report)
        self.assertIn('项目复核报告',text)
        self.assertIn(run['run_id'],text)
        self.assertIn('不等于独立实验验证',text)
        self.assertIn('doi:test',text)

    def test_report_routes_return_downloadable_json_and_markdown(self):
        from fastapi import FastAPI
        from fastapi.testclient import TestClient
        from marine_agent.candidate_ledger import register_candidate_routes
        app=FastAPI();register_candidate_routes(app,lambda:self.store)
        client=TestClient(app)
        response=client.get('/projects/'+self.pid+'/report?format=json')
        self.assertEqual(response.status_code,200)
        self.assertIn('attachment',response.headers['content-disposition'])
        self.assertEqual(response.json()['report_version'],'project-review-report-v1')
        markdown=client.get('/projects/'+self.pid+'/report?format=markdown')
        self.assertIn('attachment',markdown.headers['content-disposition'])
        self.assertIn('项目复核报告',markdown.text)
        self.assertEqual(client.get('/projects/'+self.pid+'/report?format=xml').status_code,422)

    def test_project_task_menu_distinguishes_reading_and_computation(self):
        dialogue=self.store.create('读取实验','general',conversation_id=self.cid);self.store.update(dialogue['run_id'],status='succeeded')
        computed=self.store.create('注释','general',conversation_id=self.cid);self.store.update(computed['run_id'],status='succeeded',execution_json=json.dumps({'workflow':'protein_domains','input':{'query_fasta':'>x\nAAAA'}}))
        items={r['id']:r for r in self.store.project_detail(self.pid)['recent_tasks']}
        self.assertFalse(items[dialogue['run_id']]['has_scientific_input']);self.assertTrue(items[computed['run_id']]['has_scientific_input'])
    def test_historical_task_can_be_explicitly_brought_into_project(self):
        from marine_agent.research_projects import ProjectLink
        old=self.store.create('历史预测','general')
        self.store.update(old['run_id'],status='succeeded',execution_json=json.dumps({'workflow':'openfold_sequence','input':{'sequence':self.items[0].sequence}}))
        self.store.link_project(old['conversation_id'],ProjectLink(project_id=self.pid))
        self.assertEqual(self.store.link_candidate_evidence(self.pid,'A',old['run_id'])['run_id'],old['run_id'])
        self.assertIsNone(self.store.run_project_context(old['run_id']))
    def test_candidate_followup_uses_evidence_agent_even_when_planner_says_conversation(self):
        from unittest.mock import patch
        from marine_agent.turn_planning import TurnPlan
        self.commit([self.row()])
        first=self.store.create('读取实验','general',conversation_id=self.cid);self.store.update(first['run_id'],status='succeeded')
        second=self.store.create('缺失数据是0活性吗','general',conversation_id=self.cid)
        seen=[]
        class Agent:
            def stream(self,payload,**kwargs):
                seen.append(payload);yield {'messages':[{'type':'ai','content':'未测量不是零活性。'}]}
        with patch.object(api,'store',self.store),patch.object(api,'plan_turn',return_value=TurnPlan(action='inspect',reason='实验追问',context_scope='conversation',inspection='activity_evidence')),patch.object(api,'build_agent',return_value=Agent()),patch.object(api,'_reply_task_status') as shortcut:
            api._controlled_execute(second['run_id'],api._execute_run,'缺失数据是0活性吗','general')
            shortcut.assert_not_called();self.assertTrue(seen)
    def test_xlsx_value_roundtrip_formula_rejected(self):
        from openpyxl import Workbook
        book=Workbook();sheet=book.active;sheet.append(COLUMNS);sheet.append([self.row()[k] for k in COLUMNS]);out=io.BytesIO();book.save(out)
        p=self.store.preview_experiments(self.pid,ImportInput(filename='experiment.xlsx',content_base64=base64.b64encode(out.getvalue()).decode()));self.assertTrue(p['can_commit'])
        sheet['C2']='=1+2';out=io.BytesIO();book.save(out)
        with self.assertRaises(HTTPException):self.store.preview_experiments(self.pid,ImportInput(filename='experiment.xlsx',content_base64=base64.b64encode(out.getvalue()).decode()))
    def test_bounds_bad_headers_and_workbooks(self):
        for filename,raw in [('x.csv',b'a,a\n1,2'),('x.xlsx',b'not zip'),('x.exe',b'test')]:
            with self.assertRaises(ValueError):read_table(filename,base64.b64encode(raw).decode())
    def test_api_download_and_anonymous_denial(self):
        from fastapi import FastAPI
        from fastapi.testclient import TestClient
        from unittest.mock import patch
        from marine_agent import access_control
        from marine_agent.candidate_ledger import register_candidate_routes
        app=FastAPI();register_candidate_routes(app,lambda:self.store);app.add_middleware(access_control.AccountMiddleware)
        def verify(cookie):
            if cookie!='session=alice':raise HTTPException(401,'login')
            return 'alice'
        with patch.object(access_control,'verify_session',side_effect=verify),TestClient(app) as client:
            url='/projects/'+self.pid
            self.assertEqual(client.get(url+'/candidates').status_code,401)
            headers={'cookie':'session=alice'}
            p=client.post(url+'/experiments/preview',headers=headers,json=self.upload([self.row()]).model_dump());self.assertEqual(p.status_code,200)
            result=client.post(url+'/experiments/commit',headers=headers,json={'token':p.json()['token']});self.assertEqual(result.status_code,200)
            original=client.get(url+'/experiments/'+result.json()['import_id']+'/file',headers=headers)
            self.assertEqual(original.content,base64.b64decode(self.upload([self.row()]).content_base64))

if __name__=='__main__':unittest.main()
