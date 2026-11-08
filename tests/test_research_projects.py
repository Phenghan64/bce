"""Reusable project isolation, provenance, recovery and dialogue acceptance cases."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from fastapi import HTTPException
from langchain_core.messages import AIMessage
from tests.test_vina_conversation import api
from marine_agent.access_control import account_scope
from marine_agent.research_projects import ProjectBrief, ProjectUpdate, ProjectLink, project_tools
from marine_agent.scientific_knowledge import compose_explanation, references, _structured_result
from marine_agent.turn_planning import TurnPlan


class ProjectAcceptanceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.store=api.RunStore(Path(self.temp.name)/'state.sqlite3')
        self.owner=account_scope('alice');self.owner.__enter__();self.addCleanup(self.owner.__exit__,None,None,None)
        self.project=self.store.create_project(ProjectBrief(name='海洋胶原蛋白酶',objective='筛选候选，实验验证活性',hypotheses='候选A可能降解胶原，尚未验证'))
        self.cid=self.store.new_project_conversation(self.project['id'])['conversation_id']
    def test_project_cross_conversation_and_restart(self):
        first=self.store.create('已测候选A','general',conversation_id=self.cid)
        self.store.update(first['run_id'],status='succeeded',result_json=json.dumps({'final_text':'候选A尚无实验活性证据'}))
        cid=self.store.new_project_conversation(self.project['id'])['conversation_id']
        second=self.store.create('继续讨论项目','general',conversation_id=cid)
        reopened=api.RunStore(self.store.path)
        context=reopened.run_project_context(second['run_id'])
        self.assertEqual(context['project']['objective'],self.project['objective'])
        self.assertEqual(context['recent_tasks'][0]['id'],first['run_id'])
        tool=project_tools(reopened,second['run_id'])[0]
        result=json.loads(tool.invoke({'source_run_id':first['run_id']}))
        self.assertIn('尚无实验',result['reported_answer'])
        self.assertEqual(len(reopened.list_projects()),1)
        self.assertEqual(len(reopened.project_detail(self.project['id'])['conversations']),2)
    def test_snapshot_is_immutable_and_stale_edit_conflicts(self):
        run=self.store.create('任务','general',conversation_id=self.cid)
        update=ProjectUpdate(name='更新项目',objective='新目标',revision=1)
        self.assertEqual(self.store.update_project(self.project['id'],update)['revision'],2)
        self.assertEqual(self.store.run_project_context(run['run_id'])['project']['revision'],1)
        with self.assertRaises(HTTPException) as exc:self.store.update_project(self.project['id'],update)
        self.assertEqual(exc.exception.status_code,409)
    def test_accounts_cannot_read_link_update_or_inspect_foreign_projects(self):
        run=self.store.create('任务','general',conversation_id=self.cid)
        with account_scope('bob'):
            other=self.store.create('其他会话','general');self.store.update(other['run_id'],status='succeeded')
            self.assertEqual(self.store.list_projects(),[])
            self.assertIsNone(self.store.run_project_context(run['run_id']))
            for op in [lambda:self.store.project_detail(self.project['id']),lambda:self.store.new_project_conversation(self.project['id']),lambda:self.store.update_project(self.project['id'],ProjectUpdate(name='修改',revision=1)),lambda:self.store.link_project(other['conversation_id'],ProjectLink(project_id=self.project['id']))]:
                with self.assertRaises(HTTPException) as exc:op()
                self.assertEqual(exc.exception.status_code,404)
    def test_link_active_and_stale_link_fail_without_changes(self):
        other=self.store.create_project(ProjectBrief(name='另一个项目'))
        run=self.store.create('运行中','general',conversation_id=self.cid)
        with self.assertRaises(HTTPException):self.store.link_project(self.cid,ProjectLink(project_id=other['id'],expected_project_id=self.project['id']))
        self.store.update(run['run_id'],status='stopped')
        with self.assertRaises(HTTPException):self.store.link_project(self.cid,ProjectLink(project_id=other['id']))
        self.store.link_project(self.cid,ProjectLink(project_id=other['id'],expected_project_id=self.project['id']))
        self.assertEqual(self.store.conversation_project(self.cid)['id'],other['id'])
        self.assertEqual(self.store.run_project_context(run['run_id'])['project']['id'],self.project['id'])
    def test_project_read_tool_cannot_read_unlisted_tasks(self):
        unrelated=self.store.create('secret unrelated','general')
        current=self.store.create('project query','general',conversation_id=self.cid)
        self.assertEqual(json.loads(project_tools(self.store,current['run_id'])[0].invoke({'source_run_id':unrelated['run_id']}))['status'],'not_found')
    def test_project_tasks_bounded_and_omissions_explicit(self):
        for i in range(22):self.store.create(f'任务{i}','general',conversation_id=self.cid)
        final=self.store.create('最新','general',conversation_id=self.cid)
        context=self.store.run_project_context(final['run_id'])
        self.assertEqual(len(context['recent_tasks']),20)
        self.assertTrue(context['older_tasks_omitted'])
        self.assertNotIn(final['run_id'],[r['id'] for r in context['recent_tasks']])
    def test_unlinked_conversation_has_no_project_context(self):
        row=self.store.create('独立问题','general')
        self.assertIsNone(self.store.run_project_context(row['run_id']))
        self.assertEqual(project_tools(self.store,row['run_id']),[])
    def test_idempotent_submission_does_not_change_snapshot_after_edit(self):
        first=self.store.create('任务','general',conversation_id=self.cid,request_id='test',fingerprint='fixed')
        self.store.update_project(self.project['id'],ProjectUpdate(name='新版本',revision=1))
        second=self.store.create('任务','general',conversation_id=self.cid,request_id='test',fingerprint='fixed')
        self.assertEqual(first['run_id'],second['run_id'])
        self.assertEqual(self.store.run_project_context(second['run_id'])['project']['revision'],1)
    def test_inspection_of_saved_project_result_does_not_claim_new_computation(self):
        from marine_agent.result_contract import execution_receipt
        row={'id':'read-only','status':'succeeded','mode':'general','operations':[], 'result':{'requested_action':'reply','agent':{'messages':[{'type':'tool','name':'inspect_project_task','content':json.dumps({'status':'succeeded','reported_answer':'旧结果'})}]}}}
        self.assertEqual(execution_receipt(row)['kind'],'dialogue')

    def test_context_budget_keeps_latest_question(self):
        row=self.store.create('old'*100,'general')
        self.store.append_message(row['conversation_id'],'assistant','a'*100)
        self.store.append_message(row['conversation_id'],'user','最新追问')
        messages=self.store.conversation_context(row['conversation_id'],max_chars=30)
        self.assertEqual(messages,[{'role':'user','content':'最新追问'}])
    def test_project_brief_is_in_model_context_but_does_not_execute(self):
        row=self.store.create('项目目标是什么','general',conversation_id=self.cid)
        contexts=[]
        class Agent:
            def stream(self,payload,**kwargs):
                contexts.append(payload)
                yield {'messages':[{'type':'ai','content':'筛选候选，实验验证活性。'}]}
        with patch.object(api,'store',self.store),patch.object(api,'plan_turn',return_value=TurnPlan(action='inspect',reason='项目讨论',context_scope='project',inspection='activity_evidence')),patch.object(api,'build_agent',return_value=Agent()):
            api._controlled_execute(row['run_id'],api._execute_run,'项目目标是什么','general')
            self.assertEqual(api.get_run_report(row['run_id'])['execution']['kind'],'dialogue')
        self.assertIn('user_confirmed_context_not_experimental_evidence',json.dumps(contexts))
        self.assertFalse(self.store.get(row['run_id'])['operations'])


class ReferenceLatencyTests(unittest.TestCase):
    def test_direct_grounded_answer_includes_followup_context_and_validates_ids(self):
        seen=[]
        def invoke(p):
            seen.append(p)
            return {'verdict':'supported','answer':'PAE 是相对位置的预测误差。','fact_ids':['pae.definition']}
        result=compose_explanation('后者是什么',[{'role':'user','content':'pLDDT 与 PAE'}],references(['pae']),invoke=invoke)
        self.assertEqual(result['review_method'],'reference_grounded_composition')
        self.assertEqual(seen[0]['context'][0]['content'],'pLDDT 与 PAE')
        self.assertEqual(result['original_text'],'')
        self.assertIsNone(result['fallback_reason'])
    def test_existing_text_json_is_parsed_without_second_provider_request(self):
        from unittest.mock import Mock
        model=Mock();model.with_structured_output.return_value.invoke.return_value={'parsed':None,'raw':AIMessage(content=json.dumps({'verdict':'supported','answer':'说明','fact_ids':['pae.definition']}))}
        result=_structured_result(model,[])
        self.assertEqual(result.answer,'说明')
        model.with_structured_output.return_value.invoke.assert_called_once()
    def test_direct_path_skips_agent_but_preserves_dialogue_receipt(self):
        with tempfile.TemporaryDirectory() as tmp,patch.object(api,'store',api.RunStore(Path(tmp)/'state.sqlite3')):
            row=api.store.create('解释PAE','general')
            response={'final_text':'参考回答','decision':'reference_composition_supported','fact_ids':['pae.definition'],'knowledge_version':'test'}
            with patch.object(api,'plan_turn',return_value=TurnPlan(action='reply',reason='定义',reference_topics=['pae'],reference_only=True)),patch.object(api,'build_agent') as graph,patch('marine_agent.scientific_knowledge.compose_explanation',return_value=response):
                api._controlled_execute(row['run_id'],api._execute_run,'解释PAE','general')
            graph.assert_not_called()
            self.assertEqual(api.get_run_report(row['run_id'])['execution']['kind'],'dialogue')
            self.assertFalse(api.store.get(row['run_id'])['operations'])


class ProjectApiTests(unittest.TestCase):
    def test_authenticated_api_ownership_revision_and_history(self):
        from fastapi import FastAPI
        from fastapi.testclient import TestClient
        from marine_agent import access_control
        from marine_agent.research_projects import register_project_routes
        with tempfile.TemporaryDirectory() as tmp:
            store=api.RunStore(Path(tmp)/'state.sqlite3');app=FastAPI()
            register_project_routes(app,lambda:store);app.add_middleware(access_control.AccountMiddleware)
            def identity(cookie):
                if cookie not in {'session=a','session=b'}:raise HTTPException(401,'login')
                return cookie[-1]
            with patch.object(access_control,'verify_session',side_effect=identity),TestClient(app) as client:
                self.assertEqual(client.get('/projects').status_code,401)
                a={'cookie':'session=a'};b={'cookie':'session=b'}
                response=client.post('/projects',headers=a,json={'name':'项目','objective':'实验验证'});self.assertEqual(response.status_code,201)
                pid=response.json()['id']
                self.assertEqual(client.get('/projects/'+pid,headers=b).status_code,404)
                self.assertEqual(client.put('/projects/'+pid,headers=b,json={'name':'覆盖','revision':1}).status_code,404)
                self.assertEqual(client.post('/projects/'+pid+'/conversations',headers=b).status_code,404)
                cid=client.post('/projects/'+pid+'/conversations',headers=a).json()['conversation_id']
                with account_scope('a'):
                    task=store.create('对话','general',conversation_id=cid)
                self.assertEqual(client.put('/conversations/'+cid+'/project',headers=a,json={'project_id':None,'expected_project_id':pid}).status_code,409)
                store.update(task['run_id'],status='succeeded')
                self.assertEqual(client.put('/projects/'+pid,headers=a,json={'name':'更新','revision':1}).status_code,200)
                self.assertEqual(client.put('/projects/'+pid,headers=a,json={'name':'旧表单','revision':1}).status_code,409)
                self.assertEqual(client.get('/projects/'+pid,headers=a).json()['conversations'][0]['id'],cid)

if __name__=='__main__':unittest.main()
