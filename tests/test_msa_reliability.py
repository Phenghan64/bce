"""Exercise route recovery, bounded retries, audit retention and fragment metadata."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch
import requests
from tests.test_vina_conversation import api, ROOT
from marine_agent.openfold_tools import _search_rows


def load(name, file):
    spec=importlib.util.spec_from_file_location(name, ROOT/file)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

transport=load('msa_transport_test','openfold3/msa_transport.py')
pipeline=load('pipeline_msa_test','openfold3/predict_accession.py')


def response(code=200):
    r=Mock(status_code=code,headers={});return r


class MsaReliabilityTests(unittest.TestCase):
    def test_connection_reset_switches_route_and_next_request_keeps_working_route(self):
        original=Mock(side_effect=[requests.ConnectionError('reset'),response(),response()])
        network=transport.MsaTransport(original,proxy='socks5h://127.0.0.1:18199',sleep=Mock())
        for path in ['ticket/msa','ticket/123']:
            network.request(object(),'GET','https://api.colabfold.com/'+path)
        self.assertEqual([r['route'] for r in network.attempts],['direct','ssh_egress','ssh_egress'])
        self.assertEqual(original.call_args.kwargs['proxies']['https'],'socks5h://127.0.0.1:18199')
        self.assertNotIn('verify',original.call_args.kwargs)
        self.assertNotIn('127.0.0.1',json.dumps(network.report()))

    def test_total_transport_failure_is_bounded_and_retryable(self):
        original=Mock(side_effect=requests.ConnectionError('reset'))
        network=transport.MsaTransport(original,proxy='',sleep=Mock())
        with self.assertRaises(transport.MsaTransportFailure) as error:
            network.request(object(),'GET','https://api.colabfold.com/ticket/1')
        self.assertEqual(original.call_count,3)
        self.assertTrue(error.exception.detail['retryable'])

    def test_certificate_failure_is_not_bypassed_or_retried(self):
        original=Mock(side_effect=requests.exceptions.SSLError('bad cert'))
        network=transport.MsaTransport(original,proxy='',sleep=Mock())
        with self.assertRaises(transport.MsaTransportFailure) as error:
            network.request(object(),'POST','https://api.colabfold.com/ticket/msa')
        self.assertEqual(original.call_count,1)
        self.assertFalse(error.exception.detail['retryable'])

    def test_uncertain_post_timeout_is_not_automatically_resubmitted(self):
        original=Mock(side_effect=requests.exceptions.ReadTimeout())
        network=transport.MsaTransport(original,proxy='',sleep=Mock())
        with self.assertRaises(transport.MsaTransportFailure) as error:
            network.request(object(),'POST','https://api.colabfold.com/ticket/msa')
        self.assertEqual(original.call_count,1)
        self.assertEqual(error.exception.detail['error_code'],'msa_submission_uncertain')

    def test_rate_limit_get_retries_but_uncertain_post_502_does_not(self):
        original=Mock(side_effect=[response(429),response()])
        network=transport.MsaTransport(original,proxy='',sleep=Mock())
        self.assertEqual(network.request(object(),'GET','https://api.colabfold.com/ticket/1').status_code,200)
        original=Mock(return_value=response(502));network=transport.MsaTransport(original,proxy='',sleep=Mock())
        with self.assertRaises(transport.MsaTransportFailure):network.request(object(),'POST','https://api.colabfold.com/ticket/msa')
        self.assertEqual(original.call_count,1)

    def test_non_msa_request_is_untouched_and_overall_wait_is_bounded(self):
        original=Mock(return_value=response());clock=Mock(side_effect=[0,601])
        network=transport.MsaTransport(original,proxy='',sleep=Mock(),clock=clock)
        network.request(object(),'GET','https://rest.uniprot.org/a')
        self.assertNotIn('proxies',original.call_args.kwargs)
        with self.assertRaises(transport.MsaTransportFailure):network.request(object(),'GET','https://api.colabfold.com/ticket/1')
        self.assertEqual(original.call_count,1)

    def test_uniprot_description_fragment_flag_is_used_in_search_and_prediction_metadata(self):
        entry={'primaryAccession':'P46802','sequence':{'value':'A'*79,'length':79},'proteinDescription':{'flag':'Fragment'}}
        self.assertTrue(_search_rows({'results':[entry]})[0]['fragment'])
        metadata,seq=pipeline._uniprot_metadata(entry,'P46802','url',cache_hit=True,cache_path=None)
        self.assertTrue(metadata['fragment']);self.assertEqual(len(seq),79)

    def test_failed_msa_diagnostics_survive_temporary_directory_cleanup(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'results') as directory:
            root=Path(directory);query=root/'query.json';pipeline.write_query(query,'test','A'*79,'test')
            def fail(command,log,timeout):
                log.write_text('retained raw network diagnostics')
                out=log.parent/'cache';out.mkdir()
                (out/'failure.json').write_text(json.dumps({'error_code':'msa_network','reason':'连接失败','retryable':True}))
                (out/'network.json').write_text('{}')
                return {'ok':False,'status':'failed','output_tail':'raw'}
            with patch.object(pipeline,'ROOT',root),patch.object(pipeline,'run_command',side_effect=fail):
                result,_=pipeline.acquire_msa(query,'A'*79,'test',root/'output/msa',cache_dir=root/'cache',refresh=False,query_name='test')
            self.assertEqual(result['error_code'],'msa_network')
            self.assertEqual((root/result['log']).read_text(),'retained raw network diagnostics')
            self.assertTrue((root/'output/msa/failure.json').exists())
            self.assertFalse(list(root.glob('msa-*')))

    def test_frontend_error_is_concise_but_raw_failure_remains_in_task_record(self):
        with tempfile.TemporaryDirectory() as directory:
            store=api.RunStore(Path(directory)/'state.sqlite3')
            created=store.create('预测','manual_openfold',manual={'accession':'P46802'})
            payload={'status':'failed','stage':'msa','msa':{'reason':'服务连接失败','output_tail':'internal traceback'},'input':{'fragment':True}}
            with patch.object(api,'store',store),patch.object(api,'run_openfold_accession',Mock(invoke=Mock(return_value=json.dumps(payload)))):
                api._controlled_execute(created['run_id'],api._execute_openfold_spec,{'workflow':'openfold_accession','input':{'accession':'P46802'}},0)
            record=store.get(created['run_id'])
            self.assertEqual(record['status'],'failed')
            self.assertIn('尚未开始结构推理',record['error'])
            self.assertNotIn('traceback',record['error'])
            self.assertEqual(record['result']['agent']['openfold_failure']['msa']['output_tail'],'internal traceback')

if __name__=='__main__':unittest.main()
