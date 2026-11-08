"""Release tooling guards without touching production."""
import importlib.util,json,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def load(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'deploy/marine-agent'/f'{name}.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
class ReleaseToolTests(unittest.TestCase):
    def test_systemd_base_uses_latest_matching_override(self):
        tool=load('publish_backend')
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'a.conf').write_text('old-image')
            (root/'z.conf').write_text('current-image')
            self.assertEqual(tool.select_base_config('current-image',root),root/'z.conf')
            with self.assertRaises(RuntimeError):tool.select_base_config('missing',root)
    def test_release_scripts_compile(self):
        for name in ['publish_backend','publish_frontend']:
            tool=load(name);self.assertTrue(callable(tool.main))
    def test_verification_step_honors_supplied_environment(self):
        import importlib.util
        spec=importlib.util.spec_from_file_location('verify',ROOT/'verify.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as directory:
            log=Path(directory)/'env.log'
            step=module.run_step('env',['/bin/sh','-c','printf "%s" "$PATH"'],ROOT,log,{'PATH':'/custom-node-path'})
            self.assertEqual(step['exit_code'],0)
            self.assertEqual(log.read_text(),'/custom-node-path')

    def test_verification_runner_compiles(self):
        import py_compile
        py_compile.compile(str(ROOT/'verify.py'),doraise=True)
if __name__=='__main__':unittest.main()
