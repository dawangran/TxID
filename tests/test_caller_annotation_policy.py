"""The caller wrapper must preserve explicit complete versus inferred annotation modes."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


class CallerAnnotationPolicyTests(unittest.TestCase):
    def test_isoquant_meta_feature_inference_preserves_supplied_models(self):
        scripts = Path(__file__).resolve().parents[1] / 'workflows/publication_benchmark/scripts'
        with patch.object(sys, 'path', [str(scripts), *sys.path]):
            spec = importlib.util.spec_from_file_location('annotation_policy_caller', scripts / 'run_caller.py')
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
        for infer in (False, True):
            with self.subTest(infer=infer), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                supplied = 'chr22\tcaller\texon\t100\t200\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
                for name, content in [('ref.fa', '>chr22\n' + 'A' * 300), ('annotation.gtf', supplied), ('in.bam', 'alignment stub')]:
                    (root / name).write_text(content)
                work = root / 'work'
                commands = []

                def fake_run(argv, **kwargs):
                    commands.append(argv)
                    (work / 'case.transcript_models.gtf').write_text(supplied)
                    return {'command': argv, 'returncode': 0}

                argv = ['run_caller', '--caller', 'isoquant', '--dataset', 'case',
                        '--reference', str(root/'ref.fa'), '--annotation', str(root/'annotation.gtf'),
                        '--bam', str(root/'in.bam'), '--read-type', 'pacbio_ccs',
                        '--output-gtf', str(root/'output.gtf'), '--work-dir', str(work),
                        '--run-json', str(root/'run.json')]
                if infer:
                    argv.append('--infer-annotation-meta-features')
                with patch.object(sys, 'argv', argv), patch.object(module, '_version', return_value='3.13.0'), patch.object(module, 'run_command', side_effect=fake_run):
                    self.assertEqual(module.main(), 0)
                self.assertEqual('--complete_genedb' in commands[0], not infer)
                self.assertEqual(commands[0][commands[0].index('--genedb') + 1], str(root/'annotation.gtf'))
                self.assertEqual((root/'output.gtf').read_text(), supplied)
                self.assertEqual((root/'annotation.gtf').read_text(), supplied)


if __name__ == '__main__':
    unittest.main()
