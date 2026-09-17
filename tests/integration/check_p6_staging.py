import importlib.util
import io
from pathlib import Path
import unittest
ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('stage', ROOT / 'tools/stage_gse80336_local.py')
stage = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stage)
class Records(unittest.TestCase):
    def test_prefix_stops_after_complete_record(self):
        out = io.BytesIO()
        stage.copy_fastq_prefix(io.BytesIO(b'@a\nAC\n+\n!!\n@b\nGG\n+\n!!\n'), out, 1, lambda: None)
        self.assertEqual(out.getvalue(), b'@a\nAC\n+\n!!\n')
    def test_rejects_lengths_and_early_eof(self):
        for data in (b'@a\nAC\n+\n!\n', b'@a\nAC\n+\n', b'bad\nAC\n+\n!!\n'):
            with self.assertRaises(ValueError):
                stage.copy_fastq_prefix(io.BytesIO(data), io.BytesIO(), 1, lambda: None)
unittest.main()
