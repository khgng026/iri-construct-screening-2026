"""Portable importer and host-monitor regression checks; fictional bytes only."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import reproduce
import resource_guard

class ReproductionTests(unittest.TestCase):
    def fixture(self, root):
        config = root / 'config/reproduction'; config.mkdir(parents=True)
        archive = root / 'input.zip'
        content = b'fictional numeric fixture'
        with zipfile.ZipFile(archive, 'w') as z:
            z.writestr('data/raw/a.txt', content)
            z.writestr('DATA_NOTICE.txt', 'fictional')
        (config / 'data_bundle_manifest.json').write_text(json.dumps({'data/raw/a.txt': hashlib.sha256(content).hexdigest()}))
        (config / 'bundle.json').write_text(json.dumps({'sha256': hashlib.sha256(archive.read_bytes()).hexdigest()}))
        return archive

    def test_import_validates_and_is_idempotent(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive = self.fixture(root)
            with patch.object(reproduce, 'ROOT', root), patch.object(reproduce, 'verify_inputs') as verify:
                reproduce.install_bundle(archive); reproduce.install_bundle(archive)
                self.assertEqual(verify.call_count, 2)
            self.assertEqual((root / 'data/raw/a.txt').read_bytes(), b'fictional numeric fixture')

    def test_corruption_and_conflicting_input_fail_before_writes(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive = self.fixture(root)
            with patch.object(reproduce, 'ROOT', root), patch.object(reproduce, 'verify_inputs'):
                original = archive.read_bytes(); archive.write_bytes(original + b'changed')
                with self.assertRaisesRegex(ValueError, 'checksum'): reproduce.install_bundle(archive)
                self.assertFalse((root / 'data').exists())
                archive.write_bytes(original)
                target = root / 'data/raw/a.txt'; target.parent.mkdir(parents=True); target.write_text('existing')
                with self.assertRaisesRegex(ValueError, 'Refusing'): reproduce.install_bundle(archive)
                self.assertEqual(target.read_text(), 'existing')

    def test_linux_memory_and_power_fail_closed(self):
        fake = SimpleNamespace(virtual_memory=lambda: SimpleNamespace(available=2**30, total=16*2**30),
                               swap_memory=lambda: SimpleNamespace(used=0),
                               sensors_battery=lambda: SimpleNamespace(power_plugged=False))
        with patch.object(resource_guard.sys, 'platform', 'linux'), patch.dict(sys.modules, {'psutil': fake}):
            sample = resource_guard.snapshot()
        self.assertIn('memory pressure is not normal', resource_guard.reasons(sample))
        self.assertIn('AC power disconnected', resource_guard.reasons(sample))
