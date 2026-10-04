import importlib.util
import os
import unittest

_SCRIPT = os.path.join(os.path.dirname(__file__), os.pardir, '.github', 'scripts', 'next_version.py')
_spec = importlib.util.spec_from_file_location('next_version', _SCRIPT)
next_version = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(next_version)


class NextVersionTest(unittest.TestCase):

    def test_first_release_ignores_legacy_tags(self):
        self.assertEqual(next_version.next_version(['V12.0', '0.14.0'], 'major', 2026), (2026, 1, 0))
        self.assertEqual(next_version.next_version([], 'bugfix', 2026), (2026, 1, 0))

    def test_major_and_bugfix_within_a_year(self):
        tags = ['v2026.1.0', 'v2026.2.0', 'v2026.2.3', 'v2026.10.0']
        self.assertEqual(next_version.next_version(tags, 'major', 2026), (2026, 11, 0))
        self.assertEqual(next_version.next_version(tags, 'bugfix', 2026), (2026, 10, 1))

    def test_new_year_restarts_numbering(self):
        tags = ['v2026.4.2']
        self.assertEqual(next_version.next_version(tags, 'major', 2027), (2027, 1, 0))
        self.assertEqual(next_version.next_version(tags, 'bugfix', 2027), (2027, 1, 0))

    def test_rejects_release_dated_in_the_future(self):
        with self.assertRaises(ValueError):
            next_version.next_version(['v2027.1.0'], 'bugfix', 2026)
