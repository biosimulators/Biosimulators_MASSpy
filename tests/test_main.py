from biosimulators_masspy import __main__
from unittest import mock
import unittest


class CliTestCase(unittest.TestCase):
    def test_raw(self):
        with mock.patch('sys.argv', ['', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegex(context.Exception, 'usage: ')
