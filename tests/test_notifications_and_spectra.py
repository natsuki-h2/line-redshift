from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from line_redshift.cli import desktop_notify, load_spectra

class NotificationAndSpectraTests(unittest.TestCase):
    def test_notification_passes_message_as_argument_not_shell_code(self):
        message='A quoted "candidate" $(not-a-command)'
        with patch('line_redshift.cli.sys.platform','darwin'), patch('line_redshift.cli.shutil.which',return_value='/usr/bin/osascript'), patch('line_redshift.cli.subprocess.run') as run:
            desktop_notify(message)
            args,kwargs=run.call_args
            self.assertEqual(args[0][-2],message)
            self.assertNotIn('shell',kwargs)
    def test_descending_spectrum_is_sorted_and_units_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'s.dat').write_text('102 3\n101 2\n100 1\n')
            (p/'manifest.csv').write_text('pair_id,path,flux_unit\nx,s.dat,mJy/beam\n')
            r=load_spectra(p/'manifest.csv')['x'][0]
            self.assertEqual(list(r[0]),[100,101,102]);self.assertEqual(list(r[1]),[1,2,3])
    def test_mixed_flux_units_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'s.dat').write_text('100 1\n101 2\n')
            (p/'manifest.csv').write_text('pair_id,path,flux_unit\nx,s.dat,Jy\nx,s.dat,mJy\n')
            with self.assertRaisesRegex(ValueError,'common flux unit'):load_spectra(p/'manifest.csv')
