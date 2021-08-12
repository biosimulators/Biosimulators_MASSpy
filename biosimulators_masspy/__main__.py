""" BioSimulators-compliant command-line interface to the
`MASSpy <https://masspy.readthedocs.io/>`_ simulation program.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2021-08-12
:Copyright: 2021, BioSimulators Team
:License: MIT
"""

from ._version import __version__
from .core import exec_sedml_docs_in_combine_archive
from biosimulators_utils.simulator.cli import build_cli
import mass

App = build_cli('biosimulators-masspy', __version__,
                'MASSpy', mass.__version__, 'https://masspy.readthedocs.io/',
                exec_sedml_docs_in_combine_archive)


def main():
    with App() as app:
        app.run()
