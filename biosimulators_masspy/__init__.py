from ._version import __version__
from .core import exec_sed_task, preprocess_sed_task, exec_sed_doc, exec_sedml_docs_in_combine_archive  # noqa: F401
import mass

__all__ = [
    '__version__',
    'get_simulator_version',
    'exec_sed_task',
    'preprocess_sed_task',
    'exec_sed_doc',
    'exec_sedml_docs_in_combine_archive',
]


def get_simulator_version():
    """ Get the version of MASSpy

    Returns:
        :obj:`str`: version
    """
    return mass.__version__
