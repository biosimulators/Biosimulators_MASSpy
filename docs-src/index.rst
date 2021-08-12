BioSimulators-MASSpy documentation
=====================================

BioSimulators-MASSpy provides a `BioSimulators <https://biosimulators.org>`_-compliant command-line interface to the `MASSpy <https://masspy.readthedocs.io/>`_ tool for kinetic simulation of metabolic reaction networks. A Docker image for this package is also available.

This command-line interface and Docker image enable users to use MASSpy to execute `COMBINE/OMEX archives <https://combinearchive.org/>`_ that describe one or more simulation experiments (in `SED-ML format <https://sed-ml.org>`_) of one or more kinetic models of metabolic networks in the `MASSpy schema <https://masspy.readthedocs.io/en/stable/tutorials/reading_writing_models.html>`_ for the SBML format.

A list of the algorithms and algorithm parameters supported by MASSpy is available at `BioSimulators <https://biosimulators.org/simulators/masspy>`_.

A simple web application and web service for using MASSpy to execute COMBINE/OMEX archives is also available at `runBioSimulations <https://run.biosimulations.org>`_.

Contents
--------

.. toctree::
   :maxdepth: 2

   installation.rst
   tutorial.rst
   API documentation <source/biosimulators_masspy.rst>
   about.rst
   genindex.rst
