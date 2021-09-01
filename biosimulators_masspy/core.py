""" BioSimulators-compliant command-line interface to the `MASSpy <https://masspy.readthedocs.io/>`_ simulation program.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2021-08-12
:Copyright: 2021, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from .data_model import KISAO_ALGORITHM_MAP
from biosimulators_utils.combine.exec import exec_sedml_docs_in_archive
from biosimulators_utils.config import get_config, Config  # noqa: F401
from biosimulators_utils.log.data_model import CombineArchiveLog, TaskLog  # noqa: F401
from biosimulators_utils.report.data_model import ReportFormat, VariableResults, SedDocumentResults  # noqa: F401
from biosimulators_utils.sedml import validation
from biosimulators_utils.sedml.data_model import (Task, ModelLanguage, ModelAttributeChange,  # noqa: F401
                                                  SteadyStateSimulation, UniformTimeCourseSimulation,
                                                  Variable, Symbol)
from biosimulators_utils.sedml.exec import exec_sed_doc
from biosimulators_utils.simulator.utils import get_algorithm_substitution_policy
from biosimulators_utils.utils.core import raise_errors_warnings, validate_str_value, parse_value
from biosimulators_utils.viz.data_model import VizFormat  # noqa: F401
from biosimulators_utils.warnings import warn, BioSimulatorsWarning
from kisao.data_model import AlgorithmSubstitutionPolicy, ALGORITHM_SUBSTITUTION_POLICY_LEVELS
from kisao.utils import get_preferred_substitute_algorithm_by_ids
import functools
import numpy
import mass
import mass.io.sbml

__all__ = ['exec_sedml_docs_in_combine_archive', 'exec_sed_task']


def exec_sedml_docs_in_combine_archive(archive_filename, out_dir, config=None):
    """ Execute the SED tasks defined in a COMBINE/OMEX archive and save the outputs

    Args:
        archive_filename (:obj:`str`): path to COMBINE/OMEX archive
        out_dir (:obj:`str`): path to store the outputs of the archive

            * CSV: directory in which to save outputs to files
              ``{ out_dir }/{ relative-path-to-SED-ML-file-within-archive }/{ report.id }.csv``
            * HDF5: directory in which to save a single HDF5 file (``{ out_dir }/reports.h5``),
              with reports at keys ``{ relative-path-to-SED-ML-file-within-archive }/{ report.id }`` within the HDF5 file

        config (:obj:`Config`, optional): BioSimulators common configuration

    Returns:
        :obj:`tuple`:

            * :obj:`SedDocumentResults`: results
            * :obj:`CombineArchiveLog`: log
    """
    sed_doc_executer = functools.partial(exec_sed_doc, exec_sed_task)
    return exec_sedml_docs_in_archive(sed_doc_executer, archive_filename, out_dir,
                                      apply_xml_model_changes=True,
                                      config=config)


def exec_sed_task(task, variables, log=None, config=None):
    """ Execute a task and save its results

    Args:
       task (:obj:`Task`): task
       variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
       log (:obj:`TaskLog`, optional): log for the task
       config (:obj:`Config`, optional): BioSimulators common configuration

    Returns:
        :obj:`tuple`:

            :obj:`VariableResults`: results of variables
            :obj:`TaskLog`: log

    Raises:
        :obj:`NotImplementedError`:

          * Task requires a time course that MASSpy doesn't support
          * Task requires an algorithm that MASSpy doesn't support
    """
    config = config or get_config()

    if config.LOG and not log:
        log = TaskLog()

    # validate task
    model = task.model
    sim = task.simulation

    if config.VALIDATE_SEDML:
        raise_errors_warnings(validation.validate_task(task),
                              error_summary='Task `{}` is invalid.'.format(task.id))
        raise_errors_warnings(validation.validate_model_language(model.language, ModelLanguage.SBML),
                              error_summary='Language for model `{}` is not supported.'.format(model.id))
        raise_errors_warnings(validation.validate_model_change_types(model.changes, ()),
                              error_summary='Changes for model `{}` are not supported.'.format(model.id))
        raise_errors_warnings(*validation.validate_model_changes(model),
                              error_summary='Changes for model `{}` are invalid.'.format(model.id))
        raise_errors_warnings(validation.validate_simulation_type(sim, (UniformTimeCourseSimulation, )),
                              error_summary='{} `{}` is not supported.'.format(sim.__class__.__name__, sim.id))
        raise_errors_warnings(*validation.validate_simulation(sim),
                              error_summary='Simulation `{}` is invalid.'.format(sim.id))
        raise_errors_warnings(*validation.validate_data_generator_variables(variables),
                              error_summary='Data generator variables for task `{}` are invalid.'.format(task.id))

    target_x_paths_to_sbml_ids = validation.validate_variable_xpaths(variables, model.source, attr='id')

    if config.VALIDATE_SEDML_MODELS:
        raise_errors_warnings(*validation.validate_model(model, [], working_dir='.'),
                              error_summary='Model `{}` is invalid.'.format(model.id),
                              warning_summary='Model `{}` may be invalid.'.format(model.id))

    # read model
    mass_model = mass.io.sbml.read_sbml_model(model.source)

    # validate variables
    met_ids = ['M_' + mass_met.id for mass_met in mass_model.metabolites]
    rxn_ids = ['R_' + mass_rxn.id for mass_rxn in mass_model.reactions]

    invalid_symbols = []
    invalid_targets = []
    for variable in variables:
        if variable.symbol:
            if variable.symbol != Symbol.time.value:
                invalid_symbols.append(variable.symbol)

        else:
            sbml_id = target_x_paths_to_sbml_ids.get(variable.target, None)

            if not sbml_id or not (
                (sbml_id.startswith('M_') and sbml_id in met_ids) or
                (sbml_id.startswith('R_') and sbml_id in rxn_ids)
            ):
                invalid_targets.append(variable.target)

    if invalid_symbols:
        msg = (
            'The following symbols are not supported:\n  - {}'
            '\n'
            '\n'
            'Only following symbols are supported:\n  - {}'
        ).format(
            '\n  - '.join(sorted(invalid_symbols)),
            '\n  - '.join(sorted([Symbol.time.value])),
        )
        raise NotImplementedError(msg)

    if invalid_targets:
        valid_targets = []
        for mass_met in mass_model.metabolites:
            valid_targets.append("/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_{}']".format(mass_met.id))
        for mass_rxn in mass_model.reactions:
            valid_targets.append("/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_{}']".format(mass_rxn.id))

        msg = (
            'The following targets are not supported:\n  - {}'
            '\n'
            '\n'
            'Only following targets are supported:\n  - {}'
        ).format(
            '\n  - '.join(sorted(invalid_targets)),
            '\n  - '.join(sorted(valid_targets)),
        )
        raise ValueError(msg)

    # instantiate simulation
    mass_sim = mass.Simulation(reference_model=mass_model, verbose=config.VERBOSE)

    # configure simulation time course
    number_of_points = sim.number_of_points * (
        sim.output_end_time - sim.initial_time
    ) / (
        sim.output_end_time - sim.output_start_time
    ) + 1
    if abs(number_of_points % 1) > 1e-8:
        msg = (
            'The number of simulation points `{}` must be an integer:'
            '\n  Initial time: {}'
            '\n  Output start time: {}'
            '\n  Output end time: {}'
            '\n  Number of points: {}'
        ).format(number_of_points, sim.initial_time, sim.output_start_time, sim.output_start_time, sim.number_of_points)
        raise NotImplementedError(msg)

    number_of_points = round(number_of_points)

    time = (sim.initial_time, sim.output_end_time, number_of_points)

    # configure simulation algorithm
    algorithm_substitution_policy = get_algorithm_substitution_policy(config=config)
    exec_kisao_id = get_preferred_substitute_algorithm_by_ids(
        sim.algorithm.kisao_id, KISAO_ALGORITHM_MAP.keys(),
        substitution_policy=algorithm_substitution_policy)

    alg_props = KISAO_ALGORITHM_MAP[exec_kisao_id]
    mass_sim.roadrunner.setIntegrator(alg_props['id'])

    # configure parameters of the simulation algorithm
    if exec_kisao_id == sim.algorithm.kisao_id:
        for change in sim.algorithm.changes:
            param_props = alg_props['parameters'].get(change.kisao_id, None)

            if param_props:

                if validate_str_value(change.new_value, param_props['type']):
                    value = parse_value(change.new_value, param_props['type'])
                    setattr(mass_sim.integrator, param_props['id'], value)

                else:
                    if (
                        ALGORITHM_SUBSTITUTION_POLICY_LEVELS[algorithm_substitution_policy]
                        <= ALGORITHM_SUBSTITUTION_POLICY_LEVELS[AlgorithmSubstitutionPolicy.NONE]
                    ):
                        msg = "'{}' is not a valid {} value for parameter {}".format(
                            change.new_value, param_props['type'].name, change.kisao_id)
                        raise ValueError(msg)
                    else:
                        msg = "'{}' was ignored because it is not a valid {} value for parameter {}".format(
                            change.new_value, param_props['type'].name, change.kisao_id)
                        warn(msg, BioSimulatorsWarning)

            else:
                if (
                    ALGORITHM_SUBSTITUTION_POLICY_LEVELS[algorithm_substitution_policy]
                    <= ALGORITHM_SUBSTITUTION_POLICY_LEVELS[AlgorithmSubstitutionPolicy.NONE]
                ):
                    msg = "".join([
                        "Algorithm parameter with KiSAO id '{}' is not supported. ".format(change.kisao_id),
                        "Parameter must have one of the following KiSAO ids:\n  - {}".format('\n  - '.join(
                            '{}: {} ({})'.format(kisao_id, param_props['id'], param_props['name'])
                            for kisao_id, param_props in alg_props['parameters'].items())),
                    ])
                    raise NotImplementedError(msg)
                else:
                    msg = "".join([
                        "Algorithm parameter with KiSAO id '{}' was ignored because it is not supported. ".format(change.kisao_id),
                        "Parameter must have one of the following KiSAO ids:\n  - {}".format('\n  - '.join(
                            '{}: {} ({})'.format(kisao_id, param_props['id'], param_props['name'])
                            for kisao_id, param_props in alg_props['parameters'].items())),
                    ])
                    warn(msg, BioSimulatorsWarning)

    if hasattr(mass_sim.integrator, 'variable_step_size'):
        mass_sim.integrator.variable_step_size = False

    # execute simulation
    met_concs, rxn_fluxes = mass_sim.simulate(mass_model, time=time)

    # check simulation executed successfully
    if numpy.any(numpy.isnan(met_concs.to_frame())):
        msg = 'Simulation failed with algorithm `{}` ({})'.format(exec_kisao_id, alg_props['id'])
        for i_param in range(mass_sim.integrator.getNumParams()):
            param_name = mass_sim.integrator.getParamName(i_param)
            msg += '\n  - {}: {}'.format(param_name, getattr(mass_sim.integrator, param_name))
        raise ValueError(msg)

    # transform simulation results
    met_ics = {met.id: met.initial_condition for met in mass_model.metabolites}

    variable_results = VariableResults()
    for variable in variables:
        if variable.symbol:
            variable_results[variable.id] = met_concs._time[-(sim.number_of_points + 1):]

        else:
            sbml_id = target_x_paths_to_sbml_ids[variable.target]

            if sbml_id.startswith('M_'):
                if sbml_id[2:] in met_concs:
                    variable_results[variable.id] = met_concs[sbml_id[2:]][-(sim.number_of_points + 1):]
                else:
                    variable_results[variable.id] = numpy.full((sim.number_of_points + 1,), met_ics[sbml_id[2:]])

            else:
                variable_results[variable.id] = rxn_fluxes[sbml_id[2:]][-(sim.number_of_points + 1):]

    # log action
    if config.LOG:
        log.algorithm = exec_kisao_id

        log.simulator_details = {}
        log.simulator_details['integrator'] = mass_sim.integrator.getName()
        for i_param in range(mass_sim.integrator.getNumParams()):
            param_name = mass_sim.integrator.getParamName(i_param)
            log.simulator_details[param_name] = getattr(mass_sim.integrator, param_name)

    ############################
    # return the result of each variable and log
    return variable_results, log
