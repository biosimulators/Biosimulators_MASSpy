""" BioSimulators-compliant command-line interface to the `MASSpy <https://masspy.readthedocs.io/>`_ simulation program.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2021-08-12
:Copyright: 2021, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from .data_model import KISAO_ALGORITHM_MAP
from biosimulators_utils.combine.exec import exec_sedml_docs_in_archive
from biosimulators_utils.config import get_config, Config  # noqa: F401
from biosimulators_utils.log.data_model import CombineArchiveLog, TaskLog, StandardOutputErrorCapturerLevel  # noqa: F401
from biosimulators_utils.report.data_model import ReportFormat, VariableResults, SedDocumentResults  # noqa: F401
from biosimulators_utils.sedml import validation
from biosimulators_utils.sedml.data_model import (Task, ModelLanguage, ModelAttributeChange,  # noqa: F401
                                                  SteadyStateSimulation, UniformTimeCourseSimulation,
                                                  Variable, Symbol)
from biosimulators_utils.sedml.exec import exec_sed_doc as base_exec_sed_doc
from biosimulators_utils.simulator.utils import get_algorithm_substitution_policy
from biosimulators_utils.utils.core import raise_errors_warnings, validate_str_value, parse_value
from biosimulators_utils.viz.data_model import VizFormat  # noqa: F401
from biosimulators_utils.warnings import warn, BioSimulatorsWarning
from kisao.data_model import AlgorithmSubstitutionPolicy, ALGORITHM_SUBSTITUTION_POLICY_LEVELS
from kisao.utils import get_preferred_substitute_algorithm_by_ids
import lxml.etree
import mass
import mass.io.sbml
import numpy
import libsbml

__all__ = ['exec_sedml_docs_in_combine_archive', 'exec_sed_doc', 'exec_sed_task', 'preprocess_sed_task']


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
    return exec_sedml_docs_in_archive(exec_sed_doc, archive_filename, out_dir,
                                      apply_xml_model_changes=True,
                                      config=config)


def exec_sed_doc(doc, working_dir, base_out_path, rel_out_path=None,
                 apply_xml_model_changes=True,
                 log=None, indent=0, pretty_print_modified_xml_models=False,
                 log_level=StandardOutputErrorCapturerLevel.c, config=None):
    """ Execute the tasks specified in a SED document and generate the specified outputs

    Args:
        doc (:obj:`SedDocument` or :obj:`str`): SED document or a path to SED-ML file which defines a SED document
        working_dir (:obj:`str`): working directory of the SED document (path relative to which models are located)

        base_out_path (:obj:`str`): path to store the outputs

            * CSV: directory in which to save outputs to files
              ``{base_out_path}/{rel_out_path}/{report.id}.csv``
            * HDF5: directory in which to save a single HDF5 file (``{base_out_path}/reports.h5``),
              with reports at keys ``{rel_out_path}/{report.id}`` within the HDF5 file

        rel_out_path (:obj:`str`, optional): path relative to :obj:`base_out_path` to store the outputs
        apply_xml_model_changes (:obj:`bool`, optional): if :obj:`True`, apply any model changes specified in the SED-ML file before
            calling :obj:`task_executer`.
        log (:obj:`SedDocumentLog`, optional): log of the document
        indent (:obj:`int`, optional): degree to indent status messages
        pretty_print_modified_xml_models (:obj:`bool`, optional): if :obj:`True`, pretty print modified XML models
        log_level (:obj:`StandardOutputErrorCapturerLevel`, optional): level at which to log output
        config (:obj:`Config`, optional): BioSimulators common configuration
        simulator_config (:obj:`SimulatorConfig`, optional): tellurium configuration

    Returns:
        :obj:`tuple`:

            * :obj:`ReportResults`: results of each report
            * :obj:`SedDocumentLog`: log of the document
    """
    return base_exec_sed_doc(exec_sed_task, doc, working_dir, base_out_path,
                             rel_out_path=rel_out_path,
                             apply_xml_model_changes=apply_xml_model_changes,
                             log=log,
                             indent=indent,
                             pretty_print_modified_xml_models=pretty_print_modified_xml_models,
                             log_level=log_level,
                             config=config)


def exec_sed_task(task, variables, preprocessed_task=None, log=None, config=None):
    """ Execute a task and save its results

    Args:
        task (:obj:`Task`): task
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
        preprocessed_task (:obj:`dict`, optional): preprocessed information about the task, including possible
            model changes and variables. This can be used to avoid repeatedly executing the same initialization
            for repeated calls to this method.
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

    if preprocessed_task is None:
        preprocessed_task = preprocess_sed_task(task, variables, config=config)

    # validate task
    model = task.model
    sim = task.simulation

    # get model
    mass_model = preprocessed_task['model']['model']

    # modify model
    if model.changes:
        raise_errors_warnings(validation.validate_model_change_types(model.changes, (ModelAttributeChange, )),
                              error_summary='Changes for model `{}` are not supported.'.format(model.id))
        model_change_target_mass_map = preprocessed_task['model']['model_change_target_mass_map']
        for change in model.changes:
            component, attr_name = model_change_target_mass_map[change.target]
            new_value = float(change.new_value)
            if isinstance(component, dict):
                component[attr_name] = new_value
            else:
                setattr(component, attr_name, new_value)

    # get simulation
    mass_sim = preprocessed_task['simulation']['simulation']

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

    # execute simulation
    met_concs, rxn_fluxes = mass_sim.simulate(mass_model, time=time)

    # check simulation executed successfully
    if numpy.any(numpy.isnan(met_concs.to_frame())):
        msg = 'Simulation failed with algorithm `{}` ({})'.format(
            preprocessed_task['simulation']['algorithm_kisao_id'],
            preprocessed_task['simulation']['algorithm_roadrunner_id'])
        for i_param in range(mass_sim.integrator.getNumParams()):
            param_name = mass_sim.integrator.getParamName(i_param)
            msg += '\n  - {}: {}'.format(param_name, getattr(mass_sim.integrator, param_name))
        raise ValueError(msg)

    # transform simulation results
    met_ics = {met.id: met.initial_condition for met in mass_model.metabolites}
    variable_target_sbml_id_map = preprocessed_task['model']['variable_target_sbml_id_map']
    variable_results = VariableResults()
    for variable in variables:
        if variable.symbol:
            variable_results[variable.id] = met_concs._time[-(sim.number_of_points + 1):]

        else:
            sbml_id = variable_target_sbml_id_map[variable.target]

            if sbml_id.startswith('M_'):
                if sbml_id[2:] in met_concs:
                    variable_results[variable.id] = met_concs[sbml_id[2:]][-(sim.number_of_points + 1):]
                else:
                    variable_results[variable.id] = numpy.full((sim.number_of_points + 1,), met_ics[sbml_id[2:]])

            else:
                variable_results[variable.id] = rxn_fluxes[sbml_id[2:]][-(sim.number_of_points + 1):]

    # log action
    if config.LOG:
        log.algorithm = preprocessed_task['simulation']['algorithm_kisao_id']

        log.simulator_details = {}
        log.simulator_details['integrator'] = mass_sim.integrator.getName()
        for i_param in range(mass_sim.integrator.getNumParams()):
            param_name = mass_sim.integrator.getParamName(i_param)
            log.simulator_details[param_name] = getattr(mass_sim.integrator, param_name)

    ############################
    # return the result of each variable and log
    return variable_results, log


def preprocess_sed_task(task, variables, config=None):
    """ Preprocess a SED task, including its possible model changes and variables. This is useful for avoiding
    repeatedly initializing tasks on repeated calls of :obj:`exec_sed_task`.

    Args:
        task (:obj:`Task`): task
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
        config (:obj:`Config`, optional): BioSimulators common configuration

    Returns:
        :obj:`dict`: preprocessed information about the task
    """
    config = config or get_config()

    # validate task
    model = task.model
    sim = task.simulation

    if config.VALIDATE_SEDML:
        raise_errors_warnings(validation.validate_task(task),
                              error_summary='Task `{}` is invalid.'.format(task.id))
        raise_errors_warnings(validation.validate_model_language(model.language, ModelLanguage.SBML),
                              error_summary='Language for model `{}` is not supported.'.format(model.id))
        raise_errors_warnings(validation.validate_model_change_types(model.changes, (ModelAttributeChange, )),
                              error_summary='Changes for model `{}` are not supported.'.format(model.id))
        raise_errors_warnings(*validation.validate_model_changes(model),
                              error_summary='Changes for model `{}` are invalid.'.format(model.id))
        raise_errors_warnings(validation.validate_simulation_type(sim, (UniformTimeCourseSimulation, )),
                              error_summary='{} `{}` is not supported.'.format(sim.__class__.__name__, sim.id))
        raise_errors_warnings(*validation.validate_simulation(sim),
                              error_summary='Simulation `{}` is invalid.'.format(sim.id))
        raise_errors_warnings(*validation.validate_data_generator_variables(variables),
                              error_summary='Data generator variables for task `{}` are invalid.'.format(task.id))

    model_etree = lxml.etree.parse(model.source)

    model_change_target_sbml_id_map = validation.validate_target_xpaths(model.changes, model_etree, attr='id')
    variable_target_sbml_id_map = validation.validate_target_xpaths(variables, model_etree, attr='id')

    if config.VALIDATE_SEDML_MODELS:
        raise_errors_warnings(*validation.validate_model(model, [], working_dir='.'),
                              error_summary='Model `{}` is invalid.'.format(model.id),
                              warning_summary='Model `{}` may be invalid.'.format(model.id))

    # read model
    # Have to convert to L3v1 by hand; automatic conversion fails.
    doc = libsbml.readSBML(model.source)
    doc.setLevelAndVersion(3, 1)
    sbml = libsbml.writeSBMLToString(doc)
    mass_model = mass.io.sbml.read_sbml_model(sbml)
    met_ids = ['M_' + mass_met.id for mass_met in mass_model.metabolites]
    rxn_ids = ['R_' + mass_rxn.id for mass_rxn in mass_model.reactions]

    # validate model changes
    model_change_target_mass_map = {}

    sbml_id_mass_parameter_map = {}
    for reaction in mass_model.reactions:
        if reaction.equilibrium_constant is not None:
            sbml_id_mass_parameter_map['Keq_R_' + reaction.id] = (reaction, 'equilibrium_constant')
            sbml_id_mass_parameter_map['Keq_' + reaction.id] = (reaction, 'equilibrium_constant')

        if reaction.forward_rate_constant is not None:
            sbml_id_mass_parameter_map['kf_R_' + reaction.id] = (reaction, 'forward_rate_constant')
            sbml_id_mass_parameter_map['kf_' + reaction.id] = (reaction, 'forward_rate_constant')

        if reaction.reverse_rate_constant is not None:
            sbml_id_mass_parameter_map['kr_R_' + reaction.id] = (reaction, 'reverse_rate_constant')
            sbml_id_mass_parameter_map['kr_' + reaction.id] = (reaction, 'reverse_rate_constant')

        if reaction.steady_state_flux is not None:
            sbml_id_mass_parameter_map['v_R_' + reaction.id] = (reaction, 'steady_state_flux')
            sbml_id_mass_parameter_map['v_' + reaction.id] = (reaction, 'steady_state_flux')

    for sbml_id in mass_model.custom_parameters.keys():
        sbml_id_mass_parameter_map[sbml_id] = (mass_model.custom_parameters, sbml_id)

    for sbml_id in mass_model.boundary_conditions.keys():
        sbml_id_mass_parameter_map[sbml_id] = (mass_model.boundary_conditions, sbml_id)

    invalid_changes = []
    for target, sbml_id in model_change_target_sbml_id_map.items():
        if sbml_id in met_ids:
            model_change_target_mass_map[target] = (mass_model.metabolites[met_ids.index(sbml_id)], 'ic')

        elif sbml_id in sbml_id_mass_parameter_map:
            model_change_target_mass_map[target] = sbml_id_mass_parameter_map[sbml_id]

        else:
            invalid_changes.append(target)

    if invalid_changes:
        msg = (
            'The following model change targets are not supported:\n  - {}'
            '\n'
            '\n'
            'Only following change targets are supported:\n  - {}'
        ).format(
            '\n  - '.join(sorted(invalid_changes)),
            '\n  - '.join(sorted(met_ids + list(sbml_id_mass_parameter_map.keys()))),
        )
        raise ValueError(msg)

    # validate variables
    invalid_symbols = []
    invalid_targets = []
    for variable in variables:
        if variable.symbol:
            if variable.symbol != Symbol.time.value:
                invalid_symbols.append(variable.symbol)

        else:
            sbml_id = variable_target_sbml_id_map.get(variable.target, None)

            if not sbml_id or not (sbml_id in met_ids or sbml_id in rxn_ids):
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

    # implemented this way due to tellurium bug: https://github.com/sys-bio/tellurium/issues/546
    if 'variable_step_size' in dir(mass_sim.integrator):
        mass_sim.integrator.variable_step_size = False

    ############################
    # return preprocessed information
    return {
        'model': {
            'model': mass_model,
            'model_change_target_mass_map': model_change_target_mass_map,
            'variable_target_sbml_id_map': variable_target_sbml_id_map,
        },
        'simulation': {
            'simulation': mass_sim,
            'algorithm_kisao_id': exec_kisao_id,
            'algorithm_roadrunner_id': alg_props['id'],
        },
    }
