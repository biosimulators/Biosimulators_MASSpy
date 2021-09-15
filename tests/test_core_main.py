""" Tests of the command-line interface

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2021-08-12
:Copyright: 2021, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from biosimulators_masspy import __main__
from biosimulators_masspy import core
from biosimulators_masspy.data_model import KISAO_ALGORITHM_MAP
from biosimulators_utils.combine import data_model as combine_data_model
from biosimulators_utils.combine.exceptions import CombineArchiveExecutionError
from biosimulators_utils.combine.io import CombineArchiveWriter
from biosimulators_utils.config import get_config
from biosimulators_utils.report import data_model as report_data_model
from biosimulators_utils.report.io import ReportReader
from biosimulators_utils.simulator.exec import exec_sedml_docs_in_archive_with_containerized_simulator
from biosimulators_utils.simulator.specs import gen_algorithms_from_specs
from biosimulators_utils.sedml import data_model as sedml_data_model
from biosimulators_utils.sedml.io import SedmlSimulationWriter
from biosimulators_utils.sedml.utils import append_all_nested_children_to_doc
from biosimulators_utils.warnings import BioSimulatorsWarning
from kisao.exceptions import AlgorithmCannotBeSubstitutedException
from unittest import mock
import copy
import datetime
import dateutil.tz
import json
import mass
import numpy
import numpy.testing
import os
import shutil
import tempfile
import unittest
import yaml


class CliTestCase(unittest.TestCase):
    EXAMPLE_MODEL_FILENAME = os.path.join(os.path.dirname(__file__), 'fixtures', 'textbook.xml')
    NAMESPACES = {
        'sbml': 'http://www.sbml.org/sbml/level3/version1/core'
    }
    SPECIFICATIONS_FILENAME = os.path.join(os.path.dirname(__file__), '..', 'biosimulators.json')
    DOCKER_IMAGE = 'ghcr.io/biosimulators/biosimulators_masspy/masspy:latest'

    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_exec_sed_task_successfully(self):
        # configure simulation
        task = sedml_data_model.Task(
            model=sedml_data_model.Model(
                source=self.EXAMPLE_MODEL_FILENAME,
                language=sedml_data_model.ModelLanguage.SBML.value,
            ),
            simulation=sedml_data_model.UniformTimeCourseSimulation(
                initial_time=0.,
                output_start_time=0.,
                output_end_time=10.,
                number_of_points=10,
                algorithm=sedml_data_model.Algorithm(
                    kisao_id='KISAO_0000019',
                    changes=[
                        sedml_data_model.AlgorithmParameterChange(
                            kisao_id='KISAO_0000209',
                            new_value='1e-8',
                        )
                    ]
                ),
            ),
        )

        variables = [
            sedml_data_model.Variable(
                id='Time',
                symbol=sedml_data_model.Symbol.time,
                task=task),
            sedml_data_model.Variable(
                id='g6p',
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_g6p_c']",
                target_namespaces=self.NAMESPACES,
                task=task),
            sedml_data_model.Variable(
                id='f6p',
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_f6p_c']",
                target_namespaces=self.NAMESPACES,
                task=task),
            sedml_data_model.Variable(
                id='HEX1',
                target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_HEX1']",
                target_namespaces=self.NAMESPACES,
                task=task),
        ]

        # execute simulation
        variable_results, log = core.exec_sed_task(task, variables)

        # check that the simulation was executed correctly
        self.assertEqual(set(variable_results.keys()), set(['Time', 'g6p', 'f6p', 'HEX1']))
        for variable_result in variable_results.values():
            self.assertFalse(numpy.any(numpy.isnan(variable_result)))
        numpy.testing.assert_allclose(
            variable_results['Time'],
            numpy.linspace(
                task.simulation.output_start_time,
                task.simulation.output_end_time,
                task.simulation.number_of_points + 1,
            ))

        # check that log can be serialized to JSON
        self.assertEqual(log.algorithm, 'KISAO_0000019')
        self.assertEqual(log.simulator_details['integrator'], 'cvode')
        self.assertEqual(log.simulator_details['relative_tolerance'], 1e-8)

        json.dumps(log.to_json())

        log.out_dir = self.dirname
        log.export()
        with open(os.path.join(self.dirname, get_config().LOG_PATH), 'rb') as file:
            log_data = yaml.load(file, Loader=yaml.Loader)
        json.dumps(log_data)

    def test_exec_sed_task_non_zero_initial_time(self):
        # configure simulation
        task = sedml_data_model.Task(
            model=sedml_data_model.Model(
                source=self.EXAMPLE_MODEL_FILENAME,
                language=sedml_data_model.ModelLanguage.SBML.value,
            ),
            simulation=sedml_data_model.UniformTimeCourseSimulation(
                initial_time=10.,
                output_start_time=20.,
                output_end_time=30.,
                number_of_points=10,
                algorithm=sedml_data_model.Algorithm(
                    kisao_id='KISAO_0000019',
                ),
            ),
        )

        variables = [
            sedml_data_model.Variable(
                id='Time',
                symbol=sedml_data_model.Symbol.time,
                task=task),
            sedml_data_model.Variable(
                id='g6p',
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_g6p_c']",
                target_namespaces=self.NAMESPACES,
                task=task),
        ]

        # execute simulation
        variable_results, log = core.exec_sed_task(task, variables)

        # check that the simulation was executed correctly
        self.assertEqual(set(variable_results.keys()), set(['Time', 'g6p']))
        for variable_result in variable_results.values():
            self.assertFalse(numpy.any(numpy.isnan(variable_result)))
        numpy.testing.assert_allclose(
            variable_results['Time'],
            numpy.linspace(
                task.simulation.output_start_time,
                task.simulation.output_end_time,
                task.simulation.number_of_points + 1,
            ))

    def test_exec_sed_task_alt_alg(self):
        # configure simulation
        task = sedml_data_model.Task(
            model=sedml_data_model.Model(
                source=self.EXAMPLE_MODEL_FILENAME,
                language=sedml_data_model.ModelLanguage.SBML.value,
            ),
            simulation=sedml_data_model.UniformTimeCourseSimulation(
                initial_time=0.,
                output_start_time=0.,
                output_end_time=10.,
                number_of_points=10,
                algorithm=sedml_data_model.Algorithm(
                    kisao_id='KISAO_0000086',
                ),
            ),
        )

        variables = [
            sedml_data_model.Variable(
                id='Time',
                symbol=sedml_data_model.Symbol.time,
                task=task),
            sedml_data_model.Variable(
                id='g6p',
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_g6p_c']",
                target_namespaces=self.NAMESPACES,
                task=task),
        ]

        # execute simulation
        variable_results, log = core.exec_sed_task(task, variables)

        # check that the simulation was executed correctly
        self.assertEqual(set(variable_results.keys()), set(['Time', 'g6p']))
        for variable_result in variable_results.values():
            self.assertFalse(numpy.any(numpy.isnan(variable_result)), variable_result)
        numpy.testing.assert_allclose(
            variable_results['Time'],
            numpy.linspace(
                task.simulation.output_start_time,
                task.simulation.output_end_time,
                task.simulation.number_of_points + 1,
            ))

    def test_exec_sed_task_alg_substitution(self):
        # configure simulation
        task = sedml_data_model.Task(
            model=sedml_data_model.Model(
                source=self.EXAMPLE_MODEL_FILENAME,
                language=sedml_data_model.ModelLanguage.SBML.value,
            ),
            simulation=sedml_data_model.UniformTimeCourseSimulation(
                initial_time=0.,
                output_start_time=0.,
                output_end_time=10.,
                number_of_points=10,
                algorithm=sedml_data_model.Algorithm(
                    kisao_id='KISAO_0000019',
                ),
            ),
        )

        variables = [
            sedml_data_model.Variable(
                id='Time',
                symbol=sedml_data_model.Symbol.time,
                task=task),
            sedml_data_model.Variable(
                id='g6p',
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_g6p_c']",
                target_namespaces=self.NAMESPACES,
                task=task),
        ]

        # execute simulation
        task_2 = copy.deepcopy(task)
        task_2.simulation.algorithm.kisao_id = 'KISAO_0000560'
        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'NONE'}):
            with self.assertRaises(AlgorithmCannotBeSubstitutedException):
                core.exec_sed_task(task_2, variables)

        task_2 = copy.deepcopy(task)
        task_2.simulation.algorithm.kisao_id = 'KISAO_0000560'
        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'SIMILAR_VARIABLES'}):
            core.exec_sed_task(task_2, variables)

        task_2 = copy.deepcopy(task)
        task_2.simulation.algorithm.changes.append(sedml_data_model.AlgorithmParameterChange(
            kisao_id='KISAO_0000488',
            new_value='1',
        ))
        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'NONE'}):
            with self.assertRaises(NotImplementedError):
                core.exec_sed_task(task_2, variables)

        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'SIMILAR_VARIABLES'}):
            with self.assertWarns(BioSimulatorsWarning):
                core.exec_sed_task(task_2, variables)

        task_2 = copy.deepcopy(task)
        task_2.simulation.algorithm.changes.append(sedml_data_model.AlgorithmParameterChange(
            kisao_id='KISAO_0000209',
            new_value='abc',
        ))
        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'NONE'}):
            with self.assertRaises(ValueError):
                core.exec_sed_task(task_2, variables)

        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'SIMILAR_VARIABLES'}):
            with self.assertWarns(BioSimulatorsWarning):
                core.exec_sed_task(task_2, variables)

    def test_exec_sed_task_with_changes(self):
        # configure simulation
        task = sedml_data_model.Task(
            model=sedml_data_model.Model(
                source=self.EXAMPLE_MODEL_FILENAME,
                language=sedml_data_model.ModelLanguage.SBML.value,
            ),
            simulation=sedml_data_model.UniformTimeCourseSimulation(
                initial_time=0.,
                output_start_time=0.,
                output_end_time=10.,
                number_of_points=10,
                algorithm=sedml_data_model.Algorithm(
                    kisao_id='KISAO_0000019',
                ),
            ),
        )

        variables = [
            sedml_data_model.Variable(
                id='Time',
                symbol=sedml_data_model.Symbol.time,
                task=task),
        ]
        mass_model = mass.io.sbml.read_sbml_model(task.model.source)
        for met in mass_model.metabolites:
            task.model.changes.append(sedml_data_model.ModelAttributeChange(
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_{}']".format(met.id),
                target_namespaces=self.NAMESPACES,
                new_value=None))
            variables.append(sedml_data_model.Variable(
                id=met.id,
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_{}']".format(met.id),
                target_namespaces=self.NAMESPACES,
                task=task))

        task.model.changes.append(sedml_data_model.ModelAttributeChange(
            target="/sbml:sbml/sbml:model/sbml:listOfParameters/sbml:parameter[@id='v_R_HEX1']",
            target_namespaces=self.NAMESPACES,
            new_value=10))
        task.model.changes.append(sedml_data_model.ModelAttributeChange(
            target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_PFK_R01']/sbml:kineticLaw/sbml:listOfLocalParameters/sbml:localParameter[@id='Keq_PFK_A']",
            target_namespaces=self.NAMESPACES,
            new_value=20))
        task.model.changes.append(sedml_data_model.ModelAttributeChange(
            target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_SK_lac__L_c']/sbml:kineticLaw/sbml:listOfLocalParameters/sbml:localParameter[@id='lac__L_b']",
            target_namespaces=self.NAMESPACES,
            new_value=25))
        task.model.changes.append(sedml_data_model.ModelAttributeChange(
            target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ADK1']/sbml:kineticLaw/sbml:listOfLocalParameters/sbml:localParameter[@id='kf_R_ADK1']",
            target_namespaces=self.NAMESPACES,
            new_value=30))

        # execute simulation
        preprocessed_task = core.preprocess_sed_task(task, variables)

        task.model.changes.append(sedml_data_model.ModelAttributeChange(
            target="/sbml:sbml/sbml:model",
            target_namespaces=self.NAMESPACES,
            new_value=None))
        with self.assertRaises(ValueError):
            core.preprocess_sed_task(task, variables)

        task.model.changes = []
        results, _ = core.exec_sed_task(task, variables, preprocessed_task=preprocessed_task)

        for met in mass_model.metabolites:
            with self.assertRaises(AssertionError):
                numpy.testing.assert_allclose(results[met.id][0:task.simulation.number_of_points + 1],
                                              results[met.id][(-task.simulation.number_of_points + 1):])

        task.simulation.output_end_time = task.simulation.output_end_time / 2
        task.simulation.number_of_points = int(task.simulation.number_of_points / 2)
        results2, _ = core.exec_sed_task(task, variables, preprocessed_task=preprocessed_task)

        for met in mass_model.metabolites:
            numpy.testing.assert_allclose(results2[met.id], results[met.id][0:task.simulation.number_of_points + 1])
            task.model.changes.append(sedml_data_model.ModelAttributeChange(
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_{}']".format(met.id),
                target_namespaces=self.NAMESPACES,
                new_value=results2[met.id][-1]))

        results3, _ = core.exec_sed_task(task, variables, preprocessed_task=preprocessed_task)
        for met in mass_model.metabolites:
            numpy.testing.assert_allclose(results3[met.id], results[met.id][-(task.simulation.number_of_points + 1):])

        # parameters
        task.model.changes.append(sedml_data_model.ModelAttributeChange(
            target="/sbml:sbml/sbml:model/sbml:listOfParameters/sbml:parameter[@id='v_R_HEX1']",
            target_namespaces=self.NAMESPACES,
            new_value=10))
        task.model.changes.append(sedml_data_model.ModelAttributeChange(
            target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_PFK_R01']/sbml:kineticLaw/sbml:listOfLocalParameters/sbml:localParameter[@id='Keq_PFK_A']",
            target_namespaces=self.NAMESPACES,
            new_value=20))
        task.model.changes.append(sedml_data_model.ModelAttributeChange(
            target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_SK_lac__L_c']/sbml:kineticLaw/sbml:listOfLocalParameters/sbml:localParameter[@id='lac__L_b']",
            target_namespaces=self.NAMESPACES,
            new_value=25))
        task.model.changes.append(sedml_data_model.ModelAttributeChange(
            target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ADK1']/sbml:kineticLaw/sbml:listOfLocalParameters/sbml:localParameter[@id='kf_R_ADK1']",
            target_namespaces=self.NAMESPACES,
            new_value=30))

        core.exec_sed_task(task, variables, preprocessed_task=preprocessed_task)
        self.assertEqual(preprocessed_task['model']['model'].parameters['v']['v_HEX1'], 10)
        self.assertEqual(preprocessed_task['model']['model'].parameters['Custom']['Keq_PFK_A'], 20)
        self.assertEqual(preprocessed_task['model']['model'].parameters['Boundary']['lac__L_b'], 25)
        self.assertEqual(preprocessed_task['model']['model'].parameters['kf']['kf_ADK1'], 30)

    def test_exec_sed_task_sim_error_handling(self):
        # configure simulation
        task = sedml_data_model.Task(
            model=sedml_data_model.Model(
                source=self.EXAMPLE_MODEL_FILENAME,
                language=sedml_data_model.ModelLanguage.SBML.value,
            ),
            simulation=sedml_data_model.UniformTimeCourseSimulation(
                initial_time=0.,
                output_start_time=0.,
                output_end_time=10.,
                number_of_points=10,
                algorithm=sedml_data_model.Algorithm(
                    kisao_id='KISAO_0000030',
                ),
            ),
        )

        variables = [
            sedml_data_model.Variable(
                id='Time',
                symbol=sedml_data_model.Symbol.time,
                task=task),
            sedml_data_model.Variable(
                id='g6p',
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_g6p_c']",
                target_namespaces=self.NAMESPACES,
                task=task),
        ]

        # execute simulation
        with self.assertRaisesRegex(ValueError, 'Simulation failed'):
            core.exec_sed_task(task, variables)

    def test_exec_sed_task_error_handling(self):
        # configure simulation
        task = sedml_data_model.Task(
            model=sedml_data_model.Model(
                source=self.EXAMPLE_MODEL_FILENAME,
                language=sedml_data_model.ModelLanguage.SBML.value,
            ),
            simulation=sedml_data_model.UniformTimeCourseSimulation(
                initial_time=0.,
                output_start_time=0.,
                output_end_time=10.,
                number_of_points=10,
                algorithm=sedml_data_model.Algorithm(
                    kisao_id='KISAO_0000019',
                ),
            ),
        )

        variables = [
            sedml_data_model.Variable(
                id='Time',
                symbol=sedml_data_model.Symbol.time,
                task=task),
            sedml_data_model.Variable(
                id='g6p',
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_g6p_c']",
                target_namespaces=self.NAMESPACES,
                task=task),
        ]

        # execute simulation
        variables_2 = copy.deepcopy(variables)
        variables_2[0].symbol = 'mass'
        with self.assertRaisesRegex(NotImplementedError, 'The following symbols are not supported'):
            core.exec_sed_task(task, variables_2)

        variables_2 = copy.deepcopy(variables)
        variables_2[1].target = '/sbml:sbml'
        with self.assertRaisesRegex(ValueError, 'The following targets are not supported'):
            core.exec_sed_task(task, variables_2)

        task_2 = copy.deepcopy(task)
        task_2.simulation.output_start_time = 1.5
        with self.assertRaisesRegex(NotImplementedError, 'must be an integer'):
            core.exec_sed_task(task_2, variables)

    def test_exec_sedml_docs_in_combine_archive_successfully(self):
        doc, archive_filename = self._build_combine_archive()

        out_dir = os.path.join(self.dirname, 'out')

        config = get_config()
        config.REPORT_FORMATS = [report_data_model.ReportFormat.h5]
        config.BUNDLE_OUTPUTS = True
        config.KEEP_INDIVIDUAL_OUTPUTS = True

        _, log = core.exec_sedml_docs_in_combine_archive(archive_filename, out_dir, config=config)
        if log.exception:
            raise log.exception

        self._assert_combine_archive_outputs(doc, out_dir)

    def _build_combine_archive(self, algorithm=None):
        doc = self._build_sed_doc(algorithm=algorithm)

        archive_dirname = os.path.join(self.dirname, 'archive')
        if not os.path.isdir(archive_dirname):
            os.mkdir(archive_dirname)

        model_filename = os.path.join(archive_dirname, 'model.xml')
        shutil.copyfile(self.EXAMPLE_MODEL_FILENAME, model_filename)

        sim_filename = os.path.join(archive_dirname, 'sim.sedml')
        SedmlSimulationWriter().run(doc, sim_filename)

        archive = combine_data_model.CombineArchive(
            contents=[
                combine_data_model.CombineArchiveContent(
                    'model.xml', combine_data_model.CombineArchiveContentFormat.SBML.value),
                combine_data_model.CombineArchiveContent(
                    'sim.sedml', combine_data_model.CombineArchiveContentFormat.SED_ML.value),
            ],
        )
        archive_filename = os.path.join(self.dirname, 'archive.omex')
        CombineArchiveWriter().run(archive, archive_dirname, archive_filename)

        return (doc, archive_filename)

    def _build_sed_doc(self, algorithm=None):
        if algorithm is None:
            algorithm = sedml_data_model.Algorithm(
                kisao_id='KISAO_0000019',
            )

        doc = sedml_data_model.SedDocument()
        doc.models.append(sedml_data_model.Model(
            id='model',
            source='model.xml',
            language=sedml_data_model.ModelLanguage.SBML.value,
        ))
        doc.simulations.append(sedml_data_model.UniformTimeCourseSimulation(
            id='sim_time_course',
            initial_time=0.,
            output_start_time=0.,
            output_end_time=10.,
            number_of_points=10,
            algorithm=algorithm,
        ))

        doc.tasks.append(sedml_data_model.Task(
            id='task_1',
            model=doc.models[0],
            simulation=doc.simulations[0],
        ))

        doc.data_generators.append(sedml_data_model.DataGenerator(
            id='data_gen_time',
            variables=[
                sedml_data_model.Variable(
                    id='var_time',
                    symbol=sedml_data_model.Symbol.time.value,
                    task=doc.tasks[0],
                ),
            ],
            math='var_time',
        ))
        doc.data_generators.append(sedml_data_model.DataGenerator(
            id='data_gen_g6p',
            variables=[
                sedml_data_model.Variable(
                    id='var_g6p',
                    target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_g6p_c']",
                    target_namespaces=self.NAMESPACES,
                    task=doc.tasks[0],
                ),
            ],
            math='var_g6p',
        ))

        doc.outputs.append(sedml_data_model.Report(
            id='report',
            data_sets=[
                sedml_data_model.DataSet(id='data_set_time', label='Time', data_generator=doc.data_generators[0]),
                sedml_data_model.DataSet(id='data_set_g6p', label='g6p', data_generator=doc.data_generators[1]),
            ],
        ))

        append_all_nested_children_to_doc(doc)

        return doc

    def _assert_combine_archive_outputs(self, doc, out_dir):
        self.assertEqual(set(['reports.h5']).difference(set(os.listdir(out_dir))), set())

        report = ReportReader().run(doc.outputs[0], out_dir, 'sim.sedml/report', format=report_data_model.ReportFormat.h5)

        self.assertEqual(sorted(report.keys()), sorted([d.id for d in doc.outputs[0].data_sets]))

        sim = doc.tasks[0].simulation
        self.assertEqual(len(report[doc.outputs[0].data_sets[0].id]), sim.number_of_points + 1)

        for data_set_result in report.values():
            self.assertFalse(numpy.any(numpy.isnan(data_set_result)))

        self.assertIn('data_set_time', report)
        numpy.testing.assert_allclose(report[doc.outputs[0].data_sets[0].id],
                                      numpy.linspace(sim.output_start_time, sim.output_end_time, sim.number_of_points + 1))

    def test_exec_sedml_docs_in_combine_archive_with_all_algorithms(self):
        failures = []
        for alg in gen_algorithms_from_specs(self.SPECIFICATIONS_FILENAME).values():
            doc, archive_filename = self._build_combine_archive(algorithm=alg)
            out_dir = os.path.join(self.dirname, alg.kisao_id)

            config = get_config()
            config.REPORT_FORMATS = [report_data_model.ReportFormat.h5]
            config.BUNDLE_OUTPUTS = True
            config.KEEP_INDIVIDUAL_OUTPUTS = True

            try:
                _, log = core.exec_sedml_docs_in_combine_archive(archive_filename, out_dir, config=config)
                if log.exception:
                    raise log.exception
                self._assert_combine_archive_outputs(doc, out_dir)
            except CombineArchiveExecutionError as exception:
                if 'Simulation failed' in str(exception):
                    failures.append(alg.kisao_id)
                else:
                    raise

        self.assertEqual(failures, ['KISAO_0000030', 'KISAO_0000032'])  # model can't be executed with these algorithms

    def test_exec_sedml_docs_in_combine_archive_with_cli(self):
        doc, archive_filename = self._build_combine_archive()
        out_dir = os.path.join(self.dirname, 'out')
        env = self._get_combine_archive_exec_env()

        with mock.patch.dict(os.environ, env):
            with __main__.App(argv=['-i', archive_filename, '-o', out_dir]) as app:
                app.run()

        self._assert_combine_archive_outputs(doc, out_dir)

    def _get_combine_archive_exec_env(self):
        return {
            'REPORT_FORMATS': 'h5'
        }

    def test_exec_sedml_docs_in_combine_archive_with_docker_image(self):
        doc, archive_filename = self._build_combine_archive()
        out_dir = os.path.join(self.dirname, 'out')
        docker_image = self.DOCKER_IMAGE
        env = self._get_combine_archive_exec_env()

        exec_sedml_docs_in_archive_with_containerized_simulator(
            archive_filename, out_dir, docker_image, environment=env, pull_docker_image=False)

        self._assert_combine_archive_outputs(doc, out_dir)
