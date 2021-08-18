from biosimulators_utils.sedml import data_model
from biosimulators_utils.sedml.io import SedmlSimulationWriter
import json
import mass
import mass.test
import matplotlib
import os

mass_model = mass.test.create_test_model("WholeCellRBC_MA_Rates")
sim = mass.Simulation(reference_model=mass_model, verbose=True)
concs, fluxes = sim.simulate(mass_model, time=(0, 1000,))

matplotlib.use('agg')
matplotlib.pyplot.plot(concs._time, concs['g6p_c'])
matplotlib.pyplot.savefig(
    'figure.png')


namespaces = {
    'sbml': "http://www.sbml.org/sbml/level3/version1/core",
}

doc = data_model.SedDocument()

model = data_model.Model(id="model", language="urn:sedml:language:sbml", source="WholeCellRBC_MA_Rates.xml")
doc.models.append(model)

model.changes.append(data_model.ModelAttributeChange(
    target="/sbml:sbml/sbml:model/sbml:listOfParameters/sbml:parameter[@id='cobra_default_lb']/@value",
    target_namespaces=namespaces,
    new_value="-1000"))

sim = data_model.UniformTimeCourseSimulation(
    id='simulation',
    initial_time=0,
    output_start_time=0.,
    output_end_time=1000.,
    number_of_steps=1000,
    algorithm=data_model.Algorithm(kisao_id='KISAO_0000019')
)
doc.simulations.append(sim)

task = data_model.Task(id='task', model=model, simulation=sim)
doc.tasks.append(task)

expected_reports = []

time_var = data_model.Variable(
    id='var_time',
    name='Time',
    symbol=data_model.Symbol.time.value,
    task=task,
)
time_data_gen = data_model.DataGenerator(
    id='data_generator_time',
    name='Time',
    variables=[time_var],
    math=time_var.id,
)
doc.data_generators.append(time_data_gen)

report = data_model.Report(id='metabolite_concentrations', name='Metabolite concentrations')
report.data_sets.append(data_model.DataSet(
    id='data_set_time_metabolite_concentrations',
    label='Time',
    name='Time',
    data_generator=time_data_gen,
))
doc.outputs.append(report)

expected_report = {
    "id": "simulation.sedml/" + report.id,
    "dataSets": [
        {"id": "data_set_time", "label": "Time"},
    ],
    "points": [sim.number_of_steps + 1],
    "values": [{
        "id": "data_set_time",
        "label": "time",
        "value": {
            "0": sim.output_start_time,
            "1000": sim.output_end_time
        }
    }]
}
expected_reports.append(expected_report)

for met in mass_model.metabolites:
    var = data_model.Variable(
        id='M_' + met.id,
        name=met.name,
        target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_{}']".format(met.id),
        target_namespaces=namespaces,
        task=task,
    )
    data_gen = data_model.DataGenerator(
        id='data_generator_M_' + met.id,
        name=met.name,
        variables=[var],
        math=var.id,
    )
    doc.data_generators.append(data_gen)

    data_set = data_model.DataSet(
        id='data_set_M_' + met.id,
        label='{} ({})'.format(met.name, met.compartment),
        name='{} ({})'.format(met.name, met.compartment),
        data_generator=data_gen)
    report.data_sets.append(data_set)

    expected_report['dataSets'].append({
        "id": data_set.id,
        "label": data_set.label,
    })

report = data_model.Report(id='reaction_fluxes', name='Reaction fluxes')
report.data_sets.append(data_model.DataSet(
    id='data_set_time_reaction_fluxes',
    label='Time',
    name='Time',
    data_generator=time_data_gen,
))
doc.outputs.append(report)

expected_report = {
    "id": "simulation.sedml/" + report.id,
    "dataSets": [
        {"id": "data_set_time", "label": "Time"},
    ],
    "points": [sim.number_of_steps + 1],
    "values": [{
        "id": "data_set_time",
        "label": "time",
        "value": {
            "0": sim.output_start_time,
            "1000": sim.output_end_time
        }
    }]
}
expected_reports.append(expected_report)

rxn_name_counts = {}
rxn_labels = {}
rxn_indices = {}
for rxn in mass_model.reactions:
    if rxn.name not in rxn_name_counts:
        rxn_name_counts[rxn.name] = 0
    rxn_name_counts[rxn.name] += 1
    rxn_indices[rxn.id] = rxn_name_counts[rxn.name]

for rxn in mass_model.reactions:
    var = data_model.Variable(
        id='R_' + rxn.id,
        name=rxn.name,
        target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_{}']".format(rxn.id),
        target_namespaces=namespaces,
        task=task,
    )
    data_gen = data_model.DataGenerator(
        id='data_generator_R_' + rxn.id,
        name=rxn.name,
        variables=[var],
        math=var.id,
    )
    doc.data_generators.append(data_gen)

    data_set = data_model.DataSet(
        id='data_set_R_' + rxn.id,
        label=rxn.name if rxn_name_counts[rxn.name] == 1 else '{} {}'.format(rxn.name, rxn_indices[rxn.id]),
        name=rxn.name,
        data_generator=data_gen)
    report.data_sets.append(data_set)

    expected_report['dataSets'].append({
        "id": data_set.id,
        "label": data_set.label,
    })

os.chdir(os.path.expanduser('~/Documents/Biosimulators_test_suite/examples/mass/Bordbar-Cell-Syst-2015-RBC-metabolism/'))
SedmlSimulationWriter().run(doc, 'simulation.sedml')
