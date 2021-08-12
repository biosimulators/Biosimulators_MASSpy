from biosimulators_masspy.data_model import KISAO_ALGORITHM_MAP
from biosimulators_utils.utils.core import parse_value
import json
import os
import unittest


class DataModelTestCase(unittest.TestCase):
    def test_data_model_matches_specs(self):
        specs_filename = os.path.join(os.path.dirname(__file__), '..', 'biosimulators.json')
        with open(specs_filename, 'r') as file:
            specs = json.load(file)

        alg_kisao_ids = [alg_specs['kisaoId']['id'] for alg_specs in specs['algorithms']]
        self.assertEqual(len(alg_kisao_ids), len(set(alg_kisao_ids)))
        self.assertEqual(set(KISAO_ALGORITHM_MAP.keys()),
                         set(alg_kisao_ids))

        for alg_specs in specs['algorithms']:
            alg_props = KISAO_ALGORITHM_MAP[alg_specs['kisaoId']['id']]
            param_kisao_ids = [param_specs['kisaoId']['id'] for param_specs in alg_specs['parameters']]
            self.assertEqual(len(param_kisao_ids), len(set(param_kisao_ids)))
            self.assertEqual(set(alg_props['parameters'].keys()),
                             set(param_kisao_ids),
                             'Algorithm: `{}`'.format(alg_specs['kisaoId']['id']))

            param_ids = [param_specs['id'] for param_specs in alg_specs['parameters'] if param_specs['id'] is not None]
            self.assertEqual(len(param_ids), len(set(param_ids)), 'Algorithm: `{}`'.format(alg_specs['kisaoId']['id']))

            for param_specs in alg_specs['parameters']:
                param_props = alg_props['parameters'][param_specs['kisaoId']['id']]
                self.assertEqual(param_props['type'].value, param_specs['type'])
                self.assertEqual(param_props['default'],
                                 None if param_specs['value'] is None else parse_value(param_specs['value'], param_specs['type']))
                self.assertEqual(param_props['id'], param_specs.get('id', None),
                                 'Algorithm: `{}`, Parameter: `{}`'.format(alg_specs['kisaoId']['id'], param_specs['kisaoId']['id'])
                                 )
