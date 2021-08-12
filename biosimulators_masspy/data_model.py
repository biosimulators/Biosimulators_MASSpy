""" Data model

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2021-08-12
:Copyright: 2021, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from biosimulators_utils.data_model import ValueType
import collections

__all__ = [
    'KISAO_ALGORITHM_MAP',
]


KISAO_ALGORITHM_MAP = collections.OrderedDict([
    ('KISAO_0000019', {
        'kisao_id': 'KISAO_0000019',
        'id': 'cvode',
        'name': "CVODE",
        'parameters': {
            'KISAO_0000209': {
                'kisao_id': 'KISAO_0000209',
                'id': 'relative_tolerance',
                'name': 'relative tolerance',
                'type': ValueType.float,
                'default': 0.000001,
            },
            'KISAO_0000211': {
                'kisao_id': 'KISAO_0000211',
                'id': 'absolute_tolerance',
                'name': 'absolute tolerance',
                'type': ValueType.float,
                'default': 1e-12,
            },
            'KISAO_0000220': {
                'kisao_id': 'KISAO_0000220',
                'id': 'maximum_bdf_order',
                'name': 'Maximum Backward Differentiation Formula (BDF) order',
                'type': ValueType.integer,
                'default': 5,
            },
            'KISAO_0000219': {
                'kisao_id': 'KISAO_0000219',
                'id': 'maximum_adams_order',
                'name': 'Maximum Adams order',
                'type': ValueType.integer,
                'default': 12,
            },
            'KISAO_0000415': {
                'kisao_id': 'KISAO_0000415',
                'id': 'maximum_num_steps',
                'name': 'Maximum number of steps',
                'type': ValueType.integer,
                'default': 20000,
            },
            'KISAO_0000467': {
                'kisao_id': 'KISAO_0000467',
                'id': 'maximum_time_step',
                'name': 'Maximum time step',
                'type': ValueType.float,
                'default': None,
            },
            'KISAO_0000485': {
                'kisao_id': 'KISAO_0000485',
                'id': 'minimum_time_step',
                'name': 'Minimum time step',
                'type': ValueType.float,
                'default': None,
            },
            'KISAO_0000332': {
                'kisao_id': 'KISAO_0000332',
                'id': 'initial_time_step',
                'name': 'Initial time step',
                'type': ValueType.float,
                'default': None,
            },
            'KISAO_0000671': {
                'kisao_id': 'KISAO_0000671',
                'id': 'stiff',
                'name': 'Stiff',
                'type': ValueType.boolean,
                'default': True,
            },
            'KISAO_0000670': {
                'kisao_id': 'KISAO_0000670',
                'id': 'multiple_steps',
                'name': 'Multiple steps',
                'type': ValueType.boolean,
                'default': False,
            },
        },
    }),
    ('KISAO_0000030', {
        'kisao_id': 'KISAO_0000030',
        'id': 'euler',
        'name': "Forward Euler method",
        'parameters': {}
    }),
    ('KISAO_0000032', {
        'kisao_id': 'KISAO_0000032',
        'id': 'rk4',
        'name': "Runge-Kutta fourth order method",
        'parameters': {}
    }),
    ('KISAO_0000086', {
        'kisao_id': 'KISAO_0000086',
        'id': 'rk45',
        'name': "Fehlberg method",
        'parameters': {
            'KISAO_0000467': {
                'kisao_id': 'KISAO_0000467',
                'id': 'maximum_time_step',
                'name': 'Maximum time step',
                'type': ValueType.float,
                'default': 1.0,
            },
            'KISAO_0000485': {
                'kisao_id': 'KISAO_0000485',
                'id': 'minimum_time_step',
                'name': 'Minimum time step',
                'type': ValueType.float,
                'default': 1e-12,
            },
            'KISAO_0000597': {
                'kisao_id': 'KISAO_0000597',
                'id': 'epsilon',
                'name': 'Epsilon',
                'type': ValueType.float,
                'default': 0.000000000001,
            },
        }
    }),
])
