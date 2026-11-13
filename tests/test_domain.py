"""验证现有领域摘要的业务边界。"""
import copy
import json
import unittest
from pathlib import Path

from warning_policy.analysis import water_ranges


class DomainTests(unittest.TestCase):
    def setUp(self):
        self.tables = json.loads((Path(__file__).resolve().parents[1]/'fixtures/dataset.json').read_text())['tables']

    def test_missing_datum_is_explicitly_unsupported(self):
        tables = copy.deepcopy(self.tables)
        tables['datums'] = []
        result = water_ranges(tables)
        self.assertFalse(result['gauges'])
        self.assertEqual(len(tables['levels']),len(result['unsupported_records']))

    def test_offset_is_added_before_range_aggregation(self):
        tables = copy.deepcopy(self.tables)
        tables['levels'] = [dict(tables['levels'][0],level_m=2,datum_version=1)]
        result = water_ranges(tables)
        self.assertAlmostEqual(2.1,result['gauges']['g0']['minimum_m'])

    def test_null_level_is_not_zero(self):
        tables = copy.deepcopy(self.tables)
        tables['levels'] = [dict(tables['levels'][0],level_m=None)]
        self.assertFalse(water_ranges(tables)['gauges'])

    def test_datum_order_does_not_change_range(self):
        tables = copy.deepcopy(self.tables)
        tables['datums'].reverse()
        self.assertEqual(water_ranges(self.tables),water_ranges(tables))
