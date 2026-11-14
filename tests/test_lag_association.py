"""验证滞后关联研究的统计边界与可解释输出。"""
import copy
import json
import unittest
from pathlib import Path

from warning_policy.lag_association import analyze_lag_association

ROOT = Path(__file__).resolve().parents[1]
FAST = {"permutations": 40}


def fixture_tables():
    return json.loads((ROOT / "fixtures/dataset.json").read_text())["tables"]


class LagAssociationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tables = fixture_tables()
        cls.study = analyze_lag_association(cls.tables, FAST)

    def test_scope_explicitly_refuses_causal_and_forecast_claims(self):
        scope = self.study["scope"]
        self.assertIn("不证明", scope)
        self.assertIn("预测精度", scope)
        self.assertTrue(any("因果" in item for item in self.study["not_claimed"]))
        self.assertTrue(any("插值" in item for item in self.study["limitations"]))

    def test_each_gauge_reports_explainable_group_and_lag_range(self):
        self.assertEqual({"g0", "g1"}, {item["gauge_id"] for item in self.study["gauges"]})
        for gauge in self.study["gauges"]:
            self.assertEqual("supported_candidate", gauge["status"])
            selection = gauge["selection"]
            self.assertEqual(["s1"], selection["stations"])
            self.assertEqual(["c0"] if gauge["gauge_id"] == "g0" else ["c1"],
                             selection["catchments"])
            low, high = selection["lag_minutes"]
            self.assertLess(low, high)
            effective = gauge["effective_lag"]
            self.assertLessEqual(effective["low_minutes"], low)
            self.assertGreaterEqual(effective["high_minutes"], high)
            self.assertGreaterEqual(gauge["search_adjusted_p"], gauge["p_resolution"])

    def test_search_corrected_p_is_above_single_lag_resolution_and_below_alpha(self):
        for gauge in self.study["gauges"]:
            self.assertGreater(gauge["null_max_median"], 0.0)
            self.assertLess(gauge["null_max_median"], gauge["selection"]["r"])
            self.assertLessEqual(gauge["search_adjusted_p"], self.study["config"]["alpha"])
            self.assertAlmostEqual(gauge["p_resolution"],
                                   1.0 / (FAST["permutations"] + 1), places=3)

    def test_selection_and_validation_use_disjoint_complete_processes(self):
        data = self.study["data"]
        self.assertFalse(set(data["selection_processes"]) & set(data["validation_processes"]))
        self.assertEqual(set(data["processes"]),
                         set(data["selection_processes"]) | set(data["validation_processes"]))
        for gauge in self.study["gauges"]:
            self.assertTrue(gauge["validation"]["stable"])
            self.assertEqual(len(gauge["validation"]["per_process"]),
                             len(gauge["validation"]["processes"]))
            self.assertEqual(1.0, gauge["validation"]["positive_share"])

    def test_upstream_candidates_ignore_downstream_station_for_selection(self):
        for gauge in self.study["gauges"]:
            self.assertEqual(["s0", "s1"], gauge["candidate_stations"])
            self.assertEqual("network.upstream", gauge["candidate_basis"])
            for candidate in gauge["candidate_groups"]:
                self.assertTrue(set(candidate["stations"]) <= {"s0", "s1"})

    def test_correlated_upstream_stations_are_reported_as_redundant(self):
        for gauge in self.study["gauges"]:
            self.assertTrue(gauge["redundant_stations"])
            value = gauge["station_correlations"]["s0|s1"]["r"]
            self.assertGreater(value, self.study["config"]["redundant_station_r"])
            self.assertLess(gauge["selection"]["r"], gauge["candidate_groups"][0]["r"] + 1e-9)

    def test_shuffled_change_levels_do_not_form_a_lag_candidate(self):
        tables = copy.deepcopy(self.tables)
        by_process = {}
        for row in tables["levels"]:
            by_process.setdefault(row["process_id"], []).append(row)
        import random
        rng = random.Random(20250909)
        for rows in by_process.values():
            # 固定种子的过程内重排：保留每场水位分布，但打散与雨量的时刻对应。
            values = [row["level_m"] for row in rows]
            rng.shuffle(values)
            for row, value in zip(rows, values):
                row["level_m"] = value
        study = analyze_lag_association(tables, FAST)
        for gauge in study["gauges"]:
            self.assertNotEqual("supported_candidate", gauge["status"])
        self.assertTrue(any(gauge["status"] in {"below_peak_threshold", "not_significant_after_search"}
                            for gauge in study["gauges"]))

    def test_missing_datum_is_excluded_rather_than_shifted(self):
        tables = copy.deepcopy(self.tables)
        tables["datums"] = [row for row in tables["datums"] if row["datum_version"] != 2]
        study = analyze_lag_association(tables, FAST)
        self.assertEqual({"missing_datum_offset": 288}, study["data"]["level_excluded"])
        # 只剩一个基准版本时，可用过程变少，但选型与验证仍各自使用完整过程。
        self.assertEqual(12, len(study["data"]["processes"]))
        self.assertEqual(6, len(study["data"]["selection_processes"]))
        self.assertEqual(6, len(study["data"]["validation_processes"]))
        for gauge in study["gauges"]:
            self.assertEqual(6, gauge["selection"]["processes"])
            self.assertEqual(["1"], sorted({
                item["datum_version"] for item in gauge["validation"]["per_process"]
            }))

    def test_short_processes_do_not_produce_a_search(self):
        tables = copy.deepcopy(self.tables)
        tables["levels"] = [row for row in tables["levels"]
                            if row["process_id"] in {"p0", "p1", "p2"}]
        study = analyze_lag_association(tables, FAST)
        self.assertEqual([], study["data"]["selection_processes"])
        for gauge in study["gauges"]:
            self.assertEqual("insufficient_processes", gauge["status"])
            self.assertIsNone(gauge["selection"])

    def test_result_is_independent_of_row_order(self):
        tables = copy.deepcopy(self.tables)
        for rows in tables.values():
            rows.reverse()
        other = analyze_lag_association(tables, FAST)
        self.assertEqual(json.dumps(self.study, sort_keys=True),
                         json.dumps(other, sort_keys=True))

    def test_common_seasonal_rise_alone_does_not_create_a_lag_candidate(self):
        tables = copy.deepcopy(self.tables)
        import random
        rng = random.Random(11)
        first = {}
        for row in tables["rain"]:
            first.setdefault(row["process_id"], row["at"])
        order = {process: index for index, process in
                 enumerate(sorted(first, key=lambda name: first[name]))}
        for row in tables["rain"]:
            row["rain_mm"] += order[row["process_id"]] * 2.0
        grouped = {}
        for row in tables["levels"]:
            grouped.setdefault(row["process_id"], []).append(row)
        for process, rows in grouped.items():
            # 季节同步上抬水位，但打散过程内的雨—水位对应。
            values = [row["level_m"] + order[process] * 0.05 for row in rows]
            rng.shuffle(values)
            for row, value in zip(rows, values):
                row["level_m"] = value
        study = analyze_lag_association(tables, FAST)
        for gauge in study["gauges"]:
            self.assertNotEqual("supported_candidate", gauge["status"])

    def test_zero_lag_reference_is_reported_and_lower_than_peak(self):
        for gauge in self.study["gauges"]:
            zero = gauge["zero_lag_reference"]
            self.assertEqual(0, gauge["lag_profile"][0]["low_minutes"])
            self.assertLess(zero["r"], gauge["selection"]["r"] - self.study["config"]["min_lag_gain"])
            self.assertIn("滞后", zero["note"])

    def test_invalid_config_is_rejected_instead_of_silently_ignored(self):
        for config in ({"unknown_option": 1}, {"bin_minutes": 0}, {"permutations": -5}):
            with self.assertRaises(ValueError):
                analyze_lag_association(self.tables, config)

    def test_null_family_uses_rotated_levels_not_fresh_search(self):
        # 零分布中位数应显著低于观测峰值，说明搜索校正确有惩罚作用。
        for gauge in self.study["gauges"]:
            self.assertLess(gauge["null_max_quantile"], gauge["selection"]["r"])


if __name__ == "__main__":
    unittest.main()
