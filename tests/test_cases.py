import unittest

from warning_policy.cases import EventCase, classify


class PolicyTests(unittest.TestCase):
    def test_unobserved_impact_is_not_a_negative(self):
        self.assertEqual("unknown", classify(EventCase("e1", 0.8, None), 0.5))
