import copy
import json
import math
from pathlib import Path
import unittest

from open_cae_evidence.model import audit_deck, make_deck, validate
from open_cae_evidence.reference import BETAS, mac, references, shape

CONFIG = Path(__file__).resolve().parents[1] / "benchmarks/cantilever.json"


class InputTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads(CONFIG.read_text())

    def test_negative_physical_inputs(self):
        cases = [("material", "density", -7800), ("material", "density", 7.8),
                 ("material", "young_modulus", 210000), ("material", "young_modulus", 0),
                 ("material", "poisson_ratio", 0.5), ("material", "poisson_ratio", -1),
                 ("material", "density", float("nan")), ("material", "young_modulus", float("inf")),
                 ("geometry", "length", -1), ("geometry", "length", 1000),
                 ("geometry", "height", 0.2), ("geometry", "width", True)]
        for group, name, value in cases:
            with self.subTest(group=group, name=name, value=value):
                corrupt = copy.deepcopy(self.config)
                corrupt[group][name] = value
                with self.assertRaises(ValueError):
                    validate(corrupt)

    def test_unit_and_boundary_corruption(self):
        for group, key, value in [("units", "length", "mm"), ("units", "stress", "MPa"),
                                  ("units", "density", "g/cm^3"), ("boundary", "plane", "x=1"),
                                  ("boundary", "dofs", [1, 2]), ("boundary", "dofs", [])]:
            with self.subTest(group=group, key=key):
                corrupt = copy.deepcopy(self.config)
                corrupt[group][key] = value
                with self.assertRaises(ValueError):
                    validate(corrupt)

    def test_bad_refinement(self):
        for meshes in [[[12, 2, 2]], [[12, 2, 2], [12, 4, 4], [48, 6, 6]],
                       [[12, 2, 2], [24, 3, 4], [48, 6, 6]], [[12, 2, 2], [24, 4, 4], [128, 128, 128]]]:
            with self.subTest(meshes=meshes):
                corrupt = copy.deepcopy(self.config)
                corrupt["meshes"] = meshes
                with self.assertRaises(ValueError):
                    validate(corrupt)

    def test_mass_and_fixed_face(self):
        deck, nodes, audit = make_deck(self.config, [12, 2, 2])
        self.assertAlmostEqual(audit["mass_kg"], 1.17, places=12)
        self.assertAlmostEqual(audit["volume_m3"], 0.00015, places=12)
        self.assertEqual(audit["fixed_nodes"], 21)
        self.assertEqual(audit["fixed_nodes"], sum(x == 0 for x, _, _ in nodes.values()))

    def test_corrupted_serialized_deck_is_detected(self):
        deck, _, _ = make_deck(self.config, [12, 2, 2])
        corruptions = [deck.replace("FIXED, 1, 3", "FIXED, 1, 2"),
                       deck.replace("FIXED, 1, 3", "ALL, 1, 3"),
                       deck.replace("*BOUNDARY\nFIXED, 1, 3\n", ""),
                       deck.replace("210000000000, 0.3", "210000, 0.3"),
                       deck.replace("*DENSITY\n7800", "*DENSITY\n-7800"),
                       deck.replace("*FREQUENCY\n6", "*FREQUENCY\n0"),
                       deck.replace("*NSET, NSET=FIXED\n1,", "*NSET, NSET=FIXED\n2,")]
        # Reverse the first element's corner orientation while keeping its IDs.
        lines = deck.splitlines()
        start = lines.index("*ELEMENT, TYPE=C3D20R, ELSET=BEAM") + 1
        row = lines[start].split(",")
        row[1], row[2] = row[2], row[1]
        lines[start] = ",".join(row)
        corruptions.append("\n".join(lines))
        for corrupt in corruptions:
            with self.subTest(corruption=corrupt[-200:]):
                with self.assertRaises(ValueError):
                    audit_deck(corrupt, self.config, [12, 2, 2])


class AnalyticalTests(unittest.TestCase):
    def test_characteristic_roots(self):
        for beta in BETAS:
            # Scaled equation avoids magnifying cancellation for higher roots.
            self.assertAlmostEqual(math.cos(beta) + 1 / math.cosh(beta), 0, places=13)

    def test_shape_end_conditions(self):
        for beta in BETAS[:3]:
            ratio = (math.cosh(beta) + math.cos(beta)) / (math.sinh(beta) + math.sin(beta))
            self.assertAlmostEqual(shape(beta, 0), 0)
            self.assertAlmostEqual(math.cosh(beta) + math.cos(beta) - ratio * (math.sinh(beta) + math.sin(beta)), 0, places=12)
            self.assertAlmostEqual(math.sinh(beta) - math.sin(beta) - ratio * (math.cosh(beta) + math.cos(beta)), 0, places=12)
            # phi'(0)=beta*(0+0-ratio*(1-1))=0 exactly.

    def test_first_frequency_from_closed_form_and_plane_ratio(self):
        config = json.loads(CONFIG.read_text())
        ref = references(config)
        expected = BETAS[0]**2 / (2 * math.pi) * 0.01 / math.sqrt(12) * math.sqrt(210e9 / 7800)
        self.assertAlmostEqual(ref[0]["frequency_hz"], expected, places=12)
        self.assertAlmostEqual(ref[1]["frequency_hz"] / ref[0]["frequency_hz"], 1.5)

    def test_mac_detects_wrong_shape_and_is_sign_invariant(self):
        a = [shape(BETAS[0], i / 24) for i in range(25)]
        b = [shape(BETAS[1], i / 24) for i in range(25)]
        self.assertAlmostEqual(mac(a, [-2 * v for v in a]), 1)
        self.assertLess(mac(a, b), 0.1)
        with self.assertRaises(ValueError):
            mac(a, [0] * len(a))


if __name__ == "__main__":
    unittest.main()
