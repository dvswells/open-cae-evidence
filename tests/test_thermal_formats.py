"""Actual captured 2.21 output prevents a regression hidden by local 2.20."""

import hashlib
import json
from pathlib import Path
import unittest

from open_cae_evidence.thermal_input import make_deck
from open_cae_evidence.thermal_output import parse_dat, verify_frd
from open_cae_evidence.thermal_verification import verify

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/thermal-ccx-2.21"


class CapturedThermalFormatTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads(
            (ROOT / "benchmarks/slab_steady_heat.json").read_text()
        )
        self.deck, self.nodes, self.elements, _ = make_deck(self.config, [2, 2, 2])
        self.dat = (FIXTURE / "slab.dat").read_text()

    def test_unmodified_real_221_fixture_passes_all_physics(self):
        metadata = json.loads((FIXTURE / "provenance.json").read_text())
        for name, digest in metadata["sha256"].items():
            self.assertEqual(
                hashlib.sha256((FIXTURE / name).read_bytes()).hexdigest(), digest
            )
        self.assertEqual((FIXTURE / "slab.inp").read_text(), self.deck)
        data = parse_dat(self.dat, self.nodes, self.elements)
        self.assertTrue(verify(self.config, self.nodes, self.elements, data)["passed"])
        self.assertTrue(
            verify_frd(
                (FIXTURE / "slab.frd").read_text(),
                self.nodes,
                self.elements,
                data,
                "2.21",
            )["passed"]
        )

    def test_wrong_or_repeated_step_increment_metadata_fails(self):
        for bad in (
            self.dat.replace("S T E P       1", "S T E P       2", 1),
            self.dat.replace("INCREMENT     1", "INCREMENT     2", 1),
            "S T E P 1\nINCREMENT 1\n" + self.dat,
        ):
            with self.subTest(length=len(bad)), self.assertRaises(ValueError):
                parse_dat(bad, self.nodes, self.elements)
