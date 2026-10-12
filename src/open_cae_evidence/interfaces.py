"""Minimal in-process interfaces; JSON contracts are the interoperability boundary."""

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol


@dataclass(frozen=True)
class ProblemSpec:
    id: str
    analysis: str
    parameters: dict
    units: dict
    boundary_conditions: list
    required_quantities: tuple[str, ...]
    unsupported_capabilities: tuple[str, ...]


@dataclass(frozen=True)
class PreparedInput:
    problem: ProblemSpec
    counts: list[int]
    files: dict[str, str]
    stem: str
    nodes: dict
    elements: dict
    audit: dict
    topology: str
    integration: str
    auxiliary_nodes: frozenset[int] = frozenset()


class SolverAdapter(Protocol):
    name: str
    version: str

    def prepare(self, problem: ProblemSpec, counts: list[int]) -> PreparedInput: ...
    def execute(
        self, prepared: PreparedInput, directory: Path, executable: str
    ) -> str: ...
    def extract(
        self, prepared: PreparedInput, directory: Path, solver_version: str
    ) -> dict: ...


class NumericalVerifier(Protocol):
    tolerances: dict[str, float]

    def verify(self, problem: ProblemSpec, normalized: dict) -> dict: ...
    def study(self, runs: list[dict]) -> dict | list: ...
