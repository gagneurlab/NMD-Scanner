"""
Run every conformance case on both strands.
Check that the drawing of each case holds its rendered layout block.
"""

import importlib
import textwrap
from pathlib import Path

import pytest

from .runner import case_params, check, render_case

CASE_MODULES = {
    path.stem: importlib.import_module(f".{path.stem}", __package__)
    for path in sorted(Path(__file__).parent.glob("cases_*.py"))
}
CASES = [case for module in CASE_MODULES.values() for case in module.CASES]
PARAMS = case_params(CASES)


@pytest.mark.parametrize(("case", "change", "strand"), PARAMS)
def test_case(case, change, strand, tmp_path):
    check(case, change, strand, tmp_path)


@pytest.mark.parametrize("case", CASES, ids=[case.name for case in CASES])
def test_drawing_holds_its_layout_block(case):
    block = render_case(case)
    lines = [line.rstrip() for line in textwrap.dedent(case.drawing).splitlines()]
    wanted = block.splitlines()
    held = any(lines[i : i + len(wanted)] == wanted for i in range(len(lines)))
    assert held, f"the drawing does not hold its layout block, at the indentation of its other lines:\n{block}"


def test_case_names_are_unique():
    names = [case.name for case in CASES]
    assert sorted(name for name in set(names) if names.count(name) > 1) == []
