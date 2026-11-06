"""Source-named tests for differential/conftest.py's strict red-team marks."""

from __future__ import annotations

import runpy
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

collect = runpy.run_path(str(Path(__file__).parents[1] / "differential" / "conftest.py"))[
    "pytest_collection_modifyitems"
]


@dataclass
class Item:
    callspec: SimpleNamespace | None = None
    markers: list[pytest.MarkDecorator] = field(default_factory=list)

    def add_marker(self, marker: pytest.MarkDecorator) -> None:
        self.markers.append(marker)


def case_item(spec: dict[str, Any]) -> Item:
    return Item(SimpleNamespace(params={"case": SimpleNamespace(id="example", spec=spec)}))


def test_redteam_mark_is_strict_and_keeps_reason() -> None:
    item = case_item({"xfail": "Known circular overlap bug"})
    collect([item])
    assert len(item.markers) == 1
    assert item.markers[0].name == "xfail"
    assert item.markers[0].kwargs == {"reason": "Known circular overlap bug", "strict": True}


def test_redteam_mark_leaves_existing_cases_unchanged() -> None:
    items = [case_item({}), Item(), Item(SimpleNamespace(params={"other": 1}))]
    collect(items)
    assert all(not item.markers for item in items)


@pytest.mark.parametrize("reason", ["", "  ", None, False])
def test_redteam_mark_requires_reason(reason: Any) -> None:
    with pytest.raises(pytest.UsageError, match="example: xfail must be a nonempty reason"):
        collect([case_item({"xfail": reason})])
