"""Keep unresolved red-team cases executable without accepting their differences."""

from __future__ import annotations

import pytest


def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    for item in items:
        callspec = getattr(item, "callspec", None)
        case = callspec.params.get("case") if callspec is not None else None
        if case is not None and "xfail" in case.spec:
            reason = case.spec["xfail"]
            if not isinstance(reason, str) or not reason.strip():
                raise pytest.UsageError(f"{case.id}: xfail must be a nonempty reason")
            item.add_marker(pytest.mark.xfail(reason=reason, strict=True))
