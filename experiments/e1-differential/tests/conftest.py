import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--refresh-r", action="store_true", default=False,
        help="re-run every differential case in R and rewrite r_cache/",
    )
