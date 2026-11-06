"""A 'port' that only raises, so the runner's traceback rule is exercised for real (stdlib only)."""


def boom(message: str = "fakeport failure"):
    raise RuntimeError(message)
