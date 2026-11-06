import sys
from pathlib import Path

EVAL = Path(__file__).resolve().parents[1]
for sub in ("runner", "adapters", "selftest"):
    sys.path.insert(0, str(EVAL / sub))
