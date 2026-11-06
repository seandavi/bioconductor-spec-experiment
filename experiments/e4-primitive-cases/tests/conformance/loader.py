"""Find case files and their R-generated goldens."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[2]
CASES_DIR = ROOT / "cases"
GOLDENS_DIR = ROOT / "goldens"
ADR_DIR = ROOT / "adr"


@dataclass(frozen=True)
class Case:
    id: str
    spec: dict[str, Any]
    golden: dict[str, Any] | None
    source: Path


def load_cases() -> list[Case]:
    out: list[Case] = []
    for f in sorted(CASES_DIR.glob("*/*.yaml")):
        gold_file = GOLDENS_DIR / f.parent.name / f"{f.stem}.json"
        goldens = json.loads(gold_file.read_text())["cases"] if gold_file.exists() else {}
        for spec in yaml.safe_load(f.read_text())["cases"]:
            out.append(Case(spec["id"], spec, goldens.get(spec["id"]), f))
    return out


def golden_ids() -> set[str]:
    ids: set[str] = set()
    for f in GOLDENS_DIR.glob("*/*.json"):
        ids |= set(json.loads(f.read_text())["cases"])
    return ids
