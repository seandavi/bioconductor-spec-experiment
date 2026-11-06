"""Self-test adapter: returns R's own projection (from expected.jsonl). Must score 100 % pass.

An R error becomes a port error raised inside `fakeport`, so it is classified `port_error` (R1 -> pass).
"""

import json
from pathlib import Path

from adapter_base import AdapterBase

import fakeport


class Adapter(AdapterBase):
    port = "r_echo"
    port_packages = ("fakeport",)

    def setup(self, ctx: dict) -> None:
        self.expected = {}
        for line in Path(ctx["expected_path"]).read_text().splitlines():
            rec = json.loads(line)
            if rec["status"] == "ok":
                self.expected[rec["id"]] = rec

    def run_case(self, case: dict, stage) -> dict:
        exp = self.expected[case["id"]]
        stage("build", lambda: None)
        if exp["result"]["kind"] == "error":
            stage("call", fakeport.boom, "R error: " + exp["result"]["message"][:100])
        return {"result": stage("call", lambda: json.loads(json.dumps(exp["result"])))}
