import json
from pathlib import Path

import pytest

from posegate.exceptions import ImmutableRecordError
from posegate.records import capture_path, write_immutable_json


def test_capture_path_is_deterministic_per_run_and_policy(tmp_path: Path) -> None:
    first = capture_path(tmp_path, "protein|pocket 1", "checkpoint/5ns")
    second = capture_path(tmp_path, "protein|pocket 1", "checkpoint/5ns")
    assert first == second
    assert first.name == "protein-pocket-1__checkpoint-5ns.json"


def test_record_cannot_be_overwritten(tmp_path: Path) -> None:
    path = tmp_path / "record.json"
    write_immutable_json(path, {"record_id": "first"})
    assert json.loads(path.read_text())["record_id"] == "first"
    with pytest.raises(ImmutableRecordError, match="already exists"):
        write_immutable_json(path, {"record_id": "second"})
    assert json.loads(path.read_text())["record_id"] == "first"
