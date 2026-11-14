#!/usr/bin/env python3
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.modeling.evidence_report import run_evidence_reporting

if __name__ == "__main__":
    print(json.dumps(run_evidence_reporting(ROOT), indent=2, sort_keys=True))
