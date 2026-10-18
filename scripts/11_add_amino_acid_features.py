#!/usr/bin/env python3
"""Calculate deterministic Grantham and Kyte-Doolittle substitution descriptors."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import (
    PROJECT_ROOT,
    append_stage_manifest,
    audit_basic_rows,
    parse_int,
    read_csv,
    write_csv,
    write_json,
)


# Grantham (1974): composition, polarity and molecular volume.
# Kyte and Doolittle (1982): hydropathy.
AA = {
    "A": (0.00, 8.1, 31.0, 1.8), "R": (0.65, 10.5, 124.0, -4.5),
    "N": (1.33, 11.6, 56.0, -3.5), "D": (1.38, 13.0, 54.0, -3.5),
    "C": (2.75, 5.5, 55.0, 2.5), "Q": (0.89, 10.5, 85.0, -3.5),
    "E": (0.92, 12.3, 83.0, -3.5), "G": (0.74, 9.0, 3.0, -0.4),
    "H": (0.58, 10.4, 96.0, -3.2), "I": (0.00, 5.2, 111.0, 4.5),
    "L": (0.00, 4.9, 111.0, 3.8), "K": (0.33, 11.3, 119.0, -3.9),
    "M": (0.00, 5.7, 105.0, 1.9), "F": (0.00, 5.2, 132.0, 2.8),
    "P": (0.39, 8.0, 32.5, -1.6), "S": (1.42, 9.2, 32.0, -0.8),
    "T": (0.71, 8.6, 61.0, -0.7), "W": (0.13, 5.4, 170.0, -0.9),
    "Y": (0.20, 6.2, 136.0, -1.3), "V": (0.00, 5.9, 84.0, 4.2),
}
NONPOLAR = set("AVLIMFWPG")
CHARGE = {**{code: 0.0 for code in AA}, "K": 1.0, "R": 1.0, "D": -1.0, "E": -1.0, "H": 0.1}


ADDED_COLUMNS = [
    "aa_feature_status",
    "aa_grantham_score",
    "aa_composition_delta",
    "aa_polarity_delta",
    "aa_volume_delta",
    "aa_hydropathy_delta",
    "aa_charge_delta",
    "aa_charge_changed",
    "aa_polarity_class_changed",
    "aa_relative_protein_position",
]


def grantham(reference: str, alternate: str) -> float:
    c1, p1, v1, _ = AA[reference]
    c2, p2, v2, _ = AA[alternate]
    return 50.723 * math.sqrt(
        1.833 * (c1 - c2) ** 2 + 0.1018 * (p1 - p2) ** 2 + 0.000399 * (v1 - v2) ** 2
    )


def values(row: dict[str, str]) -> dict[str, Any]:
    result = {column: "" for column in ADDED_COLUMNS}
    reference = row.get("reconciled_reference_amino_acid", "").upper()
    alternate = row.get("reconciled_alternate_amino_acid", "").upper()
    if reference not in AA or alternate not in AA:
        result["aa_feature_status"] = "INVALID_OR_NONSTANDARD_AMINO_ACID"
        return result
    c1, p1, v1, h1 = AA[reference]
    c2, p2, v2, h2 = AA[alternate]
    position = parse_int(row.get("reconciled_protein_position"))
    length = parse_int(row.get("uniprot_protein_length"))
    relative = position / length if position and length and position <= length else None
    result.update(
        {
            "aa_feature_status": "COMPLETE" if relative is not None else "RELATIVE_POSITION_MISSING",
            "aa_grantham_score": round(grantham(reference, alternate), 6),
            "aa_composition_delta": round(c2 - c1, 6),
            "aa_polarity_delta": round(p2 - p1, 6),
            "aa_volume_delta": round(v2 - v1, 6),
            "aa_hydropathy_delta": round(h2 - h1, 6),
            "aa_charge_delta": round(CHARGE[alternate] - CHARGE[reference], 6),
            "aa_charge_changed": int(CHARGE[alternate] != CHARGE[reference]),
            "aa_polarity_class_changed": int((reference in NONPOLAR) != (alternate in NONPOLAR)),
            "aa_relative_protein_position": "" if relative is None else round(relative, 8),
        }
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(PROJECT_ROOT / "data/interim/10_dbnsfp.csv"))
    parser.add_argument("--output", default=str(PROJECT_ROOT / "data/interim/11_amino_acid_features.csv"))
    args = parser.parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    rows = read_csv(input_path)
    enriched = []
    for row in rows:
        combined = dict(row)
        combined.update(values(row))
        enriched.append(combined)
    fieldnames = list(rows[0]) + [column for column in ADDED_COLUMNS if column not in rows[0]]
    write_csv(output_path, enriched, fieldnames)
    audit = {
        **audit_basic_rows(enriched),
        "feature_status_counts": dict(Counter(row["aa_feature_status"] for row in enriched)),
        "sources": {
            "composition_polarity_volume_and_distance": "Grantham (1974)",
            "hydropathy": "Kyte and Doolittle (1982)",
        },
        "note": "These are deterministic physicochemical descriptors, not clinical evidence by themselves.",
    }
    audit_path = PROJECT_ROOT / "audits/11_amino_acid_features_audit.json"
    write_json(audit_path, audit)
    append_stage_manifest("11_add_amino_acid_features", [input_path], [output_path, audit_path], audit["sources"])
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
