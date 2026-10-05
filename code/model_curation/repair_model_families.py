#!/usr/bin/env python3
"""Repair and package the four publication model families.

The script never overwrites source models.  It reconstructs stale annotation
vectors by identifier from the older, internally aligned model generation,
uses the local KEGG snapshot only for identifiers absent from that generation,
rebuilds ``rules`` and ``rxnGeneMat`` from human-readable ``grRules``, and
removes genes that are unused by every GPR. The source archive is read-only;
all repaired models are written to a separate publication repository.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat, savemat
from scipy.sparse import csc_matrix, issparse


METABOLITE_FIELDS = (
    "metAbbreviation",
    "metAliases",
    "metBiGGID",
    "metInChIKey",
    "metMetaNetXID",
    "metMiriams",
    "metPka",
    "metPkb",
    "metSBOTerms",
    "metSmiles",
)

REACTION_FIELDS = (
    "rxnKEGGID",
    "rxnMiriams",
    "rxnRheaID",
    "rxnSBOTerms",
)

GENE_FIELDS = ("genes", "geneSBOTerms", "EnzymeDescription")

TARGET_FIELDS = set(METABOLITE_FIELDS) | set(REACTION_FIELDS) | {
    "genes",
    "geneSBOTerms",
    "EnzymeDescription",
    "rules",
    "rxnGeneMat",
}

GPR_TOKEN_PATTERN = re.compile(r"\(|\)|\band\b|\bor\b|[^\s()]+")
COMPARTMENT_PATTERN = re.compile(r"_[cep]$")
KEGG_REACTION_PATTERN = re.compile(r"^R\d{5}$")
KEGG_COMPOUND_PATTERN = re.compile(r"^[CGD]\d{5}$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    script_dir = Path(__file__).resolve().parent
    repository_root = script_dir.parents[1]
    parser.add_argument(
        "--family-manifest",
        type=Path,
        default=script_dir / "model_families.tsv",
        help="TSV describing source, donor, and output paths for each family",
    )
    parser.add_argument(
        "--code-root",
        type=Path,
        required=True,
        help="Path to the untouched SoilSynCom Code directory",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=repository_root / "models" / "final",
        help="Destination root for repaired model families",
    )
    parser.add_argument(
        "--model-manifest",
        type=Path,
        default=repository_root / "models" / "model_manifest.tsv",
    )
    parser.add_argument(
        "--removed-genes",
        type=Path,
        default=repository_root / "models" / "removed_unused_genes.tsv",
    )
    parser.add_argument(
        "--validation-report",
        type=Path,
        default=repository_root / "models" / "VALIDATION.md",
    )
    parser.add_argument(
        "--validation-json",
        type=Path,
        default=repository_root / "models" / "validation_summary.json",
    )
    parser.add_argument(
        "--kegg-mets",
        type=Path,
        help="Defaults to DraftModelReconstruction/keggdb/keggMets.mat",
    )
    parser.add_argument(
        "--kegg-rxns",
        type=Path,
        help="Defaults to DraftModelReconstruction/keggdb/keggRxns.mat",
    )
    parser.add_argument(
        "--strain-crosswalk",
        type=Path,
        help="Optional strain_identity_crosswalk.tsv with GTDB assignments",
    )
    parser.add_argument(
        "--existing-rich-tfa",
        type=Path,
        help="Optional earlier validated repaired rich/TFA directory for semantic comparison",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_model(path: Path) -> tuple[str, dict[str, Any]]:
    payload = loadmat(path, simplify_cells=True)
    candidates = [
        (key, value)
        for key, value in payload.items()
        if not key.startswith("__")
        and isinstance(value, dict)
        and "rxns" in value
        and "mets" in value
    ]
    if len(candidates) != 1:
        raise ValueError(f"{path}: expected one model struct, found {len(candidates)}")
    return candidates[0]


def values(value: Any) -> list[Any]:
    if isinstance(value, list):
        return list(value)
    if isinstance(value, np.ndarray):
        if value.ndim == 0:
            return [value.item()]
        return list(value.reshape(-1))
    return [value]


def text(value: Any) -> str:
    if isinstance(value, np.ndarray) and value.size == 0:
        return ""
    if value is None:
        return ""
    result = str(value).strip()
    return "" if result in {"[]", "None", "nan"} else result


def object_vector(items: list[Any]) -> np.ndarray:
    result = np.empty(len(items), dtype=object)
    for index, item in enumerate(items):
        result[index] = copy.deepcopy(item)
    return result


def empty_cell() -> np.ndarray:
    return np.empty((0, 0), dtype=float)


def first_model(path: Path) -> dict[str, Any]:
    payload = loadmat(path, simplify_cells=True)
    models = [
        value
        for key, value in payload.items()
        if not key.startswith("__") and isinstance(value, dict)
    ]
    if len(models) != 1:
        raise ValueError(f"{path}: expected one reference model struct")
    return models[0]


def canonical_metabolite(identifier: str) -> str:
    return COMPARTMENT_PATTERN.sub("", identifier)


def preferred_metabolite_indices(identifiers: list[str]) -> dict[str, int]:
    rank = {"_c": 0, "_e": 1, "_p": 2}
    result: dict[str, tuple[int, int]] = {}
    for index, identifier in enumerate(identifiers):
        canonical = canonical_metabolite(identifier)
        suffix = identifier[-2:] if len(identifier) >= 2 else ""
        candidate = (rank.get(suffix, 3), index)
        if canonical not in result or candidate < result[canonical]:
            result[canonical] = candidate
    return {key: candidate[1] for key, candidate in result.items()}


def build_reference_maps(
    kegg_mets_path: Path, kegg_rxns_path: Path
) -> tuple[dict[str, Any], dict[str, Any]]:
    kegg_mets = first_model(kegg_mets_path)
    metabolite_miriams = {
        text(identifier): copy.deepcopy(miriam)
        for identifier, miriam in zip(
            values(kegg_mets["mets"]), values(kegg_mets["metMiriams"])
        )
        if text(identifier)
    }

    kegg_rxns = first_model(kegg_rxns_path)
    reaction_miriams = {
        text(identifier): copy.deepcopy(miriam)
        for identifier, miriam in zip(
            values(kegg_rxns["rxns"]), values(kegg_rxns["rxnMiriams"])
        )
        if text(identifier)
    }
    return metabolite_miriams, reaction_miriams


def miriam_rhea(miriam: Any) -> Any:
    if not isinstance(miriam, dict):
        return empty_cell()
    names = [text(item) for item in values(miriam.get("name", []))]
    vals = values(miriam.get("value", []))
    matches = [copy.deepcopy(vals[index]) for index, name in enumerate(names) if name == "rhea"]
    if not matches:
        return empty_cell()
    return matches[0] if len(matches) == 1 else object_vector(matches)


def minimal_metabolite_miriam(kegg_id: str) -> dict[str, np.ndarray]:
    return {
        "name": object_vector(["kegg.compound"]),
        "value": object_vector([kegg_id]),
    }


def minimal_reaction_miriam(kegg_id: str) -> dict[str, np.ndarray]:
    return {
        "name": object_vector(["kegg.reaction"]),
        "value": object_vector([kegg_id]),
    }


def repair_metabolite_fields(
    model: dict[str, Any],
    donor: dict[str, Any],
    kegg_miriams: dict[str, Any],
) -> dict[str, int]:
    identifiers = [text(item) for item in values(model["mets"])]
    donor_identifiers = [text(item) for item in values(donor["mets"])]
    donor_exact = {identifier: index for index, identifier in enumerate(donor_identifiers)}
    donor_canonical = preferred_metabolite_indices(donor_identifiers)

    donor_kegg_ids = [text(item) for item in values(donor["metKEGGID"])]
    donor_by_kegg: dict[str, int] = {}
    for index, kegg_id in enumerate(donor_kegg_ids):
        if not kegg_id:
            continue
        current = donor_by_kegg.get(kegg_id)
        if current is None:
            donor_by_kegg[kegg_id] = index
            continue
        old_id = donor_identifiers[current]
        new_id = donor_identifiers[index]
        if new_id.endswith("_c") and not old_id.endswith("_c"):
            donor_by_kegg[kegg_id] = index

    current_kegg_ids = [text(item) for item in values(model["metKEGGID"])]
    if len(current_kegg_ids) != len(identifiers):
        raise ValueError("metKEGGID is not aligned with mets")

    donor_fields = {field: values(donor[field]) for field in METABOLITE_FIELDS}
    for field, field_values in donor_fields.items():
        if len(field_values) != len(donor_identifiers):
            raise ValueError(f"donor field {field} is not aligned with donor mets")

    repaired = {field: [] for field in METABOLITE_FIELDS}
    counts = {"exact": 0, "kegg_donor": 0, "canonical": 0, "kegg_database": 0, "minimal": 0}

    for identifier, kegg_id in zip(identifiers, current_kegg_ids):
        canonical_identifier = canonical_metabolite(identifier)
        effective_kegg_id = kegg_id
        if not effective_kegg_id and KEGG_COMPOUND_PATTERN.fullmatch(
            canonical_identifier
        ):
            effective_kegg_id = canonical_identifier
        donor_index: int | None = None
        route = "minimal"
        if identifier in donor_exact:
            donor_index = donor_exact[identifier]
            route = "exact"
        elif effective_kegg_id and effective_kegg_id in donor_by_kegg:
            donor_index = donor_by_kegg[effective_kegg_id]
            route = "kegg_donor"
        elif canonical_identifier in donor_canonical:
            donor_index = donor_canonical[canonical_identifier]
            route = "canonical"

        if donor_index is not None:
            for field in METABOLITE_FIELDS:
                repaired[field].append(copy.deepcopy(donor_fields[field][donor_index]))
            counts[route] += 1
            continue

        for field in METABOLITE_FIELDS:
            repaired[field].append(empty_cell())

        if effective_kegg_id:
            if effective_kegg_id in kegg_miriams:
                repaired["metMiriams"][-1] = copy.deepcopy(
                    kegg_miriams[effective_kegg_id]
                )
                counts["kegg_database"] += 1
            else:
                repaired["metMiriams"][-1] = minimal_metabolite_miriam(
                    effective_kegg_id
                )
                counts["minimal"] += 1
            repaired["metSBOTerms"][-1] = "SBO:0000247"
        else:
            counts["minimal"] += 1

    for field in METABOLITE_FIELDS:
        model[field] = object_vector(repaired[field])
    return counts


def repair_reaction_fields(
    model: dict[str, Any],
    donor: dict[str, Any],
    kegg_miriams: dict[str, Any],
) -> dict[str, int]:
    identifiers = [text(item) for item in values(model["rxns"])]
    donor_identifiers = [text(item) for item in values(donor["rxns"])]
    donor_exact = {identifier: index for index, identifier in enumerate(donor_identifiers)}
    donor_fields = {field: values(donor[field]) for field in REACTION_FIELDS}
    for field, field_values in donor_fields.items():
        if len(field_values) != len(donor_identifiers):
            raise ValueError(f"donor field {field} is not aligned with donor rxns")

    repaired = {field: [] for field in REACTION_FIELDS}
    counts = {"exact": 0, "kegg_database": 0, "custom_empty": 0}
    for identifier in identifiers:
        if identifier in donor_exact:
            donor_index = donor_exact[identifier]
            for field in REACTION_FIELDS:
                repaired[field].append(copy.deepcopy(donor_fields[field][donor_index]))
            counts["exact"] += 1
            continue

        for field in REACTION_FIELDS:
            repaired[field].append(empty_cell())
        if KEGG_REACTION_PATTERN.fullmatch(identifier):
            miriam = copy.deepcopy(
                kegg_miriams.get(identifier, minimal_reaction_miriam(identifier))
            )
            repaired["rxnKEGGID"][-1] = identifier
            repaired["rxnMiriams"][-1] = miriam
            repaired["rxnRheaID"][-1] = miriam_rhea(miriam)
            repaired["rxnSBOTerms"][-1] = "SBO:0000176"
            counts["kegg_database"] += 1
        else:
            counts["custom_empty"] += 1

    for field in REACTION_FIELDS:
        model[field] = object_vector(repaired[field])
    return counts


def expression_tokens(expression: Any) -> list[str]:
    if isinstance(expression, np.ndarray) and expression.size == 0:
        return []
    value = text(expression)
    if not value:
        return []
    return GPR_TOKEN_PATTERN.findall(value)


def gene_tokens(expression: Any) -> list[str]:
    return [
        token
        for token in expression_tokens(expression)
        if token not in {"(", ")", "and", "or"}
    ]


def build_rule(expression: Any, gene_index: dict[str, int]) -> str:
    converted: list[str] = []
    for token in expression_tokens(expression):
        if token == "and":
            converted.append("&")
        elif token == "or":
            converted.append("|")
        elif token in {"(", ")"}:
            converted.append(token)
        else:
            if token not in gene_index:
                raise ValueError(f"GPR references undeclared gene {token!r}")
            converted.append(f"x({gene_index[token] + 1})")
    rule = " ".join(converted)
    return rule.replace("( ", "(").replace(" )", ")")


def repair_genes(model: dict[str, Any]) -> tuple[list[str], int, int]:
    genes = [text(item) for item in values(model["genes"])]
    original_gene_index = {gene: index for index, gene in enumerate(genes)}
    if len(original_gene_index) != len(genes):
        raise ValueError("duplicate gene identifiers prevent safe repair")

    gr_rules = values(model["grRules"])
    rxns = values(model["rxns"])
    if len(gr_rules) != len(rxns):
        raise ValueError("grRules is not aligned with rxns")

    used: set[str] = set()
    for expression in gr_rules:
        for gene in gene_tokens(expression):
            if gene not in original_gene_index:
                raise ValueError(f"GPR references undeclared gene {gene!r}")
            used.add(gene)

    keep_indices = [index for index, gene in enumerate(genes) if gene in used]
    removed = [gene for gene in genes if gene not in used]
    kept_genes = [genes[index] for index in keep_indices]
    new_gene_index = {gene: index for index, gene in enumerate(kept_genes)}

    trimmed_gene_metadata_entries = 0
    for field in GENE_FIELDS:
        field_values = values(model[field])
        if len(field_values) != len(genes):
            trailing = field_values[len(genes) :]
            safe_trailing_duplicate = (
                field == "EnzymeDescription"
                and len(field_values) > len(genes)
                and all(text(item) in original_gene_index for item in trailing)
            )
            if not safe_trailing_duplicate:
                raise ValueError(f"gene field {field} is not aligned with genes")
            trimmed_gene_metadata_entries += len(trailing)
            field_values = field_values[: len(genes)]
        model[field] = object_vector([field_values[index] for index in keep_indices])

    row_indices: list[int] = []
    column_indices: list[int] = []
    rules: list[str] = []
    for row, expression in enumerate(gr_rules):
        tokens = gene_tokens(expression)
        for gene in sorted(set(tokens), key=lambda item: new_gene_index[item]):
            row_indices.append(row)
            column_indices.append(new_gene_index[gene])
        rules.append(build_rule(expression, new_gene_index))

    original_dtype = getattr(model["rxnGeneMat"], "dtype", np.int8)
    data = np.ones(len(row_indices), dtype=original_dtype)
    model["rxnGeneMat"] = csc_matrix(
        (data, (row_indices, column_indices)),
        shape=(len(rxns), len(kept_genes)),
    )
    model["rules"] = object_vector(rules)
    return removed, len(used), trimmed_gene_metadata_entries


def equal_value(left: Any, right: Any) -> bool:
    """Recursively compare values after a scipy MAT round-trip."""
    if issparse(left) or issparse(right):
        if not (issparse(left) and issparse(right)) or left.shape != right.shape:
            return False
        difference = left != right
        return difference.nnz == 0
    if isinstance(left, dict) or isinstance(right, dict):
        if not (isinstance(left, dict) and isinstance(right, dict)):
            return False
        return left.keys() == right.keys() and all(
            equal_value(left[key], right[key]) for key in left
        )
    if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
        left_array = np.asarray(left)
        right_array = np.asarray(right)
        if left_array.shape != right_array.shape:
            return False
        if left_array.dtype == object or right_array.dtype == object:
            return all(
                equal_value(left_item, right_item)
                for left_item, right_item in zip(
                    left_array.reshape(-1), right_array.reshape(-1)
                )
            )
        try:
            return bool(np.array_equal(left_array, right_array, equal_nan=True))
        except TypeError:
            return bool(np.array_equal(left_array, right_array))
    if isinstance(left, (list, tuple)) or isinstance(right, (list, tuple)):
        if not (isinstance(left, (list, tuple)) and isinstance(right, (list, tuple))):
            return False
        return len(left) == len(right) and all(
            equal_value(left_item, right_item)
            for left_item, right_item in zip(left, right)
        )
    if isinstance(left, (float, np.floating)) and isinstance(
        right, (float, np.floating)
    ):
        if math.isnan(float(left)) and math.isnan(float(right)):
            return True
    try:
        result = left == right
        return bool(result) if np.isscalar(result) else bool(np.all(result))
    except (TypeError, ValueError):
        return repr(left) == repr(right)


def validate_repair(
    original: dict[str, Any], repaired: dict[str, Any], variable_name: str
) -> dict[str, int]:
    n_rxns = len(values(repaired["rxns"]))
    n_mets = len(values(repaired["mets"]))
    n_genes = len(values(repaired["genes"]))

    for field in METABOLITE_FIELDS:
        if len(values(repaired[field])) != n_mets:
            raise ValueError(f"{variable_name}: repaired {field} is not aligned")
    for field in REACTION_FIELDS:
        if len(values(repaired[field])) != n_rxns:
            raise ValueError(f"{variable_name}: repaired {field} is not aligned")
    for field in GENE_FIELDS:
        if len(values(repaired[field])) != n_genes:
            raise ValueError(f"{variable_name}: repaired {field} is not aligned")
    if repaired["rxnGeneMat"].shape != (n_rxns, n_genes):
        raise ValueError(f"{variable_name}: rxnGeneMat has the wrong shape")

    matrix = repaired["rxnGeneMat"].tocsr()
    if np.any(np.asarray(matrix.getnnz(axis=0)).reshape(-1) == 0):
        raise ValueError(f"{variable_name}: orphan genes remain")

    rules = values(repaired["rules"])
    for row, rule in enumerate(rules):
        rule_indices = {
            int(index) - 1 for index in re.findall(r"x\((\d+)\)", text(rule))
        }
        matrix_indices = set(matrix.indices[matrix.indptr[row] : matrix.indptr[row + 1]])
        if rule_indices != matrix_indices:
            raise ValueError(f"{variable_name}: rule/matrix mismatch at row {row}")
        if any(index < 0 or index >= n_genes for index in rule_indices):
            raise ValueError(f"{variable_name}: rule contains an out-of-range gene")

    if original.keys() != repaired.keys():
        added = sorted(repaired.keys() - original.keys())
        removed = sorted(original.keys() - repaired.keys())
        raise ValueError(
            f"{variable_name}: model fields changed; added={added}, removed={removed}"
        )
    for field in original:
        if field not in TARGET_FIELDS and not equal_value(original[field], repaired[field]):
            raise ValueError(f"{variable_name}: non-target field {field} changed")

    return {"rxns": n_rxns, "mets": n_mets, "genes": n_genes}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def strain_from_variable(variable_name: str) -> str:
    return variable_name.removeprefix("ttmodel_")


def model_files(path: Path) -> list[Path]:
    return [
        item
        for item in sorted(path.glob("*.mat"))
        if item.name != "AllWorkspace.mat"
    ]


def index_donors(path: Path) -> dict[str, tuple[Path, str, dict[str, Any]]]:
    result: dict[str, tuple[Path, str, dict[str, Any]]] = {}
    for donor_path in model_files(path):
        variable_name, donor = load_model(donor_path)
        strain = strain_from_variable(variable_name)
        if strain in result:
            raise ValueError(f"{path}: duplicate donor for strain {strain}")
        result[strain] = (donor_path, variable_name, donor)
    return result


def load_crosswalk(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None or not path.is_file():
        return {}
    rows = read_tsv(path)
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        strain = row.get("model_code", "").strip()
        if strain:
            result[strain] = row
    return result


def orphan_gene_names(model: dict[str, Any]) -> list[str]:
    genes = [text(item) for item in values(model["genes"])]
    matrix = model["rxnGeneMat"]
    if not issparse(matrix):
        matrix = csc_matrix(matrix)
    nonzero = np.asarray(matrix.tocsc().getnnz(axis=0)).reshape(-1)
    if len(nonzero) != len(genes):
        raise ValueError("rxnGeneMat gene dimension does not match genes")
    return [gene for gene, count in zip(genes, nonzero) if count == 0]


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    code_root = args.code_root.resolve()
    output_root = args.output_root.resolve()
    family_manifest_path = args.family_manifest.resolve()
    model_manifest_path = args.model_manifest.resolve()
    removed_genes_path = args.removed_genes.resolve()
    validation_report_path = args.validation_report.resolve()
    validation_json_path = args.validation_json.resolve()
    kegg_mets_path = (
        args.kegg_mets.resolve()
        if args.kegg_mets
        else code_root / "DraftModelReconstruction" / "keggdb" / "keggMets.mat"
    )
    kegg_rxns_path = (
        args.kegg_rxns.resolve()
        if args.kegg_rxns
        else code_root / "DraftModelReconstruction" / "keggdb" / "keggRxns.mat"
    )
    default_crosswalk = (
        code_root.parent
        / "SoilSynCom_supmat"
        / "V4_assets"
        / "source_data"
        / "strain_identity_crosswalk.tsv"
    )
    crosswalk_path = (
        args.strain_crosswalk.resolve()
        if args.strain_crosswalk
        else default_crosswalk
    )
    default_existing_rich_tfa = (
        code_root.parent
        / "Model_Curation_Python"
        / "FinalGEMsCuratedExchanges_PythonFixed"
    )
    existing_rich_tfa = (
        args.existing_rich_tfa.resolve()
        if args.existing_rich_tfa
        else default_existing_rich_tfa
    )

    if not code_root.is_dir():
        raise FileNotFoundError(code_root)
    if code_root == output_root or code_root in output_root.parents:
        raise ValueError("output root must not be inside the untouched Code archive")
    for destination in (
        model_manifest_path,
        removed_genes_path,
        validation_report_path,
        validation_json_path,
    ):
        if destination.exists():
            raise FileExistsError(f"refusing to overwrite existing file: {destination}")

    family_rows = read_tsv(family_manifest_path)
    required_family_columns = {
        "family_id",
        "medium",
        "carbon_source",
        "formulation",
        "model_type",
        "source_relative",
        "donor_relative",
        "output_relative",
        "expected_count",
        "expected_reconnected_gpr_genes",
    }
    if len(family_rows) != 4:
        raise ValueError(f"expected four model families, found {len(family_rows)}")
    if not required_family_columns.issubset(family_rows[0]):
        missing = sorted(required_family_columns - family_rows[0].keys())
        raise ValueError(f"family manifest is missing columns: {missing}")
    expected_total = sum(int(row["expected_count"]) for row in family_rows)
    if expected_total != 80:
        raise ValueError(f"family manifest must describe exactly 80 models, found {expected_total}")

    for row in family_rows:
        family_output = output_root / row["output_relative"]
        if family_output.exists():
            raise FileExistsError(
                f"refusing to overwrite existing family directory: {family_output}"
            )

    kegg_met_miriams, kegg_rxn_miriams = build_reference_maps(
        kegg_mets_path, kegg_rxns_path
    )
    crosswalk = load_crosswalk(crosswalk_path)
    manifest_rows: list[dict[str, Any]] = []
    removed_gene_rows: list[dict[str, Any]] = []
    family_summaries: list[dict[str, Any]] = []
    model_validation_rows: list[dict[str, Any]] = []
    rich_tfa_semantic_matches = 0

    for family in family_rows:
        family_id = family["family_id"]
        source_dir = code_root / family["source_relative"]
        donor_dir = code_root / family["donor_relative"]
        output_dir = output_root / family["output_relative"]
        expected_count = int(family["expected_count"])
        expected_reconnected = int(family["expected_reconnected_gpr_genes"])
        if not source_dir.is_dir() or not donor_dir.is_dir():
            raise FileNotFoundError(
                f"{family_id}: missing source or donor directory: {source_dir}, {donor_dir}"
            )
        source_files = model_files(source_dir)
        if len(source_files) != expected_count:
            raise ValueError(
                f"{family_id}: expected {expected_count} source models, found {len(source_files)}"
            )
        donor_index = index_donors(donor_dir)
        output_dir.mkdir(parents=True, exist_ok=False)
        family_reconnected: set[tuple[str, str]] = set()
        family_removed = 0
        family_reaction_rebuilds = 0

        for source_path in source_files:
            variable_name, original = load_model(source_path)
            strain = strain_from_variable(variable_name)
            if strain not in donor_index:
                raise FileNotFoundError(
                    f"{family_id}: no donor model found for strain {strain}"
                )
            donor_path, donor_variable, donor = donor_index[strain]
            if donor_variable != variable_name:
                raise ValueError(
                    f"{family_id}: model variable mismatch for strain {strain}"
                )

            repaired = copy.deepcopy(original)
            met_counts = repair_metabolite_fields(
                repaired, donor, kegg_met_miriams
            )
            n_rxns = len(values(repaired["rxns"]))
            reaction_lengths_before = {
                field: len(values(repaired[field])) for field in REACTION_FIELDS
            }
            reaction_rebuild = any(
                length != n_rxns for length in reaction_lengths_before.values()
            )
            rxn_counts = {"exact": 0, "kegg_database": 0, "custom_empty": 0}
            if reaction_rebuild:
                rxn_counts = repair_reaction_fields(
                    repaired, donor, kegg_rxn_miriams
                )
                family_reaction_rebuilds += 1

            matrix_before = original["rxnGeneMat"]
            if not issparse(matrix_before):
                matrix_before = csc_matrix(matrix_before)
            matrix_before = matrix_before.tocsr()
            rules_before = values(original["rules"])
            mismatch_before = 0
            for row_index, rule in enumerate(rules_before):
                rule_indices = {
                    int(index) - 1
                    for index in re.findall(r"x\((\d+)\)", text(rule))
                }
                matrix_indices = set(
                    matrix_before.indices[
                        matrix_before.indptr[row_index] : matrix_before.indptr[row_index + 1]
                    ]
                )
                mismatch_before += rule_indices != matrix_indices

            genes_before = [text(item) for item in values(original["genes"])]
            orphan_before = set(orphan_gene_names(original))
            removed_genes, used_gene_count, trimmed_gene_metadata_entries = repair_genes(
                repaired
            )
            removed_set = set(removed_genes)
            reconnected_genes = sorted(orphan_before - removed_set)
            family_reconnected.update((strain, gene) for gene in reconnected_genes)
            family_removed += len(removed_genes)
            validation = validate_repair(original, repaired, variable_name)

            output_path = output_dir / source_path.name
            if output_path.exists():
                raise FileExistsError(output_path)
            savemat(
                output_path,
                {variable_name: repaired},
                do_compression=True,
                long_field_names=True,
                oned_as="column",
            )
            reloaded_variable, reloaded = load_model(output_path)
            if reloaded_variable != variable_name:
                raise ValueError(
                    f"{output_path}: model variable changed during serialization"
                )
            validate_repair(original, reloaded, variable_name)

            semantic_match = "not_applicable"
            if family_id == "rich_tfa":
                reference_path = existing_rich_tfa / source_path.name
                if not reference_path.is_file():
                    raise FileNotFoundError(
                        f"missing prior validated rich/TFA repair: {reference_path}"
                    )
                reference_variable, reference_model = load_model(reference_path)
                if reference_variable != variable_name or not equal_value(
                    reference_model, reloaded
                ):
                    raise ValueError(
                        f"{source_path.name}: new repair is not semantically identical "
                        "to the prior validated rich/TFA repair"
                    )
                semantic_match = "true"
                rich_tfa_semantic_matches += 1

            for original_index, gene in enumerate(genes_before, start=1):
                if gene in removed_set:
                    removed_gene_rows.append(
                        {
                            "family_id": family_id,
                            "strain_code": strain,
                            "gene_id": gene,
                            "original_index_1_based": original_index,
                            "reason": "not referenced by any grRules expression",
                        }
                    )

            identity = crosswalk.get(strain, {})
            carbon_source = family["carbon_source"]
            if family["medium"] == "M9" and strain in {"cxy", "lur"}:
                carbon_source = "xylose"
            source_archive_path = Path("Code") / family["source_relative"] / source_path.name
            donor_archive_path = Path("Code") / family["donor_relative"] / donor_path.name
            output_repo_path = Path("models") / "final" / family["output_relative"] / source_path.name
            manifest_rows.append(
                {
                    "family_id": family_id,
                    "strain_code": strain,
                    "gtdb_assignment": identity.get("GTDB_assignment", ""),
                    "assembly_identifier": identity.get("assembly_identifier", ""),
                    "assembly_accession": identity.get("assembly_accession", ""),
                    "medium": family["medium"],
                    "carbon_source": carbon_source,
                    "formulation": family["formulation"],
                    "model_type": family["model_type"],
                    "model_variable": variable_name,
                    "source_file": source_path.name,
                    "source_archive_path": source_archive_path.as_posix(),
                    "source_sha256": sha256_file(source_path),
                    "donor_file": donor_path.name,
                    "donor_archive_path": donor_archive_path.as_posix(),
                    "donor_sha256": sha256_file(donor_path),
                    "output_file": output_path.name,
                    "output_repo_path": output_repo_path.as_posix(),
                    "output_sha256": sha256_file(output_path),
                    "rxns": validation["rxns"],
                    "mets": validation["mets"],
                    "genes_before": len(genes_before),
                    "genes_after": validation["genes"],
                    "unused_genes_removed": len(removed_genes),
                    "gpr_genes_retained": used_gene_count,
                    "gpr_genes_reconnected": len(reconnected_genes),
                    "rule_matrix_mismatches_repaired": mismatch_before,
                    "trimmed_trailing_gene_metadata_entries": trimmed_gene_metadata_entries,
                    "met_exact_donors": met_counts["exact"],
                    "met_kegg_donors": met_counts["kegg_donor"],
                    "met_canonical_donors": met_counts["canonical"],
                    "met_kegg_database_records": met_counts["kegg_database"],
                    "met_minimal_records": met_counts["minimal"],
                    "reaction_fields_rebuilt": str(reaction_rebuild).lower(),
                    "rxn_exact_donors": rxn_counts["exact"],
                    "rxn_kegg_database_records": rxn_counts["kegg_database"],
                    "rxn_custom_empty_records": rxn_counts["custom_empty"],
                    "prior_rich_tfa_semantic_match": semantic_match,
                    "validation_status": "PASS",
                }
            )
            model_validation_rows.append(
                {
                    "family_id": family_id,
                    "strain_code": strain,
                    "output_repo_path": output_repo_path.as_posix(),
                    "rxns": validation["rxns"],
                    "mets": validation["mets"],
                    "genes": validation["genes"],
                    "metabolite_annotation_lengths": {
                        field: len(values(reloaded[field]))
                        for field in METABOLITE_FIELDS
                    },
                    "reaction_annotation_lengths": {
                        field: len(values(reloaded[field]))
                        for field in REACTION_FIELDS
                    },
                    "gene_annotation_lengths": {
                        field: len(values(reloaded[field])) for field in GENE_FIELDS
                    },
                    "rules_length": len(values(reloaded["rules"])),
                    "gr_rules_length": len(values(reloaded["grRules"])),
                    "rxn_gene_matrix_shape": list(reloaded["rxnGeneMat"].shape),
                    "metabolite_annotation_lengths_match_mets": True,
                    "reaction_annotation_lengths_match_rxns": True,
                    "gene_annotation_lengths_match_genes": True,
                    "rxn_gene_matrix_shape_valid": True,
                    "rxn_gene_matrix_orphan_columns": 0,
                    "rules_match_rxn_gene_matrix": True,
                    "mat_reload_succeeded": True,
                    "non_target_fields_preserved": True,
                    "prior_rich_tfa_semantic_match": (
                        True if semantic_match == "true" else None
                    ),
                    "status": "PASS",
                }
            )
            print(
                f"{family_id}/{strain}: PASS; {len(removed_genes)} unused genes "
                f"removed; {len(reconnected_genes)} GPR genes reconnected; "
                f"reaction annotations rebuilt={reaction_rebuild}",
                flush=True,
            )

        if len(family_reconnected) != expected_reconnected:
            raise ValueError(
                f"{family_id}: expected {expected_reconnected} distinct GPR genes to be "
                f"reconnected, found {len(family_reconnected)}"
            )
        family_summaries.append(
            {
                "family_id": family_id,
                "count": len(source_files),
                "removed": family_removed,
                "reconnected": len(family_reconnected),
                "reaction_rebuilds": family_reaction_rebuilds,
            }
        )

    if len(manifest_rows) != 80:
        raise ValueError(f"expected 80 repaired outputs, found {len(manifest_rows)}")
    if rich_tfa_semantic_matches != 21:
        raise ValueError(
            f"expected 21 prior rich/TFA semantic matches, found {rich_tfa_semantic_matches}"
        )

    write_tsv(model_manifest_path, list(manifest_rows[0]), manifest_rows)
    write_tsv(
        removed_genes_path,
        [
            "family_id",
            "strain_code",
            "gene_id",
            "original_index_1_based",
            "reason",
        ],
        removed_gene_rows,
    )
    validation_report_path.parent.mkdir(parents=True, exist_ok=True)
    report_lines = [
        "# Model package validation",
        "",
        "All 80 publication models passed structural and serialization validation.",
        "The untouched source archive was read only; corrected copies are stored under",
        "`models/final/`.",
        "",
        "| Family | Models | Unused genes removed | GPR genes reconnected | Reaction annotation rebuilds |",
        "|---|---:|---:|---:|---:|",
    ]
    for summary in family_summaries:
        report_lines.append(
            f"| {summary['family_id']} | {summary['count']} | {summary['removed']} | "
            f"{summary['reconnected']} | {summary['reaction_rebuilds']} |"
        )
    report_lines.extend(
        [
            "",
            "Validation checks:",
            "",
            "- All ten repaired metabolite annotation arrays match `mets`.",
            "- KEGG, Rhea, MIRIAM, and SBO reaction annotation arrays match `rxns`.",
            "- `rxnGeneMat` has shape reactions-by-genes and has no orphan columns.",
            "- Every numeric `rules` index exactly matches its `rxnGeneMat` row.",
            "- Every field outside the explicitly repaired metadata/GPR fields is unchanged.",
            "- Every output reloads successfully from its MAT file.",
            "- The 21 rich/TFA models are semantically identical to the earlier validated Python-fixed set.",
            "- Five distinct rich-medium and four distinct M9 GPR genes were reconnected in both FBA and TFA formulations.",
            "",
            f"KEGG metabolite snapshot SHA-256: `{sha256_file(kegg_mets_path)}`",
            "",
            f"KEGG reaction snapshot SHA-256: `{sha256_file(kegg_rxns_path)}`",
            "",
        ]
    )
    with validation_report_path.open("x", encoding="utf-8") as handle:
        handle.write("\n".join(report_lines))

    validation_payload = {
        "status": "PASS",
        "expected_model_count": 80,
        "actual_model_count": len(model_validation_rows),
        "rich_tfa_semantic_identity": {
            "expected": 21,
            "actual": rich_tfa_semantic_matches,
            "status": "PASS",
        },
        "invariants": {
            "metabolite_annotation_lengths_match_mets": True,
            "reaction_annotation_lengths_match_rxns": True,
            "gene_annotation_lengths_match_genes": True,
            "rxn_gene_matrix_shape_valid": True,
            "rxn_gene_matrix_orphan_columns": 0,
            "rules_match_rxn_gene_matrix": True,
            "mat_reload_succeeded": True,
            "non_target_fields_preserved": True,
        },
        "families": family_summaries,
        "models": model_validation_rows,
    }
    with validation_json_path.open("x", encoding="utf-8") as handle:
        json.dump(validation_payload, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"Wrote and validated 80 repaired models under {output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
