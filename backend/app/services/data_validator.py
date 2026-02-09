from __future__ import annotations

from pathlib import Path


def validate_h5ad(
    data_path: Path,
    use_drug_structure: bool,
) -> dict[str, object]:
    errors: list[str] = []
    warnings: list[str] = []

    if not data_path.exists():
        return {
            "valid": False,
            "errors": [f"Input file does not exist: {data_path}"],
            "warnings": warnings,
        }

    try:
        import scanpy as sc
    except ModuleNotFoundError as exc:  # noqa: PERF203
        return {
            "valid": False,
            "errors": [
                "scanpy is not installed in current environment. "
                "Install dependencies before dataset validation."
            ],
            "warnings": warnings,
            "detail": str(exc),
        }

    adata = sc.read_h5ad(str(data_path))
    obs_columns = set(adata.obs.columns.tolist())

    if "Group" not in obs_columns:
        errors.append("Missing required obs column: Group")

    if use_drug_structure:
        if "SMILES" not in obs_columns:
            errors.append("Missing required obs column for drug mode: SMILES")
        if "dose" not in obs_columns:
            errors.append("Missing required obs column for drug mode: dose")

    if adata.n_obs <= 0:
        errors.append("No cells in dataset (n_obs <= 0).")
    if adata.n_vars <= 0:
        errors.append("No genes in dataset (n_vars <= 0).")

    if adata.n_vars < 50:
        warnings.append(
            "Low gene count detected (n_vars < 50). This may reduce model quality."
        )

    return {
        "valid": len(errors) == 0,
        "errors": errors,
        "warnings": warnings,
        "summary": {
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "obs_columns": sorted(obs_columns),
        },
    }
