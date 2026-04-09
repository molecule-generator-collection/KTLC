#!/usr/bin/env python3
"""Predict the KTLC optimal range-separated parameter μ from SMILES.

This module mirrors the feature-generation logic used in the Colab notebook so
that predictions remain consistent after moving the workflow into a standalone
Python script inside the KTLC repository.

Expected repository layout:
KTLC/
├── KTLC_pred_mu.py
└── model/
    ├── mordred_index_mask.pkl
    └── model_lightGBM_All_Desc_0.pkl
    └── ...
"""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from mordred import Calculator, descriptors
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.warning')

DEFAULT_P_MODEL = "lightGBM"
DEFAULT_FEATURE = "All_Desc"
N_FOLDS = 5

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_MODEL_DIR = BASE_DIR / "model"


def _validate_smiles(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return mol


def _bitvect_to_numpy(bitvect) -> np.ndarray:
    return np.fromiter((int(x) for x in bitvect.ToBitString()), dtype=np.int8)


def _mols_from_smiles(smiles_list: Sequence[str]) -> list[Chem.Mol]:
    return [_validate_smiles(smi) for smi in smiles_list]


def smi2descriptor(
    smi_list: Sequence[str],
    mordred_index_list: Sequence[int],
) -> np.ndarray:
    """Generate descriptors exactly as in the Colab notebook."""
    mol_list = _mols_from_smiles(smi_list)

    morgan_fp_list = np.vstack(
        [_bitvect_to_numpy(AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)) for mol in mol_list]
    )

    maccs_fp_list = np.vstack(
        [_bitvect_to_numpy(AllChem.GetMACCSKeysFingerprint(mol)) for mol in mol_list]
    )

    params = FingerprintMols.FingerprinterDetails()
    params.tgtDensity = 0
    topol_fp_list = np.vstack(
        [_bitvect_to_numpy(FingerprintMols.FingerprintMol(mol, **params.__dict__)) for mol in mol_list]
    )

    descriptor_names = [descriptor_name[0] for descriptor_name in Descriptors._descList]
    descriptor_calculation = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)
    rdkit_desc_list = np.asarray(
        [descriptor_calculation.CalcDescriptors(mol_temp) for mol_temp in mol_list],
        dtype=float,
    )

    calc = Calculator(descriptors, ignore_3D=True)
    df_descriptors_mordred = calc.pandas(mol_list)

    # Keep the filtering logic aligned with the notebook:
    # convert to string -> drop columns that contain alphabetic characters
    # anywhere -> cast back to float -> apply the stored mask.
    df_descriptors = df_descriptors_mordred.astype(str)
    masks = df_descriptors.apply(lambda column: column.str.contains(r"[a-zA-Z]", na=False))
    df_descriptors = df_descriptors[~masks]
    df_descriptors = df_descriptors.astype(float)
    df_descriptors = df_descriptors.values[:, list(mordred_index_list)]

    all_desc_list = np.concatenate(
        [morgan_fp_list, topol_fp_list, maccs_fp_list, rdkit_desc_list, df_descriptors],
        axis=1,
    )
    return all_desc_list


def load_selected_index_mask(model_dir: str | Path = DEFAULT_MODEL_DIR):
    model_dir = Path(model_dir)
    mask_path = model_dir / "mordred_index_mask.pkl"
    if not mask_path.exists():
        raise FileNotFoundError(f"Missing mask file: {mask_path}")
    with mask_path.open("rb") as f:
        return pickle.load(f)


def _load_model(model_dir: Path, p_model: str, feature: str, fold: int):
    model_path = model_dir / f"model_{p_model}_{feature}_{fold}.pkl"
    if not model_path.exists():
        raise FileNotFoundError(f"Missing model file: {model_path}")
    with model_path.open("rb") as f:
        return pickle.load(f)


def smiles2mu(
    smi_list: Sequence[str],
    model_dir: str | Path = DEFAULT_MODEL_DIR,
    p_model: str = DEFAULT_P_MODEL,
    feature: str = DEFAULT_FEATURE,
):
    model_dir = Path(model_dir)
    selected_index_mask = load_selected_index_mask(model_dir)
    desc = smi2descriptor(smi_list, selected_index_mask)

    predicted_mu_list = []
    for fold in range(N_FOLDS):
        loaded_model = _load_model(model_dir, p_model, feature, fold)
        predicted_mu = loaded_model.predict(desc)
        predicted_mu_list.append(predicted_mu)

    predicted_mu_array = np.asarray(predicted_mu_list, dtype=float)
    return np.mean(predicted_mu_array, axis=0), np.std(predicted_mu_array, axis=0)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Predict KTLC optimal μ values from SMILES."
    )
    parser.add_argument(
        "smiles",
        nargs="*",
        help="One or more SMILES strings. Use --input-file for batch prediction.",
    )
    parser.add_argument(
        "-i",
        "--input-file",
        type=Path,
        help="Text file with one SMILES per line.",
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=DEFAULT_MODEL_DIR,
        help="Directory containing mordred_index_mask.pkl and model_*.pkl files.",
    )
    parser.add_argument(
        "--p-model",
        default=DEFAULT_P_MODEL,
        help=f"Model name prefix (default: {DEFAULT_P_MODEL}).",
    )
    parser.add_argument(
        "--feature",
        default=DEFAULT_FEATURE,
        help=f"Feature name prefix (default: {DEFAULT_FEATURE}).",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print predictions as JSON.",
    )
    return parser


def _read_smiles_from_args(args: argparse.Namespace) -> list[str]:
    smiles_list: list[str] = []
    if args.input_file is not None:
        if not args.input_file.exists():
            raise FileNotFoundError(f"Input file not found: {args.input_file}")
        with args.input_file.open() as f:
            smiles_list.extend(line.strip() for line in f if line.strip())

    smiles_list.extend(args.smiles)

    if not smiles_list:
        raise ValueError("No SMILES were provided. Pass SMILES directly or via --input-file.")
    return smiles_list


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()
    smiles_list = _read_smiles_from_args(args)

    mu_mean, mu_std = smiles2mu(
        smiles_list,
        model_dir=args.model_dir,
        p_model=args.p_model,
        feature=args.feature,
    )

    rows = [
        {
            "smiles": smi,
            "predicted_mu": float(mean),
            "predicted_mu_std": float(std),
        }
        for smi, mean, std in zip(smiles_list, mu_mean, mu_std)
    ]

    if args.json:
        print(json.dumps(rows, ensure_ascii=False, indent=2))
    else:
        for row in rows:
            print(f'SMILES: {row["smiles"]}')
            print(f'Predicted μ: {row["predicted_mu"]:.5f}')
            print(f'Predicted μ (std): {row["predicted_mu_std"]:.5f}')
            print("")


if __name__ == "__main__":
    main()
