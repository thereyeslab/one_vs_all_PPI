from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml


VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class Protein:
    protein_id: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)


def extract_fasta_id(header: str) -> str:
    """
    Extract clean IDs from UniProt or simple FASTA headers.

    Examples:
    >sp|P27636|RFA1_YEAST ... GN=RFA1 ...
    -> RFA1

    >RFA1
    -> RFA1
    """
    header = header.strip().lstrip(">")

    gene_match = re.search(r"\bGN=([A-Za-z0-9_.-]+)", header)
    if gene_match:
        return gene_match.group(1)

    first_token = header.split()[0]

    if "|" in first_token:
        parts = first_token.split("|")
        if len(parts) >= 3:
            entry_name = parts[2]
            return entry_name.split("_")[0] if "_" in entry_name else entry_name

    return first_token


def validate_sequence(protein_id: str, sequence: str) -> str:
    """
    Validate sequence for AlphaFold Server.
    Only 20 standard amino acids allowed.
    """
    sequence = sequence.upper().replace(" ", "").replace("*", "")

    invalid = sorted(set(sequence) - VALID_AA)

    if invalid:
        raise ValueError(
            f"Protein '{protein_id}' contains invalid amino acid characters: {invalid}"
        )

    if not sequence:
        raise ValueError(f"Protein '{protein_id}' has empty sequence.")

    return sequence


def read_fasta(fasta_path: Path) -> dict[str, Protein]:
    proteins: dict[str, Protein] = {}
    current_id: str | None = None
    current_seq: list[str] = []

    with fasta_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    proteins[current_id] = Protein(current_id, "".join(current_seq))

                current_id = extract_fasta_id(line)
                current_seq = []
            else:
                current_seq.append(line)

    if current_id is not None:
        seq = validate_sequence(current_id, "".join(current_seq))
        proteins[current_id] = Protein(current_id, seq)

    if not proteins:
        raise ValueError(f"No proteins found in FASTA: {fasta_path}")

    return proteins


def make_chain_ids(n: int) -> list[str]:
    """
    Generate AF3 entity IDs.

    For now supports up to 26 chains: A-Z.
    That is enough for pooled AF3 jobs with ~5000 aa total.
    """
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    if n > len(alphabet):
        raise ValueError(
            f"This script currently supports up to {len(alphabet)} chains per pool. "
            f"Requested {n}."
        )

    return list(alphabet[:n])


def parse_protein_ids_from_pool_row(row: pd.Series) -> list[str]:
    """
    Our pool generator writes protein_ids as semicolon-separated IDs.
    Example:
    RFA1;RAD52;RAD24
    """
    if "protein_ids" not in row:
        raise ValueError("pools.tsv must contain a 'protein_ids' column.")

    raw = str(row["protein_ids"]).strip()

    if not raw:
        raise ValueError(f"Empty protein_ids for pool row: {row.to_dict()}")

    return [pid.strip() for pid in raw.split(";") if pid.strip()]


def build_af3_json(
    pool_name: str,
    protein_ids: list[str],
    proteins: dict[str, Protein],
    model_seeds: list[int],
    dialect: str,
    version: int,
) -> tuple[dict, list[dict]]:
    chain_ids = make_chain_ids(len(protein_ids))

    sequences = []
    mapping_rows = []

    for chain_id, protein_id in zip(chain_ids, protein_ids):
        if protein_id not in proteins:
            raise ValueError(
                f"Protein ID '{protein_id}' from pools.tsv not found in FASTA."
            )

        protein = proteins[protein_id]

        sequences.append(
            {
                "proteinChain": {
                    "id": chain_id,
                    "sequence": protein.sequence,
                }
            }
        )

        mapping_rows.append(
            {
                "pool_name": pool_name,
                "chain_id": chain_id,
                "protein_id": protein_id,
                "protein_length": protein.length,
            }
        )

    af3_input = {
        "name": pool_name,
        "modelSeeds": model_seeds,
        "sequences": sequences,
        "dialect": dialect,
        "version": version,
    }

    return af3_input, mapping_rows


def write_json(path: Path, payload: dict, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"File already exists and overwrite=false: {path}")

    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def load_config(config_path: Path) -> dict:
    with config_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def json_builder(config: dict) -> None:


    fasta_path = Path(config["input"]["fasta_path"])
    pools_tsv = Path(config["input"]["pools_tsv"])
    json_dir = Path(config["output"]["json_dir"])

    model_seeds = config["alphafold3"].get("model_seeds", [1])
    server_cfg = config.get("alphafold3", {})
    options_cfg = config.get("json_preparation_options", {})
    dialect = config["alphafold3"].get("dialect", "alphafold3")
    version = int(config["alphafold3"].get("version", 1))
    overwrite = bool(config["options"].get("overwrite", True))

    json_dir.mkdir(parents=True, exist_ok=True)

    proteins = read_fasta(fasta_path)
    pools_df = pd.read_csv(pools_tsv, sep="\t")

    all_mapping_rows = []

    for _, row in pools_df.iterrows():
        pool_id = str(row["pool_id"])
        protein_ids = parse_protein_ids_from_pool_row(row)

        af3_input, mapping_rows = build_af3_json(
            pool_name=pool_id,
            protein_ids=protein_ids,
            proteins=proteins,
            model_seeds=model_seeds,
            dialect=dialect,
            version=version,
        )

        output_json = json_dir / f"{pool_id}.json"
        write_json(output_json, af3_input, overwrite=overwrite)

        all_mapping_rows.extend(mapping_rows)

    mapping_df = pd.DataFrame(all_mapping_rows)
    mapping_path = json_dir.parent / "pool_chain_mapping.tsv"
    mapping_df.to_csv(mapping_path, sep="\t", index=False)

    print("Done.")
    print(f"Proteins loaded: {len(proteins)}")
    print(f"Pools loaded: {len(pools_df)}")
    print(f"JSON files written: {len(pools_df)}")
    print(f"JSON directory: {json_dir}")
    print(f"Chain mapping written: {mapping_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()

    config = load_config(args.config)
    json_builder(config)