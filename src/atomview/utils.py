"""Utility helpers for working with mmCIF/BinaryCIF structures."""

from __future__ import annotations

import io
from pathlib import Path
from string import ascii_uppercase
from typing import Literal

import biotite.structure as struc
import numpy as np
from biotite.structure import AtomArray, AtomArrayStack
from biotite.structure.io import pdbx


def to_cif_string(
    structure: AtomArray,
    *,
    id: str = "unknown_id",
    include_nan_coords: bool = True,
    include_bonds: bool = True,
    extra_fields: list[str] | Literal["all"] = "all",
    _allow_ambiguous_bond_annotations: bool = False,
) -> str:
    """Convert an ``AtomArray`` into mmCIF text."""
    structure = structure.copy()

    if not include_nan_coords:
        structure = structure[~np.isnan(structure.coord).any(axis=1)]

    if include_bonds and structure.bonds is not None:
        mask = structure.bonds._bonds[:, 2] == struc.bonds.BondType.COORDINATION
        structure.bonds._bonds[mask, 2] = struc.bonds.BondType.SINGLE

    cif_file = pdbx.CIFFile()

    if extra_fields == "all":
        standard_fields = {
            "chain_id",
            "res_id",
            "res_name",
            "atom_name",
            "atom_id",
            "element",
            "ins_code",
            "hetero",
            "altloc_id",
            "charge",
            "occupancy",
            "b_factor",
        }
        extra_fields = [
            field
            for field in structure.get_annotation_categories()
            if field not in standard_fields
        ]

    pdbx.set_structure(
        cif_file,
        structure,
        data_block=id,
        include_bonds=include_bonds,
        # extra_fields=extra_fields,
    )

    buffer = io.StringIO()
    cif_file.write(buffer)
    buffer.seek(0)
    return buffer.getvalue()


def load_structure(
    source: str | Path | io.StringIO | io.BytesIO,
    *,
    include_bonds: bool = True,
    model: int | None = None,
) -> AtomArray:
    """Load an ``AtomArray`` from a path, buffer, or CIF string."""
    if isinstance(source, Path):
        if not source.exists():
            raise FileNotFoundError(f"File not found: {source}")
        return _read_structure_file(source, include_bonds=include_bonds, model=model)

    if isinstance(source, str):
        try:
            path = Path(source)
            if path.exists():
                return _read_structure_file(path, include_bonds=include_bonds, model=model)
        except OSError:
            pass
        return _structure_from_cif_file(
            pdbx.CIFFile.read(io.StringIO(source)),
            include_bonds=include_bonds,
            model=model,
        )

    if isinstance(source, io.StringIO):
        source.seek(0)
        return _structure_from_cif_file(
            pdbx.CIFFile.read(source),
            include_bonds=include_bonds,
            model=model,
        )

    if isinstance(source, io.BytesIO):
        source.seek(0)
        try:
            cif_file = pdbx.BinaryCIFFile.read(source)
        except (TypeError, ValueError, RuntimeError):
            source.seek(0)
            cif_file = pdbx.CIFFile.read(source)
        return _structure_from_cif_file(cif_file, include_bonds=include_bonds, model=model)

    raise ValueError(f"Unsupported source type: {type(source)!r}")


def _structure_from_cif_file(
    cif_file: pdbx.CIFFile | pdbx.BinaryCIFFile,
    *,
    include_bonds: bool,
    model: int | None,
) -> AtomArray:
    """Return a single model ``AtomArray`` from a loaded CIF object."""
    structure = pdbx.get_structure(
        cif_file,
        model=model if model is not None else 1,
        include_bonds=include_bonds,
    )
    if isinstance(structure, AtomArrayStack):
        return structure[0] if model is None else structure[model - 1]
    return structure


def _read_structure_file(
    path: Path,
    *,
    include_bonds: bool,
    model: int | None,
) -> AtomArray:
    """Load a structure from a supported on-disk CIF format."""
    suffixes = [suffix.lower() for suffix in path.suffixes]
    if suffixes[-1:] == [".bcif"] or suffixes[-2:] == [".bcif", ".gz"]:
        cif_file = pdbx.BinaryCIFFile.read(str(path))
    elif suffixes[-1:] == [".cif"] or suffixes[-2:] == [".cif", ".gz"]:
        cif_file = pdbx.CIFFile.read(str(path))
    else:
        raise ValueError(f"Unsupported file format for {path.name}. Use .cif or .bcif")
    return _structure_from_cif_file(cif_file, include_bonds=include_bonds, model=model)


_load_structure = load_structure


def reassign_chain_ids(structure: AtomArray) -> AtomArray:
    """Reassign chain IDs so polymer and hetero groups do not collide."""
    structure = structure.copy()

    is_hetero = (
        structure.hetero
        if "hetero" in structure.get_annotation_categories()
        else np.zeros(len(structure), dtype=bool)
    )
    hetero_mask = is_hetero.astype(bool)
    non_hetero_mask = ~hetero_mask
    new_chain_id = structure.chain_id.copy()

    letters = list(ascii_uppercase)
    non_hetero_chain_count = 0

    if np.any(non_hetero_mask):
        non_hetero_structure = structure[non_hetero_mask]
        chain_starts = list(struc.get_chain_starts(non_hetero_structure))
        chain_starts.append(len(non_hetero_structure))
        non_hetero_indices = np.where(non_hetero_mask)[0]

        for i in range(len(chain_starts) - 1):
            start_idx = chain_starts[i]
            end_idx = chain_starts[i + 1]
            structure_indices = non_hetero_indices[start_idx:end_idx]
            new_chain_id[structure_indices] = letters[i]
            non_hetero_chain_count += 1

    if np.any(hetero_mask):
        hetero_structure = structure[hetero_mask]
        chain_starts = list(struc.get_chain_starts(hetero_structure))
        chain_starts.append(len(hetero_structure))
        hetero_indices = np.where(hetero_mask)[0]

        for i in range(len(chain_starts) - 1):
            start_idx = chain_starts[i]
            end_idx = chain_starts[i + 1]
            structure_indices = hetero_indices[start_idx:end_idx]
            new_chain_id[structure_indices] = letters[non_hetero_chain_count + i]

    structure.chain_id = new_chain_id
    return structure
