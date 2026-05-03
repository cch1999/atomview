"""Utility functions for working with CIF files using pure biotite.

Adapted from atomworks.io.utils.visualize:
https://github.com/RosettaCommons/atomworks/blob/production/src/atomworks/io/utils/visualize.py
"""

import io
from pathlib import Path
from string import ascii_uppercase
from typing import Literal

import biotite.structure as struc
import numpy as np
from biotite.structure import AtomArray, AtomArrayStack
from biotite.structure.io import pdbx
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdb import get_structure as pdb_get_structure


def to_cif_string(
    structure: AtomArray,
    *,
    id: str = "unknown_id",
    include_nan_coords: bool = True,
    include_bonds: bool = True,
    extra_fields: list[str] | Literal["all"] = "all",
    _allow_ambiguous_bond_annotations: bool = False,
) -> str:
    """Convert an AtomArray structure to a CIF formatted string for visualization.

    This function uses pure biotite to convert an AtomArray to CIF format,
    optimized for use with molecular viewers like py3Dmol.

    Args:
        structure (AtomArray): The atomic structure to be converted.
        id (str): The ID of the entry. This will be used as the data block name.
        include_nan_coords (bool): Whether to write NaN coordinates in the CIF file.
        include_bonds (bool): Whether to write bonds in the CIF file.
        extra_fields (list[str] | Literal["all"]): Additional atom_array annotations to include
            in the CIF file.
        _allow_ambiguous_bond_annotations (bool): Private argument, not meant for public use.
            If True, allows ambiguous bond annotations.

    Returns:
        str: The CIF formatted string representation of the structure.
    """
    # Create a copy to avoid modifying the original structure
    structure = structure.copy()

    # Handle NaN coordinates
    if not include_nan_coords:
        # Filter out atoms with NaN coordinates
        valid_coords = ~np.isnan(structure.coord).any(axis=1)
        structure = structure[valid_coords]

    # Handle bonds - convert coordination bonds to single bonds if needed
    if include_bonds and structure.bonds is not None:
        mask = structure.bonds._bonds[:, 2] == struc.bonds.BondType.COORDINATION
        structure.bonds._bonds[mask, 2] = struc.bonds.BondType.SINGLE

    # Create CIF file
    cif_file = pdbx.CIFFile()

    # Handle extra_fields="all" case
    if extra_fields == "all":
        # Get all annotation categories except standard ones
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
            field for field in structure.get_annotation_categories() if field not in standard_fields
        ]
    # Set the structure in the CIF file (this creates the block)
    pdbx.set_structure(
        cif_file,
        structure,
        data_block=id,
        include_bonds=include_bonds,
        # extra_fields=extra_fields
    )

    # Write to buffer and return as string
    buffer = io.StringIO()
    cif_file.write(buffer)
    buffer.seek(0)
    return buffer.getvalue()


def load_structure(  # noqa: PLR0912, PLR0915
    source: str | Path | io.StringIO | io.BytesIO,
    *,
    include_bonds: bool = True,
    model: int | None = None,
) -> AtomArray:
    """Load an AtomArray structure from mmCIF, BinaryCIF, or PDB files (or CIF text).

    This function uses biotite I/O, optimised for use with molecular viewers.

    Args:
        source: The source to load from. Can be:
            - A file path (str or Path) to a .cif, .bcif, or .pdb file
            - A StringIO/BytesIO buffer containing CIF data
            - A CIF string
        include_bonds (bool): Whether to include bonds in the structure. Defaults to True.
        model (int | None): The model number to load. If None and multiple models exist,
            returns the first model. Defaults to None.

    Returns:
        AtomArray: The loaded atomic structure.

    Raises:
        FileNotFoundError: If the file path does not exist.
        ValueError: If the source format is not supported.
    """
    # Determine the source type and load CIF file
    if isinstance(source, Path):
        # Path object - treat as file path
        path = source
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")
    elif isinstance(source, str):
        # String - check if it's a file path or CIF string
        path = Path(source)
        if path.exists():
            # It's a file path
            source = path
        else:
            # Assume it's a CIF string
            cif_file = pdbx.CIFFile.read(io.StringIO(source))
            structure = pdbx.get_structure(
                cif_file,
                model=model if model is not None else 1,
                include_bonds=include_bonds,
            )
            if isinstance(structure, AtomArrayStack):
                structure = structure[0] if model is None else structure[model - 1]
            return structure
    elif isinstance(source, io.StringIO):
        # StringIO buffer - assume regular CIF
        source.seek(0)
        cif_file = pdbx.CIFFile.read(source)
        structure = pdbx.get_structure(
            cif_file,
            model=model if model is not None else 1,
            include_bonds=include_bonds,
        )
        if isinstance(structure, AtomArrayStack):
            structure = structure[0] if model is None else structure[model - 1]
        return structure
    elif isinstance(source, io.BytesIO):
        # BytesIO buffer - try BinaryCIF first, fall back to regular CIF
        source.seek(0)
        try:
            cif_file = pdbx.BinaryCIFFile.read(source)
        except Exception:  # noqa: BLE001
            source.seek(0)
            cif_file = pdbx.CIFFile.read(source)
        structure = pdbx.get_structure(
            cif_file,
            model=model if model is not None else 1,
            include_bonds=include_bonds,
        )
        if isinstance(structure, AtomArrayStack):
            structure = structure[0] if model is None else structure[model - 1]
        return structure
    # Handle file paths (Path or string converted to Path)
    # At this point, source is either a Path or was converted to Path above
    if isinstance(source, Path):
        path = source
        # Determine file type from extension
        suffix = path.suffix.lower()
        if suffix == ".bcif" or (
            suffix == ".gz" and len(path.suffixes) > 1 and path.suffixes[-2].lower() == ".bcif"
        ):
            # Binary CIF file
            cif_file = pdbx.BinaryCIFFile.read(str(path))
        elif suffix == ".cif" or (
            suffix == ".gz" and len(path.suffixes) > 1 and path.suffixes[-2].lower() == ".cif"
        ):
            # Regular CIF file
            cif_file = pdbx.CIFFile.read(str(path))
        elif suffix == ".pdb" or (
            suffix == ".gz" and len(path.suffixes) > 1 and path.suffixes[-2].lower() == ".pdb"
        ):
            pdb_file = PDBFile.read(str(path))
            structure = pdb_get_structure(
                pdb_file,
                model=model if model is not None else 1,
                include_bonds=include_bonds,
            )
            if isinstance(structure, AtomArrayStack):
                if model is not None:
                    structure = structure[model - 1]
                else:
                    structure = structure[0]
            return structure
        else:
            raise ValueError(f"Unsupported file format: {suffix}. Use .cif, .bcif, or .pdb")

        # Load structure from CIF file
        structure = pdbx.get_structure(
            cif_file,
            model=model if model is not None else 1,
            include_bonds=include_bonds,
        )

        # Ensure we return an AtomArray, not AtomArrayStack
        if isinstance(structure, AtomArrayStack):
            if model is not None:
                structure = structure[model - 1]  # Convert to 0-based index
            else:
                structure = structure[0]  # Use first model

        return structure

    raise ValueError(f"Unexpected source type after processing: {type(source)}")


_load_structure = load_structure


def reassign_chain_ids(structure: AtomArray) -> AtomArray:
    """Reassign chain IDs to separate hetero and non-hetero atoms.

    This function ensures that:
    - Non-hetero atoms (polymers) get unique chain IDs starting from A, B, C...
    - Hetero atoms (ligands) get unique chain IDs starting after non-hetero chains
    - Hetero and non-hetero chains never share the same chain ID

    This prevents rendering style conflicts when polymers and ligands share chain IDs.

    Args:
        structure (AtomArray): The atomic structure to reassign chain IDs for.

    Returns:
        AtomArray: A copy of the structure with reassigned chain IDs.
    """
    # Create a copy to avoid modifying the original structure
    structure = structure.copy()

    # Separate hetero and non-hetero atoms
    is_hetero = (
        structure.hetero
        if "hetero" in structure.get_annotation_categories()
        else np.zeros(len(structure), dtype=bool)
    )
    hetero_mask = is_hetero.astype(bool)
    non_hetero_mask = ~hetero_mask

    # Create a copy of chain_id array to modify
    new_chain_id = structure.chain_id.copy()

    # Handle non-hetero atoms: use biotite's chain detection to assign unique chain IDs
    letters = list(ascii_uppercase)
    non_hetero_chain_count = 0

    if np.any(non_hetero_mask):
        non_hetero_structure = structure[non_hetero_mask]
        # Use biotite's chain detection on non-hetero atoms only
        # This properly separates chains even if original chain_id was wrong
        chain_starts = list(struc.get_chain_starts(non_hetero_structure))
        chain_starts.append(len(non_hetero_structure))

        # Assign sequential chain IDs (A, B, C, ...) to detected non-hetero chains
        non_hetero_indices = np.where(non_hetero_mask)[0]

        for i in range(len(chain_starts) - 1):
            start_idx = chain_starts[i]
            end_idx = chain_starts[i + 1]
            # Map indices from non-hetero subset back to full structure
            structure_indices = non_hetero_indices[start_idx:end_idx]
            new_chain_id[structure_indices] = letters[i]
            non_hetero_chain_count += 1

    # Handle hetero atoms: use biotite's chain detection to assign unique chain IDs
    # Start from letters after non-hetero chains to avoid conflicts
    if np.any(hetero_mask):
        hetero_structure = structure[hetero_mask]
        # Use biotite's chain detection on hetero atoms only
        # This properly separates hetero chains even if original chain_id was wrong
        chain_starts = list(struc.get_chain_starts(hetero_structure))
        chain_starts.append(len(hetero_structure))

        # Assign sequential chain IDs starting after non-hetero chains
        # This ensures hetero and non-hetero chains don't overlap
        hetero_indices = np.where(hetero_mask)[0]

        for i in range(len(chain_starts) - 1):
            start_idx = chain_starts[i]
            end_idx = chain_starts[i + 1]
            # Map indices from hetero subset back to full structure
            structure_indices = hetero_indices[start_idx:end_idx]
            # Use letters starting after non-hetero chains
            new_chain_id[structure_indices] = letters[non_hetero_chain_count + i]

    # Update structure with new chain IDs
    structure.chain_id = new_chain_id

    return structure
