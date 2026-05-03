"""Tests for atomview structure utility helpers."""

from pathlib import Path

import biotite.structure as struc
import numpy as np
import pytest
from biotite.structure import AtomArray

from atomview import load_structure
from atomview.utils import reassign_chain_ids, to_cif_string

DEFAULT_ATOM_COUNT = 2


def make_atom_array(*, include_nan: bool = False, include_hetero: bool = False) -> AtomArray:
    """Create a minimal AtomArray for utility tests."""
    atom_count = 3 if include_hetero else 2
    structure = struc.AtomArray(atom_count)
    structure.coord = np.array(
        [
            [0.0, 0.0, 0.0],
            [np.nan if include_nan else 1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ][:atom_count]
    )
    structure.chain_id = np.array(["X", "X", "X"][:atom_count])
    structure.res_id = np.array([1, 1, 2][:atom_count])
    structure.res_name = np.array(["GLY", "GLY", "LIG"][:atom_count])
    structure.atom_name = np.array(["N", "CA", "C1"][:atom_count])
    structure.element = np.array(["N", "C", "C"][:atom_count])
    structure.hetero = np.array([False, False, include_hetero][:atom_count])
    return structure


def test_to_cif_string_returns_mmcif_content() -> None:
    """CIF conversion should produce a non-empty mmCIF string."""
    cif = to_cif_string(make_atom_array(), id="example")

    assert cif.startswith("data_example")
    assert "_atom_site." in cif
    assert "GLY" in cif


def test_to_cif_string_can_skip_nan_coordinates() -> None:
    """Atoms with NaN coordinates should be removed when requested."""
    cif = to_cif_string(make_atom_array(include_nan=True), include_nan_coords=False)

    assert "nan" not in cif.lower()


def test_load_structure_reads_cif_path(tmp_path: Path) -> None:
    """Structure loading should accept a path to a CIF file."""
    cif_path = tmp_path / "structure.cif"
    cif_path.write_text(to_cif_string(make_atom_array()), encoding="utf-8")

    loaded = load_structure(cif_path)

    assert isinstance(loaded, AtomArray)
    assert len(loaded) == DEFAULT_ATOM_COUNT


def test_load_structure_reads_pdb_path() -> None:
    """Structure loading should accept a path to a PDB file."""
    pdb_path = Path(__file__).resolve().parent / "data" / "4tos.pdb"
    loaded = load_structure(pdb_path)

    assert isinstance(loaded, AtomArray)
    assert len(loaded) > 0


def test_load_structure_rejects_missing_path(tmp_path: Path) -> None:
    """Missing Path inputs should raise FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        load_structure(tmp_path / "missing.cif")


def test_load_structure_rejects_unsupported_extension(tmp_path: Path) -> None:
    """Existing files with unsupported extensions should raise ValueError."""
    bad_path = tmp_path / "structure.txt"
    bad_path.write_text("not a cif", encoding="utf-8")

    with pytest.raises(ValueError, match="Unsupported file format"):
        load_structure(bad_path)


def test_reassign_chain_ids_returns_copy_and_splits_hetero_atoms() -> None:
    """Hetero and non-hetero atoms should not share reassigned chain IDs."""
    structure = make_atom_array(include_hetero=True)

    reassigned = reassign_chain_ids(structure)

    assert reassigned is not structure
    assert structure.chain_id.tolist() == ["X", "X", "X"]
    assert reassigned.chain_id.tolist() == ["A", "A", "B"]
