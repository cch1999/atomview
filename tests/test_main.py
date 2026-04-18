"""Smoke tests for the public atomview API."""

from pathlib import Path

import py3Dmol

from atomview import load_structure, to_cif_string, view

FIXTURES_DIR = Path(__file__).resolve().parent.parent
STRUCTURE_PATH = FIXTURES_DIR / "4hhb.cif"
XYZ_DIMENSIONS = 3


def test_load_structure_from_path() -> None:
    """Load a bundled mmCIF structure from disk."""
    structure = load_structure(STRUCTURE_PATH)

    assert structure.array_length() > 0
    assert structure.coord.shape[1] == XYZ_DIMENSIONS


def test_load_structure_from_cif_string() -> None:
    """Load a structure from an in-memory CIF string."""
    cif_text = STRUCTURE_PATH.read_text()
    structure = load_structure(cif_text)

    assert structure.array_length() > 0
    assert structure.res_name[0]


def test_to_cif_string_returns_mmcif_text() -> None:
    """Serialize a structure back to mmCIF text."""
    structure = load_structure(STRUCTURE_PATH)
    cif_text = to_cif_string(structure, id="4hhb")

    assert cif_text.startswith("data_4hhb")
    assert "_atom_site.Cartn_x" in cif_text


def test_view_returns_py3dmol_view() -> None:
    """Build a viewer object from the package root API."""
    structure = load_structure(STRUCTURE_PATH)
    viewer = view(structure, show_surface=False, show_hover=False)

    assert isinstance(viewer, py3Dmol.view)
