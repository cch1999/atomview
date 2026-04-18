"""Smoke tests for the public atomview API."""

from pathlib import Path

import py3Dmol
import pytest

from atomview import load_structure, to_cif_string, view

DATA_DIR = Path(__file__).resolve().parent / "data"
CIF_PATHS = sorted(DATA_DIR.glob("*.cif"))
PRIMARY_STRUCTURE_PATH = CIF_PATHS[0]
XYZ_DIMENSIONS = 3
EXPECTED_CIF_COUNT = 11


def test_test_data_contains_expected_cif_fixtures() -> None:
    """Keep the bundled structure fixtures discoverable and complete."""
    assert len(CIF_PATHS) == EXPECTED_CIF_COUNT


@pytest.mark.parametrize("structure_path", CIF_PATHS, ids=lambda path: path.stem.upper())
def test_load_structure_from_path(structure_path: Path) -> None:
    """Load every bundled mmCIF structure from disk."""
    structure = load_structure(structure_path)

    assert structure.array_length() > 0
    assert structure.coord.shape[1] == XYZ_DIMENSIONS


def test_load_structure_from_cif_string() -> None:
    """Load a structure from an in-memory CIF string."""
    cif_text = PRIMARY_STRUCTURE_PATH.read_text()
    structure = load_structure(cif_text)

    assert structure.array_length() > 0
    assert structure.res_name[0]


def test_to_cif_string_returns_mmcif_text() -> None:
    """Serialize a structure back to mmCIF text."""
    structure = load_structure(PRIMARY_STRUCTURE_PATH)
    cif_text = to_cif_string(structure, id="test_structure")

    assert cif_text.startswith("data_test_structure")
    assert "_atom_site.Cartn_x" in cif_text


def test_view_smoke_for_all_test_cifs() -> None:
    """Build a viewer for every bundled CIF without crashing."""
    for structure_path in CIF_PATHS:
        structure = load_structure(structure_path)
        viewer = view(structure, show_surface=False, show_hover=False)

        assert isinstance(viewer, py3Dmol.view), structure_path.name
