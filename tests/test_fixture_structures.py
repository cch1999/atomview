"""Smoke tests using bundled mmCIF/PDB fixtures under tests/data/."""

from __future__ import annotations

from pathlib import Path

from atomview import load_structure, view

_DATA_DIR = Path(__file__).resolve().parent / "data"
_STRUCTURE_SUFFIXES = frozenset({".cif", ".pdb"})


def _fixture_structure_paths() -> list[Path]:
    if not _DATA_DIR.is_dir():
        return []
    return sorted(
        p for p in _DATA_DIR.iterdir() if p.is_file() and p.suffix.lower() in _STRUCTURE_SUFFIXES
    )


def test_view_does_not_crash_on_fixture_structures() -> None:
    """Building the full py3Dmol view pipeline must not raise for bundled fixtures."""
    paths = _fixture_structure_paths()
    assert paths, f"expected structure files under {_DATA_DIR}"
    for structure_path in paths:
        structure = load_structure(structure_path)
        view(structure, show_hover=False, show_surface=False)
