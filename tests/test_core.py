"""Tests for the atomview py3Dmol viewer entry point."""

import biotite.structure as struc
import numpy as np
import pytest
from biotite.structure import AtomArray

from atomview import core


class FakeViewer:
    """Small py3Dmol.view stand-in that records viewer calls."""

    def __init__(self) -> None:
        self.models: list[tuple[str, str, str]] = []
        self.styles: list[tuple[dict, dict]] = []
        self.surfaces: list[tuple[object, dict]] = []
        self.hoverables: list[tuple[dict, bool, str, str]] = []
        self.zoom_calls: list[dict[str, int | str] | None] = []

    def addModel(self, cif: str, model_type: str, *, format: str) -> None:
        """Record model additions."""
        self.models.append((cif, model_type, format))

    def setStyle(self, selection: dict, style: dict) -> None:
        """Record style changes."""
        self.styles.append((selection, style))

    def addSurface(self, surface_type: object, style: dict) -> None:
        """Record surface additions."""
        self.surfaces.append((surface_type, style))

    def setHoverable(self, selection: dict, hoverable: bool, on_enter: str, on_leave: str) -> None:
        """Record hover setup."""
        self.hoverables.append((selection, hoverable, on_enter, on_leave))

    def zoomTo(self, selection: dict[str, int | str] | None = None) -> None:
        """Record zoom calls."""
        self.zoom_calls.append(selection)


def make_atom_array() -> AtomArray:
    """Create a minimal AtomArray for viewer tests."""
    structure = struc.AtomArray(2)
    structure.coord = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    structure.chain_id = np.array(["A", "A"])
    structure.res_id = np.array([1, 1])
    structure.res_name = np.array(["GLY", "GLY"])
    structure.atom_name = np.array(["N", "CA"])
    structure.element = np.array(["N", "C"])
    structure.hetero = np.array([False, False])
    return structure


def test_view_returns_viewer_and_honors_zoom_selection(monkeypatch: pytest.MonkeyPatch) -> None:
    """Viewer creation should add a model and zoom to an optional selection."""
    fake_viewer = FakeViewer()
    monkeypatch.setattr(core.py3Dmol, "view", lambda **_: fake_viewer)
    selection = {"chain": "A", "resi": 1}

    result = core.view(
        make_atom_array(),
        zoom_to_selection=selection,
        show_hover=False,
        show_surface=False,
    )

    assert result is fake_viewer
    assert fake_viewer.models[0][1:] == ("structure", "mmcif")
    assert fake_viewer.zoom_calls == [None, selection]


def test_view_warns_and_uses_first_model_for_atom_array_stack(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """AtomArrayStack inputs should warn and still produce a viewer."""
    fake_viewer = FakeViewer()
    monkeypatch.setattr(core.py3Dmol, "view", lambda **_: fake_viewer)
    structure = make_atom_array()
    stack = struc.stack([structure, structure.copy()])

    with pytest.warns(UserWarning, match="using the first model"):
        result = core.view(stack, show_hover=False, show_surface=False)

    assert result is fake_viewer
    assert len(fake_viewer.models) == 1
