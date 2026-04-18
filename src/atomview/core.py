"""Utilities to visualize atom arrays with py3Dmol in Jupyter notebooks."""

from __future__ import annotations

import warnings
from itertools import cycle

import biotite.structure as struc
import numpy as np
import py3Dmol
from biotite.structure import AtomArray, AtomArrayStack

from atomview.constants import CRYSTALLIZATION_AIDS, METAL_ELEMENTS, VIEWER_COLORS
from atomview.utils import reassign_chain_ids, to_cif_string

__all__ = ["view"]


def view(  # noqa: PLR0913
    structure: AtomArray | AtomArrayStack,
    *,
    zoom_to_selection: dict[str, int | str] | None = None,
    show_hover: bool = True,
    show_unoccupied: bool = False,
    show_cartoon: bool = True,
    show_surface: bool = True,
    width: int = 600,
    height: int = 400,
    ligand_linewidth: float = 0.2,
    polymer_sidechain_linewidth: float = 0.05,
    min_polymer_size: int = 1,
    hide_solvent: bool = True,
    hide_crystallization_aids: bool = True,
    colors: list[str] | None = None,
) -> py3Dmol.view:
    """Render a Biotite structure inline with py3Dmol.

    Args:
        structure: Structure to visualize.
        zoom_to_selection: Optional py3Dmol selection dict to zoom into.
        show_hover: Enable hover labels with atom details.
        show_unoccupied: Keep atoms with zero occupancy.
        show_cartoon: Show cartoon representation for polymers.
        show_surface: Add a semi-transparent VDW surface.
        width: Viewer width in pixels.
        height: Viewer height in pixels.
        ligand_linewidth: Stick radius for ligands.
        polymer_sidechain_linewidth: Stick radius for polymer sidechains.
        min_polymer_size: Minimum chain size to classify as polymer.
        hide_solvent: Remove solvent molecules before rendering.
        hide_crystallization_aids: Remove common crystallization additives.
        colors: Optional list of per-chain colors.

    Returns:
        A ``py3Dmol.view`` object suitable for Jupyter display.
    """
    if isinstance(structure, AtomArrayStack):
        warnings.warn("AtomArrayStack is not supported; using the first model.", stacklevel=2)
        structure = structure[0]

    viewer = py3Dmol.view(width=width, height=height)

    if not show_unoccupied and "occupancy" in structure.get_annotation_categories():
        structure = structure[structure.occupancy > 0]

    if hide_solvent:
        structure = structure[~struc.filter_solvent(structure)]
    if hide_crystallization_aids:
        structure = structure[~np.isin(structure.res_name, CRYSTALLIZATION_AIDS)]

    structure = reassign_chain_ids(structure)
    viewer.addModel(
        to_cif_string(structure, _allow_ambiguous_bond_annotations=True),
        "structure",
        format="mmcif",
    )

    palette = colors or VIEWER_COLORS
    for chain_id, color in zip(struc.get_chains(structure), cycle(palette)):
        chain = structure[structure.chain_id == chain_id]
        is_protein = np.all(
            struc.filter_polymer(chain, pol_type="peptide", min_size=min_polymer_size)
            & struc.filter_amino_acids(chain)
        )
        is_nucleic = np.all(
            struc.filter_polymer(chain, pol_type="nucleotide", min_size=min_polymer_size)
        )
        is_ion = len(chain.element) > 0 and np.all(np.isin(chain.element, METAL_ELEMENTS))

        if is_protein or is_nucleic:
            style = {"stick": {"radius": polymer_sidechain_linewidth, "style": "outline"}}
            if show_cartoon:
                style["cartoon"] = {"color": color, "arrows": True}
            viewer.setStyle({"chain": chain_id}, style)
        elif is_ion:
            viewer.setStyle({"chain": chain_id}, {"sphere": {"scale": 0.8}})
        else:
            viewer.setStyle(
                {"chain": chain_id, "elem": "C"},
                {"stick": {"color": color, "radius": ligand_linewidth}},
            )
            viewer.setStyle(
                {"chain": chain_id, "not": {"elem": "C"}},
                {"stick": {"colorscheme": "element", "radius": ligand_linewidth}},
            )

    if show_surface:
        viewer.addSurface(py3Dmol.VDW, {"opacity": 0.4, "color": "gray"})

    if show_hover:
        js_script = """function(atom,viewer) {
                    if(!atom.label) {
                        atom.label = viewer.addLabel(
                            atom.chain + ':' +
                            atom.resn + '(' + atom.resi + '):' +
                            atom.atom + '(idx' + atom.serial + ')',
                            {position: atom, backgroundColor:"white", fontColor:"black"}
                        );
                    }
                }"""
        viewer.setHoverable(
            {},
            True,
            js_script,
            """function(atom,viewer) {
                    if(atom.label) {
                        viewer.removeLabel(atom.label);
                        delete atom.label;
                    }
                    }""",
        )

    viewer.zoomTo()
    if zoom_to_selection is not None:
        viewer.zoomTo(zoom_to_selection)

    return viewer
