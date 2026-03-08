"""Utility functions to visualize atom arrays with py3Dmol in Jupyter notebooks."""

__all__ = ["view"]

import warnings
from itertools import cycle

from biotite.sequence.alphabet import Alphabet
import biotite.structure as struc
import numpy as np
import py3Dmol
from biotite.structure import AtomArray, AtomArrayStack
from biotite.structure.io import mol, pdb, pdbx

from atomview.constants import VIEWER_COLORS, METAL_ELEMENTS, CRYSTALLIZATION_AIDS
from atomview.utils import reassign_chain_ids, to_cif_string


def view(
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
    colors: list[str] = VIEWER_COLORS,
) -> py3Dmol.view:
    """Visualize an AtomArray structure using py3Dmol for display in jupyter notebooks.

    Args:
        - structure (AtomArray): The atomic structure to be visualized.
        - zoom_to_selection (dict[str, int | str] | None, optional): A dictionary specifying the
            selection to zoom into. Defaults to None. Here are some examples:
                - `{'serial': 35}` - will zoom to the atom with index 35 in the atom array
                - `{'chain': 'A', 'resi': 35}` - will zoom to the residue id 35 in chain A
                - `{'chain': 'C'} - will zoom to the entire chain C
            !WARNING! If the selection is wrong, the visualization will be empty.
        - show_hover (bool, optional): Whether to enable hover functionality to display atom details.
            Defaults to True.
        - show_unoccupied (bool, optional): Whether to show unoccupied atoms. Defaults to False.
        - show_cartoon (bool, optional): Whether to show the cartoon. Defaults to True.
        - show_surface (bool, optional): Whether to show the surface. Defaults to False.
        - width (int, optional): The width of the visualization window. Defaults to 400.
        - height (int, optional): The height of the visualization window. Defaults to 300.
        - ligand_linewidth (float, optional): The linewidth for ligand representation. Defaults to 0.2.
        - polymer_sidechain_linewidth (float, optional): The linewidth for polymer sidechain representation. Defaults to 0.05.
        - min_polymer_size (int, optional): The minimum size for a chain to be displayed as a polymer. Defaults to 1.
        - hide_solvent (bool, optional): Whether to hide solvent. Defaults to True.
        - hide_crystallization_aids (bool, optional): Whether to hide crystallization aids. Defaults to True.
        - colors (list[str], optional): A list of colors to cycle through for different chains. Defaults to IPD_PYMOL_COLORS.

    Returns:
        py3Dmol.view: The py3Dmol view object for the structure visualization.
    """
    if isinstance(structure, AtomArrayStack):
        warnings.warn("AtomArrayStack is not supported; using the first model.", stacklevel=2)
        structure = structure[0]

    # Initialize the py3Dmol view with specified width and height
    view = py3Dmol.view(width=width, height=height)

    # Handle unoccupied atoms
    if not show_unoccupied and ("occupancy" in structure.get_annotation_categories()):
        structure = structure[structure.occupancy > 0]

    if hide_solvent:
        structure = structure[~struc.filter_solvent(structure)]
    if hide_crystallization_aids:
        structure = structure[~np.isin(structure.res_name, CRYSTALLIZATION_AIDS)]

    # Reassign chain IDs to separate hetero and non-hetero atoms
    structure = reassign_chain_ids(structure)

    # Convert the structure to a temporary CIF string for interacting with py3Dmol
    _tmp_cif_str = to_cif_string(
        structure,
        _allow_ambiguous_bond_annotations=True,
    )
    # ... add the structure model to the view in mmCIF format
    view.addModel(_tmp_cif_str, "structure", format="mmcif")
    # Get the chain IDs from the structure
    chain_ids = struc.get_chains(structure)

    # Iterate over each chain and assign styles based on the type of polymer
    for chain_id, color in zip(chain_ids, cycle(colors)):
        is_protein = np.all(
            struc.filter_polymer(
                structure[structure.chain_id == chain_id], pol_type="peptide", min_size=min_polymer_size
            )
            & struc.filter_amino_acids(structure[structure.chain_id == chain_id])
        )
        print(struc.filter_polymer(
               structure[structure.chain_id == chain_id], pol_type="nucleotide", min_size=min_polymer_size
             ))
        print(struc.filter_nucleotides(structure[structure.chain_id == chain_id]))
        is_nucleic = np.all(
            struc.filter_polymer(
                structure[structure.chain_id == chain_id], pol_type="nucleotide", min_size=min_polymer_size
            )
            # struc.filter_nucleotides(structure[structure.chain_id == chain_id])
        )
        # Check if all atoms in chain are metals
        chain_elements = structure[structure.chain_id == chain_id].element
        is_ion = len(chain_elements) > 0 and np.all(np.isin(chain_elements, METAL_ELEMENTS))

        print(chain_id, is_protein, is_nucleic, is_ion)

        if is_protein or is_nucleic:
            # Apply protein or nucleic acid style
            style = {"stick": {"radius": polymer_sidechain_linewidth, "style": "outline"}}
            if show_cartoon:
                style["cartoon"] = {"color": color, "arrows": True}
            view.setStyle({"chain": chain_id}, style)

        # elif is_ion:
        #     view.setStyle(
        #         {"chain": chain_id},
        #         {"stick": {"radius": polymer_sidechain_linewidth, "style": "outline"}},
        #     )

        elif is_ion:
            # Apply ion style
            view.setStyle(
                {"chain": chain_id},
                {"sphere": {"scale": 0.8}},
            )
        else:
            # Apply ligand style
            # ... first, set the style for carbon atoms colored by chain
            view.setStyle(
                {"chain": chain_id, "elem": "C"},
                {"stick": {"color": color, "radius": ligand_linewidth}},
            )
            # ... then, set the style for all other atoms based on the element
            view.setStyle(
                {"chain": chain_id, "not": {"elem": "C"}},
                {"stick": {"colorscheme": "element", "radius": ligand_linewidth}},
            )

    if show_surface:
        view.addSurface(py3Dmol.VDW, {"opacity": 0.4, "color": "gray"})

    # Add hover functionality to display atom details on hover
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
        view.setHoverable(
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

    # Zoom to the entire structure or to a specific selection if provided
    view.zoomTo()
    if zoom_to_selection is not None:
        view.zoomTo(zoom_to_selection)

    return view