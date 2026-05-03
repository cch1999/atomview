# atomview

A lightweight, standalone molecular structure viewer for Jupyter notebooks. Visualize [biotite](https://www.biotite-python.org/) `AtomArray` structures interactively using [py3Dmol](https://3dmol.csb.pitt.edu/) (3Dmol.js) — with chain-coloured cartoons, ligand sticks, ion spheres, surfaces, and hover labels out of the box.

**atomview** extracts and refactors the excellent `view()` visualisation utility from [atomworks](https://github.com/RosettaCommons/atomworks) (`atomworks.io.utils.visualize`) into a focused, minimal package. If all you need is a quick way to render structures in a notebook, you no longer have to pull in atomworks and its full dependency tree (RDKit, ML tooling, dataset infrastructure, etc.).

## Installation

```bash
pip install atomview
```

Requires Python >= 3.12.

## Quickstart

```python
from atomview import load_structure, view

# Load a structure from a local mmCIF file
structure = load_structure("4hhb.cif")

# Render it — that's it
view(structure)
```

## Features

- **One-liner visualisation** of any biotite `AtomArray` or `AtomArrayStack`
- **Automatic chain colouring** with a curated 18-colour palette
- **Smart representation** — proteins and nucleic acids get cartoon + outline sticks; ligands get element-coloured sticks; metal ions get spheres
- **Hover labels** showing chain, residue, atom name, and index
- **Optional VDW surface** overlay
- **Zoom to selection** — focus on a specific chain, residue, or atom
- **Solvent & crystallisation aid filtering** (SO4, GOL, EDO, PO4, etc.) enabled by default
- **Lightweight** — only depends on `biotite`, `py3Dmol`, `numpy`, and `matplotlib`

## Gallery

| Default protein view | Surface overlay |
|---|---|
| ![Chain-coloured cartoon rendering of hemoglobin](docs/images/default-protein.png) | ![Hemoglobin with translucent molecular surface](docs/images/surface-overlay.png) |
| One-line rendering with automatic chain colouring. | Optional translucent VDW surface for shape context. |

| Ligand focus | Mixed complex |
|---|---|
| ![4TOS protein-ligand rendering with ligand sticks](docs/images/ligand-focus.png) | ![2CV5 protein-DNA complex rendering](docs/images/mixed-complex.png) |
| Zoom to a ligand or binding-site selection. | Proteins, nucleic acids, ligands, and ions styled together. |

## API

### `view()`

The main entry point. Returns a `py3Dmol.view` object that renders inline in Jupyter.

```python
from atomview import view

v = view(
    structure,
    show_cartoon=True,            # cartoon for polymers
    show_surface=False,           # VDW surface overlay
    show_hover=True,              # hover labels
    hide_solvent=True,            # remove water
    hide_crystallization_aids=True,
    zoom_to_selection=None,       # e.g. {"chain": "A", "resi": 35}
    width=600,
    height=400,
    colors=None,                  # custom colour list, or use defaults
)
```

## Why not just use atomworks?

[atomworks](https://github.com/RosettaCommons/atomworks) is a fantastic, comprehensive framework for biomolecular modelling maintained by RosettaCommons. But it's a big library — it bundles I/O, transforms, ML dataset tooling, RDKit integration, and more. If you just want to call `view(structure)` in a notebook, that's a lot of overhead.

**atomview** gives you exactly that one thing: a great molecular viewer with sensible defaults, in a package that installs in seconds.

## License

MIT
