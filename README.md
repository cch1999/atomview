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
from atomview.core import view
from atomview.utils import _load_structure

# Load a structure from a local mmCIF file
structure = _load_structure("4hhb.cif")

# Render it — that's it
view(structure)
```

![3D viewer rendering of hemoglobin (4HHB)](https://img.shields.io/badge/py3Dmol-interactive-blue)

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

| Default protein view | Ligand focus |
|---|---|
| ![Chain-coloured cartoon rendering of hemoglobin](docs/images/default-protein.png) | ![Zoomed protein-ligand rendering with ligand sticks](docs/images/ligand-focus.png) |
| One-line rendering with automatic chain colouring. | Zoom to a ligand or binding-site selection. |

| Surface overlay | Mixed complex |
|---|---|
| ![Hemoglobin with translucent molecular surface](docs/images/surface-overlay.png) | ![Mixed biomolecular complex with ligand sticks](docs/images/mixed-complex.png) |
| Optional translucent VDW surface for shape context. | Proteins, nucleic acids, ligands, and ions styled together. |

## API

### `view()`

The main entry point. Returns a `py3Dmol.view` object that renders inline in Jupyter.

```python
from atomview.core import view

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

**Key parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `structure` | `AtomArray \| AtomArrayStack` | *required* | Biotite structure to visualise |
| `zoom_to_selection` | `dict \| None` | `None` | Zoom to `{"serial": 35}`, `{"chain": "A", "resi": 35}`, or `{"chain": "C"}` |
| `show_cartoon` | `bool` | `True` | Cartoon representation for proteins/nucleic acids |
| `show_surface` | `bool` | `True` | Semi-transparent VDW surface |
| `show_hover` | `bool` | `True` | Hover labels with atom details |
| `hide_solvent` | `bool` | `True` | Remove water molecules |
| `hide_crystallization_aids` | `bool` | `True` | Remove common crystallisation artefacts |
| `width` / `height` | `int` | `600` / `400` | Viewer dimensions in pixels |

### `to_cif_string()`

Convert a biotite `AtomArray` to a CIF-formatted string (useful for passing structures to other viewers or tools).

```python
from atomview.utils import to_cif_string

cif_str = to_cif_string(structure, id="my_protein")
```

### `_load_structure()`

Load a structure from a `.cif` or `.bcif` file, a file path string, or a `StringIO`/`BytesIO` buffer.

```python
from atomview.utils import _load_structure

structure = _load_structure("7gaw.cif")
```

## Why not just use atomworks?

[atomworks](https://github.com/RosettaCommons/atomworks) is a fantastic, comprehensive framework for biomolecular modelling maintained by RosettaCommons. But it's a big library — it bundles I/O, transforms, ML dataset tooling, RDKit integration, and more. If you just want to call `view(structure)` in a notebook, that's a lot of overhead.

**atomview** gives you exactly that one thing: a great molecular viewer with sensible defaults, in a package that installs in seconds.

## License

MIT
