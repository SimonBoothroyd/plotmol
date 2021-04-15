import functools

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from plotmol.styles import MoleculeStyle


@functools.lru_cache(1024)
def smiles_to_svg(
    smiles: str,
    style: MoleculeStyle,
) -> str:
    """Renders a 2D representation of a molecule based on its SMILES representation as
    an SVG string.

    Parameters
    ----------
    smiles
        The SMILES pattern.
    style
        Options which control how the structure should be rendered as an image.

    Returns
    -------
        The 2D SVG representation.
    """

    # Parse the SMILES into an RDKit molecule
    smiles_parser = Chem.rdmolfiles.SmilesParserParams()
    smiles_parser.removeHs = False

    rdkit_molecule = Chem.MolFromSmiles(smiles, smiles_parser)

    # look for any tagged atom indices
    tagged_atoms = (
        []
        if not style.highlight_tagged_atoms
        else [
            atom.GetIdx()
            for atom in rdkit_molecule.GetAtoms()
            if atom.GetAtomMapNum() != 0
        ]
    )
    tagged_bonds = (
        []
        if not style.highlight_tagged_bonds
        else [
            bond.GetIdx()
            for bond in rdkit_molecule.GetBonds()
            if bond.GetBeginAtom().GetAtomMapNum() != 0
            and bond.GetEndAtom().GetAtomMapNum() != 0
        ]
    )

    # Generate a set of 2D coordinates.
    if not rdkit_molecule.GetNumConformers():
        Chem.rdDepictor.Compute2DCoords(rdkit_molecule)

    drawer = rdMolDraw2D.MolDraw2DSVG(style.image_width, style.image_height)
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        rdkit_molecule,
        highlightAtoms=tagged_atoms,
        highlightBonds=tagged_bonds,
    )
    drawer.FinishDrawing()

    svg_content = drawer.GetDrawingText()
    return svg_content
