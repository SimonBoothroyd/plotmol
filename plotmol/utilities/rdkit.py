import functools
import itertools

from rdkit.Chem.Draw import rdMolDraw2D

from plotmol.styles import MoleculeStyle
from rdkit import Chem


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
    if not style.show_all_hydrogens:
        # updateExplicitCount: Keep a record of the hydrogens we remove.
        # This is used in visualization to distinguish eg radicals from normal species
        rdkit_molecule = Chem.rdmolops.RemoveHs(rdkit_molecule, updateExplicitCount=True)
    # highlight substructure
    substruct_smarts = style.substruct_smarts
    # list_of_atoms = style.list_of_atoms

    if substruct_smarts:
        atom_matches = list(itertools.chain(*rdkit_molecule.GetSubstructMatches(Chem.MolFromSmarts(substruct_smarts))))
    # elif list_of_atoms:
    #     atom_matches = list_of_atoms
    else:
        atom_matches = []
    # look for any tagged atom indices
    tagged_atoms = atom_matches

    tagged_bonds = (
        []
        if not style.highlight_tagged_bonds
        else [
            bond.GetIdx()
            for bond in rdkit_molecule.GetBonds()
            if bond.GetBeginAtom().GetIdx() in atom_matches
               and bond.GetEndAtom().GetIdx() in atom_matches
        ]
    )

    for atom in rdkit_molecule.GetAtoms():
        atom.SetAtomMapNum(0)

    Chem.Draw.rdDepictor.SetPreferCoordGen(True)
    Chem.Draw.rdDepictor.Compute2DCoords(rdkit_molecule)

    rdkit_molecule = rdMolDraw2D.PrepareMolForDrawing(rdkit_molecule)

    drawer = rdMolDraw2D.MolDraw2DSVG(style.image_width, style.image_height)
    drawer.drawOptions().setHighlightColour((0.3, 1.0, 0.5))
    drawer.drawOptions().addAtomIndices = False
    drawer.drawOptions().addBondIndices = False
    drawer.DrawMolecule(rdkit_molecule,
                        highlightAtoms=tagged_atoms,
                        highlightBonds=tagged_bonds,
                        )
    drawer.FinishDrawing()

    svg_content = drawer.GetDrawingText()
    return svg_content
