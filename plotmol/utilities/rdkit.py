#CORRECT FILE
import functools

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


@functools.lru_cache(1024)
def smiles_to_svg(smiles: str, torsion_indices: (int, int), image_width: int = 200, image_height: int = 200) -> str:
    """Renders a 2D representation of a molecule based on its SMILES representation as
    an SVG string.

    Parameters
    ----------
    smiles
        The SMILES pattern.
    image_width
        The width to make the final SVG.
    image_height
        The height to make the final SVG.

    Returns
    -------
        The 2D SVG representation.
    """

    # Parse the SMILES into an RDKit molecule
    smiles_parser = Chem.rdmolfiles.SmilesParserParams()
    smiles_parser.removeHs = False

    rdkit_molecule = Chem.MolFromSmiles(smiles, smiles_parser)

    # Generate a set of 2D coordinates.
    if not rdkit_molecule.GetNumConformers():
        Chem.rdDepictor.Compute2DCoords(rdkit_molecule)

    drawer = rdMolDraw2D.MolDraw2DSVG(image_width, image_height)
    
    # Locate the bond to be highlighted using the central torsion indices for the mol
    if rdkit_molecule.GetBondBetweenAtoms(torsion_indices[0], torsion_indices[1]) != None:
        bond = [rdkit_molecule.GetBondBetweenAtoms(torsion_indices[0], torsion_indices[1]).GetIdx()]
    else:
        bond = []
        
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, rdkit_molecule, highlightBonds = bond)
                
    drawer.FinishDrawing()

    svg_content = drawer.GetDrawingText()
    return svg_content
