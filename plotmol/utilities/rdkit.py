#CORRECT FILE
import functools

import pickle

from openeye import oechem, oeomega
from openforcefield.topology import Molecule
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
    
def smiles2oemol(smiles):
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    #initialize omega
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetIncludeInput(True)
    omega.SetCanonOrder(True)
    omega.SetSampleHydrogens(True)
    omega.SetStrictStereo(False)
    omega.SetStrictAtomTypes(True)
    omega.SetIncludeInput(False)
    status = omega(mol)

    return (mol, status)

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
    
    oe_molecule, status = smiles2oemol(smiles)
    openff_molecule = Molecule.from_openeye(oe_molecule)
    rdkit_molecule = openff_molecule.to_rdkit()

    #rdkit_molecule = Chem.MolFromSmiles(smiles, smiles_parser)
    
    # Generate a set of 2D coordinates.
    if not rdkit_molecule.GetNumConformers():
        Chem.rdDepictor.Compute2DCoords(rdkit_molecule)

    drawer = rdMolDraw2D.MolDraw2DSVG(image_width, image_height)

    if rdkit_molecule.GetBondBetweenAtoms(torsion_indices[0], torsion_indices[1]) != None:
        bonds = [rdkit_molecule.GetBondBetweenAtoms(torsion_indices[0], torsion_indices[1]).GetIdx()]
    else:
        bond = []
        
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, rdkit_molecule, highlightBonds = bonds)
                
    drawer.FinishDrawing()

    svg_content = drawer.GetDrawingText()
    return svg_content
