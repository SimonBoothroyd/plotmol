import plotmol.styles
import plotmol.utilities


def test_smiles_to_svg():
    # It's difficult to test this as it's a graphical function. For now make sure
    # it doesn't error and something is produced.
    output = plotmol.utilities.smiles_to_svg("CO", plotmol.styles.MoleculeStyle())
    assert len(output) > 0
