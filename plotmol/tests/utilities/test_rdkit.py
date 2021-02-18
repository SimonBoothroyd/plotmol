from plotmol.utilities.rdkit import smiles_to_svg


def test_smiles_to_svg():
    # It's difficult to test this as it's a graphical function. For now make sure
    # it doesn't error and something is produced.
    output = smiles_to_svg("CO")
    assert len(output) > 0
