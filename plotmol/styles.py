from dataclasses import dataclass


@dataclass(frozen=True)
class MoleculeStyle:
    """Options to control the appearance of the 2D structure of the drawn molecule that
    appears when a data point is hovered over.

    Attributes:

        image_width: The width [px] to make the image of the 2D structure.

        image_height: The height [px] to make the image of the 2D structure.

        highlight_tagged_atoms: Whether to highlight atoms which have been tagged with
            a map index, e.g. whether the oxygen atom should be highlighted for
            ``"C[O:1]"``.

        highlight_tagged_bonds: Whether to highlight bonds between atoms which have
            been tagged with a map index, e.g. whether the carbon-oxygen bond should be
            highlighted for ``"[C:1][O:2]"``.

        substruct_smarts: to highlight atoms in molecules based on smarts pattern

        show_all_hydrogens: to turn on/off displaying hydrogens in the final image

    """

    image_width: int = 200
    image_height: int = 200

    highlight_tagged_atoms: bool = True
    highlight_tagged_bonds: bool = True

    substruct_smarts: str = None
    show_all_hydrogens: bool = False
    turn_off_atom_index: bool = False
