import base64
from typing import Any, Callable, Dict, List, Literal, Optional, Union

from bokeh.models import ColumnDataSource
from bokeh.plotting import Figure

from plotmol.styles import MoleculeStyle
from plotmol.utilities.rdkit import smiles_to_svg

DataType = Union[float, int, str]

Label = Dict[str, str]
Marker = Literal["x", "o"]

MoleculeToImageFunction = Callable[[str, MoleculeStyle], str]


class InputSizeError(ValueError):
    """An error raised when inputs of different sizes are provided to a plot function."""

    def __init__(self, sizes: Dict[str, int]):

        sizes_str = ", ".join(f"{label} ({size})" for label, size in sizes.items())

        super(InputSizeError, self).__init__(
            f"The input data arrays must all have the same length: {sizes_str}."
        )


def default_tooltip_template() -> str:
    """Returns the default html based template to use for the figure hover pop-up."""

    return """
    <div>
        <div>
            <img src="@image" ></img>
        </div>
    </div>
"""


def scatter(
    figure: Figure,
    x: List[DataType],
    y: List[DataType],
    smiles: List[str],
    legend_label: Optional[str] = None,
    marker: Optional[str] = None,
    marker_size: Optional[int] = None,
    marker_color: Optional[str] = None,
    molecule_style: Optional[MoleculeStyle] = None,
    molecule_to_image_function: Optional[MoleculeToImageFunction] = None,
    **kwargs: Dict[str, Any],
):
    """Adds a scatter series to a bokeh figure which will show the molecular
    structure associated with a data point when the user hovers over it.

    Args:
        figure: The bokeh figure to plot the scatter data on.
        x: An array of the x values to plot.
        y: An array of the y values to plot.
        smiles: An array of the SMILES patterns associated with each (x, y) pair.
        legend_label: The label to show in the legend for this data series.
        marker: The marker style.
        marker_size: The size of the marker to draw.
        marker_color: The marker color.
        molecule_style: Options which control how the 2D structure which is shown when a
            data point is hovered over should be rendered.
        molecule_to_image_function: The function which should be used to render the
            molecule as a 2D SVG image. This function should accept a string SMILES
            pattern and a molecule style object and return a valid SVG string (e.g.
            ``"<svg>...</svg>"``). By default the ``plotmol.rdkit.smiles_to_svg``
            function will be used.
        kwargs: Extra keyword arguments to pass to the underlying bokeh ``scatter``
            function.
    """

    # Set the default mutable inputs.
    molecule_style = molecule_style if molecule_style is not None else MoleculeStyle()

    molecule_to_image_function = (
        molecule_to_image_function
        if molecule_to_image_function is not None
        else smiles_to_svg
    )

    # Validate the sizes of the input arrays.
    data_sizes = {"x": len(x), "y": len(y), "smiles": len(smiles)}

    if len({*data_sizes.values()}) != 1:
        raise InputSizeError(data_sizes)

    # Generate an image for each SMILES pattern.
    raw_images = [
        base64.b64encode(
            molecule_to_image_function(smiles_pattern, molecule_style).encode()
        ).decode()
        for smiles_pattern in smiles
    ]
    images = [f"data:image/svg+xml;base64,{raw_image}" for raw_image in raw_images]

    # Create a custom data source.
    source = ColumnDataSource(data={"x": x, "y": y, "smiles": smiles, "image": images})

    # Add the scatter data.
    scatter_kwargs = {**kwargs}

    if marker is not None:
        scatter_kwargs["marker"] = marker
    if marker_size is not None:
        scatter_kwargs["size"] = marker_size
    if marker_color is not None:
        scatter_kwargs["color"] = marker_color
    if legend_label is not None:
        scatter_kwargs["legend_label"] = legend_label

    figure.scatter(x="x", y="y", source=source, **scatter_kwargs)
