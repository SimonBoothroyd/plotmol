import base64
import typing

import bokeh.models
import bokeh.plotting

import plotmol.styles
import plotmol.utilities

DataType = float | int | str

MoleculeToImageFunction = typing.Callable[[str, plotmol.styles.MoleculeStyle], str]


class InputSizeError(ValueError):
    """An error raised when inputs of different sizes are provided to a plot
    function."""

    def __init__(self, sizes: dict[str, int] | None = None):
        sizes = {} if sizes is None else sizes
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
    figure: bokeh.plotting.figure,
    x: list[DataType],
    y: list[DataType],
    smiles: list[str],
    legend_label: str | None = None,
    marker: str | None = None,
    marker_size: int | None = None,
    marker_color: str | None = None,
    molecule_style: plotmol.styles.MoleculeStyle | None = None,
    molecule_to_image_function: MoleculeToImageFunction | None = None,
    custom_column_data: dict[str, list[DataType]] | None = None,
    **kwargs: dict[str, typing.Any],
) -> bokeh.models.GlyphRenderer:
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
        custom_column_data: An optional dictionary of extra entries to include in the
            ``ColumnDataSource`` used to render the plot. These custom entries can be
            accessed from a figures tooltip and used, for example, to add extra outputs
            such as a tooltip title or additional properties associated with a data
            point.

            By default the data included are the ``x``, ``y``, and ``smiles`` lists.

            Each key will correspond to a '@key' that will be made available to
            the tooltip template (see `the Bokeh documentation
            <https://docs.bokeh.org/en/latest/docs/user_guide/tools.html#custom-tooltip>`_
            for more details) and each value must be a list of the corresponding values
            with a length equal to ``x`` and ``y``.
        kwargs: Extra keyword arguments to pass to the underlying bokeh ``scatter``
            function.
    """

    molecule_style = (
        molecule_style if molecule_style is not None else plotmol.styles.MoleculeStyle()
    )

    molecule_to_image_function = (
        molecule_to_image_function
        if molecule_to_image_function is not None
        else plotmol.utilities.smiles_to_svg
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
    column_data = {"x": x, "y": y, "smiles": smiles, "image": images}

    if custom_column_data is not None:
        invalid_sizes = {
            key: len(entry)
            for key, entry in custom_column_data.items()
            if len(entry) != len(x)
        }

        if len(invalid_sizes) > 0:
            raise InputSizeError(invalid_sizes)

        column_data.update(custom_column_data)

    source = bokeh.models.ColumnDataSource(data=column_data)

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

    return figure.scatter(x="x", y="y", source=source, **scatter_kwargs)
