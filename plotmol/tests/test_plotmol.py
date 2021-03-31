import pytest
from bokeh.plotting import Figure

from plotmol import plotmol
from plotmol.plotmol import InputSizeError, default_tooltip_template


@pytest.fixture()
def bokeh_figure():
    return Figure(tooltips=default_tooltip_template())


def test_scatter(bokeh_figure):
    # It's difficult to test this as it's a graphical function. For now make sure
    # it doesn't error and something is produced.

    # Plot the data as an interactive scatter plot
    plotmol.scatter(
        bokeh_figure,
        x=[0.0, 1.0],
        y=[0.0, 1.0],
        smiles=["C", "CCO"],
        marker="x",
        marker_size=15,
        legend_label="Series A",
    )


def test_scatter_tagged_atoms(bokeh_figure):
    # again difficult to test but make sure a valid plot is made
    plotmol.scatter(
        bokeh_figure,
        x=[0.0],
        y=[0.0],
        smiles=["[H][c]1[n][c]([H])[c:2]([O:3][C:4]([H])([H])[H])[c:1]([H])[c]1[H]"],
        marker="x",
        marker_size=15,
        legend_label="Series A",
    )


def test_scatter_invalid_input_size(bokeh_figure):

    with pytest.raises(InputSizeError):
        plotmol.scatter(bokeh_figure, x=[0.0], y=[0.0, 1.0], smiles=["C", "CCO"])
