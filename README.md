# Interactive Plots with Molecular Annotations

[![tests](https://github.com/SimonBoothroyd/plotmol/workflows/CI/badge.svg?branch=main)](https://github.com/SimonBoothroyd/plotmol/actions?query=workflow%3ACI)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/SimonBoothroyd/plotmol.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/SimonBoothroyd/plotmol/context:python)
[![codecov](https://codecov.io/gh/SimonBoothroyd/plotmol/branch/main/graph/badge.svg?token=Aa8STE8WBZ)](https://codecov.io/gh/SimonBoothroyd/plotmol)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This framework provides a very simple set of functionalities for creating interactive plots that are annotated with 
molecular structures using the [`bokeh` library](https://docs.bokeh.org/en/latest/index.html).

## Installation

The package and its dependencies can be installed using the `conda` package manager:

```shell
conda install -c conda-forge -c simonboothroyd plotmol
```

## Getting Started

*We recommend viewing the getting started example in a Jupyter notebook. [This full example can be found here](
https://github.com/SimonBoothroyd/plotmol/blob/main/examples/scatter-plot.ipynb)*. 

In order to use the plotting functionality within a Jupyter notebook we first need to configure `bokeh` to run in 
notebook mode:

```python
from bokeh.io import output_notebook
output_notebook()
```

Next, similar to with `matplotlib`, we create a figure that we will render to.

```python
from bokeh.plotting import figure
from plotmol.plotmol import default_tooltip_template

figure = figure(
    # Required to show the molecule structure pop-up.
    tooltips=default_tooltip_template(),
    
    # Plot options. See the bokeh documentation for more information.
    title="My Dummy Data",
    x_axis_label="Dummy X",
    y_axis_label="Dummy Y",
    
    width=500,
    height=500,
)
```

Data can then be rendered on the figure similar to `matplotlib` using the `plotmol.scatter` function:

```python
from bokeh.palettes import Spectral4
import plotmol

# Define a simple color palette.
palette = Spectral4

# Plot some dummy data to the figure.
plotmol.scatter(
    figure,
    x=[0.0, 1.0, 2.0, 3.0],
    y=[0.0, 1.0, 2.0, 3.0],
    smiles=["C", "CO", "CC", "CCO"],
    marker="x",
    marker_size=15,
    marker_color=palette[0],
    legend_label="Series A"
)
plotmol.scatter(
    figure,
    x=[0.0, 1.0, 2.0, 3.0],
    y=[1.0, 2.0, 3.0, 4.0],
    smiles=["C=O", "CC=O", "COC", "CCCO"],
    marker="o",
    marker_size=15,
    marker_color=palette[1],
    legend_label="Series B"
)
```

We can optionally configure the legend which will be added to the plot. Here we move it to the top left corner
of the plot and enable the option to toggle data series when they are clicked in the plot legend:

```python
figure.legend.location = "top_left"
figure.legend.click_policy = "hide"
```

Finally, we show the figure. Hovering over each data point with the mouse should reveal a pop-up containing the 
2D molecular structure associated with a data point. 

```python
from bokeh.plotting import show
show(figure)
```

## Styling the Molecule Structure

The molecule structure which is shown when a data point is hovered over can be customised in two ways: by passing a
style object to the `scatter` function or by using a custom function to render the SVG image containing the 2D 
structure entirely.

### Changing the molecule style

The `scatter` function accepts a `molecule_style` argument which will control certain aspects of how the 2D
structure will be rendered:

```python
from plotmol.styles import MoleculeStyle

molecule_style = MoleculeStyle(
    image_width=200,
    image_height=200,
    
    highlight_tagged_atoms=True,
    highlight_tagged_bonds=True,
)

plotmol.scatter(
    figure,
    x=[0.0, 1.0, 2.0, 3.0],
    y=[1.0, 2.0, 3.0, 4.0],
    smiles=["C=O", "CC=O", "COC", "CCCO"],
    molecule_style=molecule_style
)
```

By default atoms which have been tagged with a map index, e.g. `"[C:1][C:2]"`, and the bonds between those atoms will
be highlighted. This behaviour can be disabled by setting `highlight_tagged_atoms` and / or `highlight_tagged_bonds`
to `False`.

### Using a custom render function

By default the 2D structure of the molecule associated with a data point is drawn using RDKit using the 
``plotmol.utilities.rdkit.smiles_to_svg`` function, however an entirely custom function can be passed to
the `scatter` function:

```python

def custom_image_function(smiles: str, style: MoleculeStyle) -> str:
    
    # Render the molecule to a SVG
    svg_contents = ...
    
    return svg_contents

plotmol.scatter(
    figure,
    x=[0.0, 1.0, 2.0, 3.0],
    y=[1.0, 2.0, 3.0, 4.0],
    smiles=["C=O", "CC=O", "COC", "CCCO"],
    molecule_to_image_function=custom_image_function
)

```

The function should take as arguments a SMILES pattern which defines the molecule that should be drawn, and a style
object which specifies options for how the molecule should be drawn.

## Customising the Tooltips

The tooltips that are displayed when a data point is hovered over can be easily customized by providing a custom 
tooltip template when creating a figure:

```python
custom_tooltip_template = """
<div>
    <div>
        <span>@title</span>
        <img src="@image" ></img>
    </div>
</div>
"""

figure = figure(
    tooltips=custom_tooltip_template,
    ...
)
```

To ensure that the molecular structure associated with a data point is shown the template must contain an image tag
which uses the `@image` key as its `src`:

```html
<img src="@image" ></img>
```

Extra data can be made available to the tooltip template by providing a ``custom_column_data`` to the ``scatter`` 
function:

```python
plotmol.scatter(
    figure,
    x=[0.0, 1.0, 2.0, 3.0],
    y=[0.0, 1.0, 2.0, 3.0],
    smiles=["C", "CO", "CC", "CCO"],
    marker="x",
    marker_size=15,
    marker_color=palette[0],
    legend_label="Series A",
    custom_column_data={"title": ["A", "B", "C", "D"]}
)
```

Each key in the data dictionary will be exposed to the tooltip template in the form `@key`. In this example the custom 
toolkit template can make use of the additional `@title` key. 

More information about custom tooltips can be found in the [Bokeh documentation](https://docs.bokeh.org/en/latest/docs/user_guide/tools.html#custom-tooltip)

### Copyright

Copyright (c) 2021, Simon Boothroyd

#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
