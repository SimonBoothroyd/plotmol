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
https://github.com/SimonBoothroyd/plotmol/examples/scatter-plot.ipynb)*. 

In order to use the plotting functionality within a Jupyter notebook we first need to configure `bokeh` to run in 
notebook mode:

```python
from bokeh.io import output_notebook
output_notebook()
```

Next, similar to with `matplotlib`, we create a figure that we will render to.

```python
from bokeh.plotting import Figure
from plotmol.plotmol import default_tooltip_template

figure = Figure(
    # Required to show the molecule structure pop-up.
    tooltips=default_tooltip_template(),
    
    # Plot options. See the bokeh documentation for more information.
    title="My Dummy Data",
    x_axis_label="Dummy X",
    y_axis_label="Dummy Y",
    
    plot_width=500,
    plot_height=500,
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

### Copyright

Copyright (c) 2021, Simon Boothroyd

#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
