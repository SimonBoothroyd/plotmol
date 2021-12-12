smiles = [
    '[C:1]([N+:2]([C:3]([H:24])([H:25])[H:26])([C:4]([H:27])([H:28])[H:29])[c:5]1[c:6]([H:30])[c:7](-[c:14]2[c:15](['
    'H:34])[c:16]([H:35])[c:17]([O-:20])[c:18]([H:36])[c:19]2[H:37])[c:8]([H:31])[c:9]2[c:10]1[C:11]([H:32])=[C:12](['
    'H:33])[O:13]2)([H:21])([H:22])[H:23]',
    '[C:1]([N+:2]([C:3]([H:24])([H:25])[H:26])([C:4]([H:27])([H:28])[H:29])[c:5]1[c:6]([H:30])[c:7](-[c:14]2[c:15](['
    'H:34])[c:16]([H:35])[c:17]([S-:20])[c:18]([H:36])[c:19]2[H:37])[c:8]([H:31])[c:9]2[c:10]1[C:11]([H:32])=[C:12](['
    'H:33])[O:13]2)([H:21])([H:22])[H:23]',
    '[C:1]([N+:2]([C:3]([H:23])([H:24])[H:25])([C:4]([H:26])([H:27])[H:28])[c:5]1[c:6]([H:29])[c:7](-[c:14]2[c:15](['
    'H:33])[c:16]([H:34])[c:17]([H:35])[c:18]([H:36])[c:19]2[H:37])[c:8]([H:30])[c:9]2[c:10]1[C:11]([H:31])=[C:12](['
    'H:32])[O:13]2)([H:20])([H:21])[H:22]',
]

points = [0, 1, 2]

smarts = "[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]"

from plotmol.styles import MoleculeStyle

molecule_style = MoleculeStyle(
    image_width=200,
    image_height=200,

    highlight_tagged_atoms=True,
    highlight_tagged_bonds=True,
    substruct_smarts=smarts,
)

from bokeh.plotting import Figure
from plotmol.plotmol import default_tooltip_template

figure = Figure(
    # Required to show the molecule structure pop-up.
    tooltips=default_tooltip_template(),

    # Plot options. See the bokeh documentation for more information.
    title="Biphenyls sample",
    x_axis_label="X-axis Label",
    y_axis_label="Y-axis Label",

    plot_width=800,
    plot_height=600,
)

from bokeh.palettes import Spectral4
import plotmol

# Define a simple color palette.
palette = Spectral4
plotmol.scatter(figure,
                        x = points,
                        y = points,
                        smiles = smiles,
                        marker="o",
                        marker_size=15,
                        marker_color=palette[1],
                        molecule_style=molecule_style,
                        )

from bokeh.plotting import show

show(figure)