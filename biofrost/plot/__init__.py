import logging
import matplotlib.pyplot as plt
from .utils import set_bar_width, set_violin_alpha, customize_cmap, get_value_colors, confidence_ellipse
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR) # Supress matplotlib font manager errors

__all__ = [
    "set_bar_width",
    "set_violin_alpha",
    "customize_cmap",
    "get_value_colors",
    "confidence_ellipse",
]

# Plotting parameters
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['axes.unicode_minus'] = False

plt.rcParams['lines.linewidth'] = .75

plt.rcParams['axes.linewidth'] = .75
plt.rcParams['grid.linewidth'] = .75
plt.rcParams['axes.labelsize'] = 7
plt.rcParams['axes.titlesize'] = 7

plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
plt.rcParams['xtick.major.size'] = 2
plt.rcParams['ytick.major.size'] = 2
plt.rcParams['xtick.major.width'] = .75
plt.rcParams['ytick.major.width'] = .75
plt.rcParams['xtick.minor.size'] = 1
plt.rcParams['ytick.minor.size'] = 1
plt.rcParams['xtick.minor.width'] = .75
plt.rcParams['ytick.minor.width'] = .75

plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 7
plt.rcParams['legend.title_fontsize'] = 7

plt.rcParams['figure.titlesize'] = 7
plt.rcParams['figure.figsize'] = 1.25, 1.25
plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.bbox'] = "tight"