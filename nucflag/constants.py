import re


PLOT_FONT_SIZE = 16
PLOT_HEIGHT = 6
PLOT_WIDTH = 16
PLOT_DPI = 600
PLOT_YLIM = 100
# TODO: Dunno if this should be configurable?
THR_MISJOIN_VALLEY_HEIGHT_PERC_BELOW = 0.000001
RGX_REGION = re.compile(r"(.+):(\d+)-(\d+)")
