import sys

import lib.plot

file = sys.argv[1]
id = sys.argv[2]
color = sys.argv[3]

lib.plot.plot_tss(file, color, "Distance from TSS (kb)", id + " count")
