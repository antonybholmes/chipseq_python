# -*- coding: utf-8 -*-
"""
Center and expand locations to make regions for heat maps

@author: antony
"""

import sys
import lib.genomic

PADDING = 6000

locations_file = sys.argv[1]
padding = int(sys.argv[2])

locations, header = lib.genomic.load_locations(locations_file, True, padding, padding)

sys.stdout.write(header + "\n")

for l in locations:
  sys.stdout.write(l.to_string() + "\n")
