# -*- coding: utf-8 -*-
"""
Filter peaks that are promoters only

@author: Antony Holmes
"""


import sys

import lib.peaks

file = sys.argv[1]

lib.peaks.filter_peaks(file, -2000, 1000)
