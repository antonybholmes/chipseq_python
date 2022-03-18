# -*- coding: utf-8 -*-
"""
Filter peaks that are promoters only

@author: Antony Holmes
"""


import sys

import pychipseq.peaks

file = sys.argv[1]

pychipseq.peaks.filter_peaks(file, -2000, 1000)
