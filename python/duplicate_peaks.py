# -*- coding: utf-8 -*-
"""
Duplicates rows of an annotation file to create unique gene and mir rows

This forms the standard peak file

Created on Tue Jul 22 15:43:36 2014

@author: Antony Holmes
"""


import sys

import pychipseq.peaks

file = sys.argv[1]

pychipseq.peaks.duplicate_peaks("Peak", file)
