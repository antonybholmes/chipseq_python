# -*- coding: utf-8 -*-
"""
Duplicates ad-hoc locations

This forms the standard peak file

Created on Tue Jul 22 15:43:36 2014

@author: Antony Holmes
"""


import sys

import peaks

file = sys.argv[1]

peaks.duplicate_peaks("Location", file)
