# -*- coding: utf-8 -*-
"""
Display peaks in a gene oriented fashion

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import lib.human.expression


file = sys.argv[1]

expression = lib.human.expression.PeakExpression()

expression.annotate(file)
