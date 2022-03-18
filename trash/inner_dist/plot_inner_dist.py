#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 12:24:32 2018

@author: antony
"""

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

import sys
sys.path.append('/ifs/scratch/cancer/Lab_RDF/abh2138/scripts/python/')
import libplot

name=sys.argv[1]

libplot.setup()

df = pd.read_csv('inner_dist.txt', sep='\t', header=0)

df2 = df[(df['inner_dist'] < 1000)]

df3 = df[(df['inner_dist'] > 1000)]

fig, ax = libplot.newfig(w=8)
sns.distplot(df2, bins=100, ax=ax)
ax.set_xlim([-200, 500])
ax.set_xlabel('Inner distance (bp)')

libplot.savefig(fig, '{}_inner_dist_lt_1000_hist.pdf'.format(name))
plt.close(fig)

fig, ax = libplot.newfig(w=8)
sns.distplot(df3, bins=100, ax=ax)
ax.set_xlim([1000, 400000])
ax.set_xlabel('Inner distance (bp)')

libplot.savefig(fig, '{}_inner_dist_gt_1000_hist.pdf'.format(name))

plt.close(fig)
