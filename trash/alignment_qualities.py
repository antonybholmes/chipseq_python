#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 09:33:33 2019

@author: antony
"""

import sys
import os
import pandas as pd
import re
import collections
import pandas as pd
import numpy as np
import time

dfq_tables = []

peak_counts = collections.defaultdict(lambda: collections.defaultdict(int))

for sample in os.listdir('.'):
    if not os.path.isdir(sample):
        continue
    
    print('Processing', sample)
    
    # find quality files
    
    for (root, dirs, files) in os.walk(sample): 
        for file in files:
            if 'alignment_quality' in file and file.endswith('txt'):
                path = root + '/' + file
                
                print(path)
                
                dfq = pd.read_csv(path, sep='\t', header=0, index_col=0)
                
                dfq_tables.append(dfq)
                
            if 'TF_targets' in file:
                path = root + '/' + file
                
                print(path)
                
                matcher = re.search(r'p(\d+)', file)
                p = int(matcher.group(1))
                
                try:
                    df = pd.read_csv(path, sep='\t', header=0)
                    peak_counts[p][sample] = df.shape[0]
                except:
                    continue
                
                

dfq = pd.concat(dfq_tables, axis=0)

counts = np.zeros((dfq.shape[0], len(peak_counts)), dtype=int)

pv = list(sorted(peak_counts))
samples = dfq.index.values   

for i in range(0, len(pv)):
    p = pv[i]
    
    for j in range(0, samples.size):
        sample = samples[j]
        
        if sample in peak_counts[p]:
            counts[j][i] = peak_counts[p][sample]

dfp = pd.DataFrame(counts)
dfp.index = samples
dfp.columns = ['p 10^-{}'.format(p) for p in pv]

dft = pd.concat((dfq, dfp), axis=1)

dft.to_csv('alignment_qc_{}_{}.txt'.format(time.strftime('%Y-%m-%d'), int(time.time())), sep='\t', header=True, index=True)
        
