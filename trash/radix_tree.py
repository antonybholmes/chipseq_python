# -*- coding: utf-8 -*-
"""
Data Structures

Created on Tue Aug 12 19:07:11 2014

@author: Antony Holmes
"""

import collections
import sys

class RadixTree(object):
  def __init__(self, word):
    #sys.stdout.write(word + "\n")
    
    char = word[0]
    
    self.char = char 
    self.map = collections.defaultdict(RadixTree)
    
    if len(word) == 1:
      return
      
    tree = RadixTree(word[1:len(word)])
    
    self.map[char] = tree
    
  def search(self, word):
    char = word[0]
    
    ret = char == self.char
    
    if ret == False:
      return ret
    
    if len(word) == 1:
      return True
    
    # the result is the and of all pipes
    return True and self.map[char].search(word[1:len(word)])
    

#tree = RadixTree("cake")

#sys.stderr.write(str(tree.search("cake")) + "\n")