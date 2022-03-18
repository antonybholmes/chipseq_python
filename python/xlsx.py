# -*- coding: utf-8 -*-
"""
Create an excel file

Created on Tue Sep 23 13:34:02 2014

@author: Antony Holmes
"""

import sys

import xlsxwriter

import excel

def create_xlsx(text_file, xlsx_file, locked):
  sys.stderr.write("Creating " + xlsx_file + " from " + text_file + "...\n")
  
  header = excel.get_header(text_file)
  
  left_indices = set()
  right_indices = set()
  center_indices = set()
  
  workbook = xlsxwriter.Workbook(xlsx_file)

  worksheet = excel.create_worksheet(text_file, workbook, left_indices, right_indices, center_indices, locked)  

  used = set()

  excel.set_default_widths(worksheet, header, used)

  workbook.close()
  
  
file = sys.argv[1]
locked = len(sys.argv) > 2 and sys.argv[2] == "locked"

xlsx_file = excel.create_xlsx_file(file, locked)

create_xlsx(file, xlsx_file, locked)
