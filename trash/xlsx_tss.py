# -*- coding: utf-8 -*-
"""
Create an excel file for the tss files

Created on Tue Sep 23 13:34:02 2014

@author: Antony Holmes
"""

import sys

import xlsxwriter

import excel

def create_xlsx(text_file, xlsx_file, locked):
  workbook = xlsxwriter.Workbook(xlsx_file)

  worksheet = workbook.add_worksheet()
  
  if locked:
    worksheet.protect(excel.PASSWORD)
  
  header = excel.get_header(text_file)
  
  excel.create_header(header, workbook, worksheet)
  
  count_style = workbook.add_format();
  count_style.set_font_size(11);
  count_style.set_font_name("Arial")
  count_style.set_align("right")
  count_style.set_align("vcenter")
  count_style.set_text_wrap()
  
  f = open(text_file, 'r')
  
  f.readline()
  
  r = 1
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    worksheet.write_number(r, 0, float(tokens[0]), count_style)
    worksheet.write_number(r, 1, float(tokens[1]), count_style)
    
    r += 1
    
  f.close()
  
  used = set()
  
  excel.set_default_widths(worksheet, header, used)

  workbook.close()
  
  
file = sys.argv[1]
locked = len(sys.argv) > 2 and sys.argv[2] == "locked"

xlsx_file = excel.create_xlsx_file(file, locked)

create_xlsx(file, xlsx_file, locked)
