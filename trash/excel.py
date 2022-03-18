# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 10:19:32 2014

@author: Antony Holmes
"""

import re
import sys

import xlsxwriter

LONG_COLUMN_WIDTH = 160
MEDIUM_COLUMN_WIDTH = 32
DEFAULT_COLUMN_WIDTH = 16
SHORT_COLUMN_WIDTH = 10

PASSWORD = "tmwrnj"

def create_xlsx_file(file, locked):
  if locked:
    return re.sub(r'\.[^\.]+$', '_L.xlsx', file)
  else:
    return re.sub(r'\.[^\.]+$', '.xlsx', file)
  

def create_workbook(file, locked):
  """
  Create a workbook with a standardized file name
  
  file      the name of the excel file to be created.
  locked    True if the file should be password locked, otherwise False.  
  """
  
  xlsx_file = create_xlsx_file(file, locked)
  
  workbook = xlsxwriter.Workbook(xlsx_file)
  
  return workbook
  
  
def find_first_column(text, header):
  """
  Given a header, find the first column matching some text.

  text      Text to search for.
  header    An array of header names  
  """
  
  for i in range(0, len(header)):
    if text in header[i]:
      return i
      
  return -1


def set_column_width(worksheet, column, width):
  if column > -1:
    worksheet.set_column(column, column, width)
  
  return column
  
def get_header(text_file):
  """
  Returns the header of a text file as an array of names.
  """
  
  f = open(text_file, 'r')
  
  header = f.readline().strip().split("\t")

  f.close()  
  
  return header
  

def create_header(header, workbook, worksheet):
  """
  Adds a formatted, standardized header to a worksheet.
  """
  
  # Use a centered, bold font for the header
  header_font = workbook.add_format();
  header_font.set_font_size(11);
  header_font.set_font_name("Arial")
  header_font.set_bold()
  header_font.set_align("center")
  header_font.set_align("vcenter")
  header_font.set_text_wrap()
  
  for i in range(0, len(header)):
    worksheet.write_string(0, i, header[i], header_font)
  
  
def create_worksheet(text_file, workbook, left_indices, right_indices, center_indices, locked):
  """
  Creates a new worksheet in a workbook from a tabular text file.
  
  text_file         A tab delimited text file to turn into an Excel file.
  workbook          The workbook to create the new sheet in.
  left_indices      indices to left align.
  """
  
  worksheet = workbook.add_worksheet()

  left_font = workbook.add_format()
  left_font.set_font_size(11);
  left_font.set_font_name("Arial")
  left_font.set_align("left")
  left_font.set_align("vcenter")
  left_font.set_text_wrap()
  
  right_font = workbook.add_format();
  right_font.set_font_size(11);
  right_font.set_font_name("Arial")
  right_font.set_align("right")
  right_font.set_align("vcenter")
  right_font.set_text_wrap()

  center_font = workbook.add_format();
  center_font.set_font_size(11);
  center_font.set_font_name("Arial")
  center_font.set_align("center")
  center_font.set_align("vcenter") 
  center_font.set_text_wrap()

  f = open(text_file, 'r')
  
  header = f.readline().strip().split("\t")
  
  columns = len(header)
  
  create_header(header, workbook, worksheet)
 
  r = 1
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")

    #sys.stderr.write("is_num " + line + "\n")    
    
    for i in range(0, columns):
      is_number = re.match(r'^-?\d+((\.\d+)([Ee][\+\-]?\d+)?)?$', tokens[i])
      is_entrez = re.match(r'.*[Ee]ntrez.*', header[i])
      
      if i in left_indices:
        font = left_font
      elif i in right_indices:
        font = right_font
      elif i in center_indices:
        font = center_font
      elif is_entrez:
        font = center_font
      elif is_number:
        # default numbers on the right if no other type of preference set
        font = right_font
      else:
        font = center_font
      
      if is_number:
        # Two decimal places for numbers
        worksheet.write_number(r, i, round(float(tokens[i]), 2), font)
      else:
        worksheet.write_string(r, i, tokens[i], font)
    
    r += 1
    
  f.close()
  
  if locked:
    worksheet.protect(PASSWORD)
  
  return worksheet
  
 
def set_default_widths(worksheet, header, used):
  for i in range(0, len(header)):
    if i in used:
      continue
    
    if "Probe ID" in header[i]:
      set_column_width(worksheet, i, MEDIUM_COLUMN_WIDTH)
    else:
      set_column_width(worksheet, i, DEFAULT_COLUMN_WIDTH)