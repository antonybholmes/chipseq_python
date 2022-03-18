# -*- coding: utf-8 -*-

import subprocess

def cmd(cmd):
  """
  Run a command line program and return the output as a string.
  
  Args:
    cmd     An array of strings to form an command line expression
            such as ['cat', 'file.txt']
  """
  return subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout.read()

def cmd_lines(cmd):
  """
  Run a command line program and return the output as a list of
  string.
  
  Args:
    cmd     An array of strings to form an command line expression
            such as ['cat', 'file.txt']
  """
  return subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout.readlines()
