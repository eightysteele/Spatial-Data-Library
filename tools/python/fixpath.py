import os
import sys

DIR_PATH = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
LIB_PATH = os.path.abspath(os.path.realpath('../../lib'))
TMG_PATH = os.path.abspath(os.path.realpath('../../app'))
EXTRA_PATHS = [
  DIR_PATH,
  LIB_PATH,
  TMG_PATH
]

def fix_sys_path():
  sys.path = EXTRA_PATHS + sys.path


