import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(SCRIPT_DIR)
# sys.path.append(os.path.dirname(PARENT_DIR))
sys.path.append(os.path.dirname(SCRIPT_DIR))