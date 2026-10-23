import numpy as np
import matplotlib.pyplot as plt
import importlib

import sys
from pathlib import Path

import Diff_System
importlib.reload(Diff_System)

import ToolBox
importlib.reload(ToolBox)

params_1 = {'name': "FIRST", 'd':50, 'K':1, 'mu':1, 'sigma':1, 'T':10, 'delta_t':0.01, 'num_meth':"FEuler"}
params_2 = {'name': "SECOND", 'd':80, 'K':0.5, 'mu':1, 'sigma':1, 'T':20, 'delta_t':0.01, 'num_meth':"FEuler"}

params = params_2