
import matplotlib as plt
import numpy as np
import os, glob, sys
os.chdir('/Users/inesgabert/Documents/LBE/Droopyv2')



module_path = os.path.abspath(os.path.join('./'))
if module_path not in sys.path:
    sys.path.append(module_path)
from functionsCleanPipeline import *

path ='/Users/inesgabert/Documents/LBE/calibration/EXP20250522_1550/analysis'


[dropMap, df, label] = loadData(path)
print(dropMap)
#print(df)
print(label)
channel = 'RFP'
runToConsider = 2
calib = calibrationFluo(df, dropMap, label, channel, runToConsider)
print(calib)