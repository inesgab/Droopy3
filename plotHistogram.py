import os, glob, sys
# 3 juin 2025
os.chdir('/Users/inesgabert/Documents/LBE/Droopyv2')



module_path = os.path.abspath(os.path.join('./'))
if module_path not in sys.path:
    sys.path.append(module_path)
from functionsCleanPipeline import *

path='/Users/inesgabert/Documents/LBE/experiences/EXP20250513_1116/analysis/'
labels = ['SBW25-WT_CAA_RFP_1','SBW25-WT_CAA_RFP_4', 'SBW25-WT_CAA_RFP_16', 'SBW25-WT_CAA_RFP_64', 'SBW25-WT_CAA_RFP_256', 'SBW25-WT_CAA_RFP_1024']
getDataHistogram(labels,path,'RFP',True)