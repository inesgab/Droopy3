#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Jun 5 2025

@author: maxime ardré, inès gabert
This script is used to plot the calibration of fluorescence / bacterial concentration.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy import inf
from numpy import nan
from scipy import odr
from functionsCleanPipeline import loadData, calibrationFluo, findIdx, func

# ---------------- Entries and parameters ------------------

# These entries are to be modified according to the experiment and the channel used
rootPath = "/Users/inesgabert/Documents/LBE/calibration/EXP20250522_1556/"  # absolute path to the calibration data folder
runToConsider = 2  # run number to consider for calibration
channel = "RFP"  # RFP, GFP, PVD
thresholdLogBlank = 3.5  # limit from which the fit starts, e.g., 3.5
tag = "RFP"  #'mSc' 'GFP'

cellCount = None  # bact. concentration (bact/mL), to be calculated or replaced by a known value to skip following entries
if cellCount is None:
    colonyCounts = (
        np.array([19, 16, 15, 20, 16, 15, 13, 20, 16]) * 1e6
    )  # number of colonies in the drops * corresponding factor of dilution
    dropVolume = 5e-3  # in mL

# ----------------- End of entries --------------------------

#                       -------

# ----------------- Code starts here ------------------------

[dropMap, dataFrame, label] = loadData(rootPath + "analysis/")
calibrationData = calibrationFluo(dataFrame, dropMap, label, channel, runToConsider)

calibrationData = calibrationData.replace(
    [-inf], nan
)  # remove the -inf value coming from the log of the raw fluo
selectedColumns = []

for i in range(len(calibrationData.columns)):
    if calibrationData.columns[i].split("_")[2] == tag:
        if "SBW25-PvdS-GFP-CAA-1.3e3" != calibrationData.columns[i]:  # bug sur la dilution??
            selectedColumns.append(calibrationData.columns[i])


channelData = calibrationData[selectedColumns]
meanFluo = channelData.mean()  # fluo area*speed/L already in log
stdFluo = channelData.std()
columnNames = channelData.columns
dilutions = [float(name.split("_")[-1]) for name in columnNames]

# calculating concentration bact/mL
if cellCount is None:
    tubeCounts = colonyCounts / dropVolume
    relStdError = np.std(tubeCounts) / np.mean(tubeCounts)
    cellCount = np.mean(tubeCounts)

# filling blank values
blankIndex = findIdx(dilutions, 0)
dilutions[blankIndex] = 1e4
fitIndices = [i for i, x in enumerate(meanFluo) if x > thresholdLogBlank]
# c = pd.Series([count/10**i  for i in d]) #when in the name there is the exponent factor of dilution
concentration = pd.Series([cellCount / i for i in dilutions])  # when in the name there is the dilution
concentration.index = list(meanFluo.index.values.tolist())

# Regression using Orthogonal Distance Regression (ODR)
quadModel = odr.Model(func)
odrData = odr.Data(
    meanFluo[fitIndices], np.log(concentration[fitIndices]), wd=1.0 / np.power(stdFluo[fitIndices] * 2, 2), we=1.0 / np.power(relStdError, 2)
)
# Set up ODR with the model and data and run
odrInstance = odr.ODR(odrData, quadModel, beta0=[1, 15])
odrOutput = odrInstance.run()

# Output fit parameters
fitParams = odrOutput.beta
fitErrors = odrOutput.sd_beta
print("fit parameter 1-sigma error")
print("———————————–")
for i in range(len(fitParams)):
    print(str(fitParams[i]) + " +- " + str(fitErrors[i]))

# Confidence level curves and other plots
nStd = 2  # to draw 1-sigma intervals
fitParamsUpper = fitParams + nStd * fitErrors
fitParamsLower = fitParams - nStd * fitErrors

a = fitParams[0]  # slope
b = fitParams[1]  # intercept
aErr = fitErrors[0]  # error on slope
bErr = fitErrors[1]  # error on intercept

xValues = np.arange(-5, 8)
fitCurveUpper = func(fitParamsUpper, xValues)
fitCurveLower = func(fitParamsLower, xValues)

# Tick marks
minorYTicks = [x * 10**e for e in range(2, 10) for x in [1, 2, 4, 6, 8]]
mainYTicks = [1 * 10**e for e in range(2, 10)]
minorXTicks = np.arange(-6, 10, 0.5)

# plot the results
fig, ax = plt.subplots()
ax.semilogy([meanFluo[blankIndex], meanFluo[blankIndex]], [1e5, 1e10], "k")
ax.semilogy(
    xValues[xValues > thresholdLogBlank],
    np.exp(func([a, b], xValues[xValues > thresholdLogBlank])),
    "--",
    color="k",
    label="log(C) ~ " + f"\n{a:.2f}({aErr:.2f})*log(F)" + f"\n+{b:.2f}({bErr:.2f})",
)
ax.errorbar(meanFluo, concentration, xerr=2 * stdFluo, yerr=(relStdError * concentration).values, fmt="k.", capsize=3)
ax.fill_between(
    xValues[xValues > thresholdLogBlank],
    np.exp(fitCurveUpper[xValues > thresholdLogBlank]),
    np.exp(fitCurveLower[xValues > thresholdLogBlank]),
    color="grey",
    alpha=1,
    label="+/- sigma interval",
)
ax.set_yscale("log")
ax.set_xlabel("log(raw fluo) (au)")
ax.set_ylabel("bacterial concentration (cell/ml)")
ax.set_yticks(
    minorYTicks,
    minor=True,
)
ax.set_yticks(mainYTicks)

ax.set_title(["CAA SBW25-WT-", channel])
ax.set_xticks(minorXTicks, minor=True)
ax.grid(True, which="major", ls="-", color="k")
ax.grid(True, which="minor", color="0.5", linestyle="-")
ax.legend(loc="lower right")
ax.set_ylim([1, 5e10])
ax.set_xlim([-6, 10])
fig.savefig(
    rootPath + "calibration_PMT" + channel + ".pdf", format="pdf", bbox_inches="tight"
)

# save the calibration data
calibrationTable = pd.DataFrame()
calibrationTable["fluo"] = meanFluo
calibrationTable["Conc"] = concentration
calibrationTable["s"] = stdFluo
calibrationTable.to_csv(rootPath + "calib.csv")