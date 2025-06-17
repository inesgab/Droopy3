#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 21:31:36 2019

@author: maxime ardré
"""

import pandas as pd
import numpy as np
from scipy.signal import chirp, find_peaks, peak_widths
import os
import shutil
import ast
from os import listdir
from os.path import isfile, join, isdir
import matplotlib.pyplot as plt
from PIL import Image
import traceback
import multiprocessing
from inspect import currentframe, getframeinfo


params = {
    "legend.fontsize": "x-large",
    "figure.figsize": (5, 5),
    "axes.labelsize": "x-large",
    "axes.titlesize": "x-large",
    "xtick.labelsize": "x-large",
    "ytick.labelsize": "x-large",
}
plt.rcParams.update(params)

import scipy.optimize
import scipy.interpolate as si
from scipy.optimize import curve_fit
from scipy.stats import norm
import seaborn as sn

# from ggplot import *
from pandas import Series as pds

import bisect
import csv
import glob

import sys

module_path = os.path.abspath(
    os.path.join("/Users/inesgabert/Documents/LBE/Droopyv2/fitderivpackage1.02")
)

if module_path not in sys.path:
    sys.path.append(module_path)
from fitderiv import fitderiv


def fluoToBact_loglog(y):
    return 0.94 * y + 8.86 * np.log(10)


def getValueInLabel(label, path):
    p = path.split("/")
    posOfFolderInSource = -4
    # try:
    if "Exp_Incub_4h30-6h30-24h-52h" == p[posOfFolderInSource]:
        # c = label.split('SBW25GFP-')
        # d = c[1]
        # b = d.split('h N=')
        c = label.split("h")
        b = c[0]

    if "RFP_inoculum" == p[posOfFolderInSource]:
        c = label.split("SBW25-WT_CAA_RFP_")
        b = c[1]

    if "Exp_Dilution_PvdS" == p[posOfFolderInSource]:
        c = label.split("CAA-PvdS-")
        d = c[1]
        e = d.split(" N")
        b = e[0]

    if "Exp_DilutionWTstock2" == p[posOfFolderInSource]:
        c = label.split("CAA-WT-")
        e = c[1]
        d = e.split(" N")
        b = d[0]

    if "Exp_Dilution_5_50_500_3500" == p[posOfFolderInSource]:
        c = label.split("CAA-WT-")
        d = c[1]
        e = d.split(" N")
        b = e[0]

    if (
        "CAA" == p[posOfFolderInSource]
        or "CAAnoGly" == p[posOfFolderInSource]
        or "CAA1" == p[posOfFolderInSource]
        or "CAAspp" == p[posOfFolderInSource]
        or "AyshaCAA1" == p[posOfFolderInSource]
    ):
        c = label.split("SBW25-WT_CAA_")
        b = c[1]

    if "M9glyFromCAA" == p[posOfFolderInSource]:
        c = label.split("SBW25-WT-M9gly-")
        b = c[1]

    if "SMM" == p[posOfFolderInSource]:
        c = label.split("SBW25-WT_SMM_")
        b = c[1]

    if "M9gly" == p[posOfFolderInSource]:
        c = label.split("SBW25-WT-M9glyStock-")
        b = c[1]

    if "Exp_IronCAA_WT" == p[posOfFolderInSource]:
        c = label.split("CAA-WT-")
        d = c[1]
        e = d.split("uM")
        b = e[0]

    if "Exp_BipyCAA_WT" == p[posOfFolderInSource]:
        c = label.split("uM")
        b = c[0]

    if "Exp_BipyCAA_PvdS" == p[posOfFolderInSource]:
        c = label.split("uM")
        b = c[0]

    if "Exp_Incubation_14h-19h-43h" == p[posOfFolderInSource]:
        c = label.split("SBW25GFP-")
        d = c[1]
        e = d.split("h")
        b = e[0]

    if "Exp_Incub_4h30-6h30-24h-52h" == p[posOfFolderInSource]:
        c = label.split("SBW25-WT_CAA_")
        d = c[1]
        e = d.split("h")
        b = e[0]

    if "Exp_Incubation_4h30_6h30_24h_43h13_1bactPerDrop" == p[posOfFolderInSource]:
        c = label.split("SBW25-WT_CAA_")
        d = c[1]
        if d[:3] == "pur":
            b = 0
        else:
            e = d.split("h")
            b = e[0]

    if "Exp_Incubation_merge" == p[posOfFolderInSource]:
        f = p[posOfFolderInSource + 1]
        ff = f.split("dilution_CAA_WT-")

        if ff[-1] == "stock24h30-6h30-24h-52h":
            c = label.split("SBW25-WT_CAA_")
            d = c[1]
            e = d.split("h")
            b = e[0]

        if ff[-1] == "stock14h33-19h-43h":
            c = label.split("SBW25GFP-")
            d = c[1]
            e = d.split("h")
            b = e[0]

    if "Exp_IronPvd" == p[posOfFolderInSource]:
        c = label.split("CAA-pur-")
        d = c[1]
        e = d.split("uM")
        b = e[0]

    if "Exp_IronCAA_PvdS" == p[posOfFolderInSource]:
        c = label.split("CAA-PvdS-")
        d = c[1]
        e = d.split("uM")
        b = e[0]

    if "Exp_Dilution_M9glycerol" == p[posOfFolderInSource]:
        c = label.split("SBW25-WT_M9gly_")
        b = c[1]

    if "CAAPvdS" == p[posOfFolderInSource]:
        c = label.split("SBW25-PvdS_CAA_")
        b = c[1]

    if "Exp_GreenRed" == p[posOfFolderInSource]:
        val = label

    if (
        "Exp_GreenRed" == p[-6]
    ):  # experiment green red to measure the yield as a function of ratio
        c = label.split("g")
        b = c[1].split("-")
        val = float(b[0])
    if "antib" in p[posOfFolderInSource]:
        c = label.split("_")
        b = float(c[-2]), float(c[-1])
        val = b

    else:   
        val = float(b)

    return val


def findIdx(your_list, item):
    lst = pd.DataFrame(your_list)
    lst2 = lst[0].values
    lst3 = lst2.tolist()
    return lst3.index(item)


def getOrderLabel(orderList, toOrganize):
    order = []
    for idx in toOrganize:
        order.append(orderList.index(idx))
    return order


def organizeLabel(order, toOrganize):
    zipped = zip(toOrganize, order)
    organized = sorted(zipped, key=lambda x: x[1])
    l = zip(*organized)
    return list(l)[0]


def loadData(path):
    if path[-1] != "/":
        path = path + "/"

    drpfiles = sorted(
        [
            path + "droplets/" + f
            for f in listdir(path + "droplets/")
            if isfile(join(path + "droplets/", f)) and f != ".DS_Store"
        ]
    )
    i = 0
    df = {}
    for file in drpfiles:
        df[i] = pd.read_csv(file)
        i += 1

    tpfile = pd.read_csv(path + "droplet.csv")
    dropMap = []
    label = list(set(tpfile["group"]))
    label2 = [l for l in label if l != "Empty" and l != "CAA"]
    dropMap = np.array(list(zip(tpfile["well"], tpfile["group"])))

    return [dropMap, df, label2]


def getfitData(
    df, j, channel, startTime, timeThresh, threshOutlayer, incCarte, noGrowthThresh
):
    """
    get the time series of a droplet for a specific channel. It returns log(fluo_area*speed_of_drop/size_of_drop) +/- cov
    The parameters are :
    df : the dataframe containing the droplet data
    j : the droplet number that you want to study
    channel : GFP, RFP, PVD, PVDoverGFP
    startTime : the time at which we want the start of the time serie. Allows to cut the weird points at the beginning of the measurment
    timeThresh : the time at which we want the end of the time serie
    threshOutlayer : the time at which if the drop grows we consider it is a contamination. The drop is flagged outliers
    incCarte : the uncertainty of the electronic card
    noGrowthThresh: the voltage under which we consider that the drop did not grow properly. The drop is flagged outliers
    """
    y = []
    x = []
    seuilDetection = incCarte
    # if incCarte<0:
    #    seuilDetection=-10
    # else:
    #    seuilDetection=-10#2e-2V checked on an empty droplet
    # seuilDetectionPVD=3e-2V checked on an empty droplet codé en dure

    try:
        # default stuff to do get the fluorescence in the drop decouple from speed and size
        [x, y, idxThresDetect] = dataUpAndDown(
            df, j, channel, seuilDetection
        )  # return the time serie with two other series corresponding to the data plus or minus the cov

        # special case
        if channel == "SpeedSize":
            p = df[j]
            idxThresDetect = 1  # find threshold of data above detection 1e-2V
            y = np.array(p.speed / p.size)
            x = np.array(p.time / 3600)
            # x=x-x[0] #reset time
            ystd = []

            # remove nan value
            idxNan = [i for i, j in enumerate(np.isnan(y)) if j]
            x = np.delete(x, idxNan)
            y = np.delete(y, idxNan)
            ystd = np.delete(ystd, idxNan)

        if channel == "PVDoverGFP":
            p = df[j]

            idxThresDetectGFP = bisect.bisect(
                p.fluo_3_median, seuilDetection
            )  # find threshold of data above detection 5e-3V
            idxThresDetectPVD = bisect.bisect(p.fluo_1_median, seuilDetection)

            idxThresDetect = max(idxThresDetectPVD, idxThresDetectGFP)

            yFluo = np.array(
                np.log(p["fluo_3_area"] * p["speed"] / p["size"]), dtype=np.float64
            )  # to screen for nonoutlayer
            y = np.array(p["fluo_1_area"] / p["fluo_3_area"])  # no need of log here
            x = np.array(p.time / 3600, dtype=np.float64)
            # x=x-x[0] #reset time
            ystd = np.array(
                np.divide(p["fluo_1_std"] + incCarte / 2, p["fluo_1_median"]),
                dtype=np.float64,
            )

            # remove nan value
            idxNan = [i for i, j in enumerate(np.isnan(y)) if j]
            x = np.delete(x, idxNan)
            y = np.delete(y, idxNan)
            ystd = np.delete(ystd, idxNan)

        # filtre nonoutlayer
        # the nonoutlayer are detected with the median of the signal
        # the threshold of detection is detected with the median of the signal. The minimum detection is 2e-2V
        # signal above empty droplets. The plotted signal is an integrated value of the fluo (it is not comparable to 2e-2V directly)

        # idxStart=max(idxThresDetect,bisect.bisect(x, startTime)) #get the first index from which the signal is above the detection 'seuilDetectionGFP' or after the time threshold applied by the operator
        idxStart = bisect.bisect(x, startTime)
        idxThres = bisect.bisect(x, min(max(x), timeThresh))  # get the index
        idxOutlayer = bisect.bisect(x, min(max(x), threshOutlayer))
        if len(y[:, 0]) > 2:
            if len(y[:, 0]) > idxOutlayer:
                if channel == "GFP" or channel == "PVD" or channel == "RFP":
                    nonoutlayer = y[idxOutlayer, 0] > y[2, 0] + 1 or threshOutlayer <= 0
                elif channel == "PVDoverGFP":
                    nonoutlayer = (
                        yFluo[idxOutlayer] > yFluo[2] + 1 or threshOutlayer <= 0
                    )
                elif channel == "SpeedSize":
                    nonoutlayer = True
                else:
                    nonoutlayer = False  # if it is an outlayer it must be discard and you set it to false
                    print("max value not high enough")

            else:
                print("BUG idxOutlayer not long enough")
                nonoutlayer = False

            if max(y[:, 0]) < noGrowthThresh:
                nonoutlayer = False
                print("did not get above noGrowthThresh= " + str(noGrowthThresh))
        else:
            print("timeserie to short drop= " + str(j) + " y[:,0]=")
            print(y[:, 0])
            nonoutlayer = False

    except Exception as inst:
        nonoutlayer = False
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)
        print("bug checkdata")
        print(type(inst))  # the exception instance
        print(inst.args)  # arguments stored in .args
        _, _, tb = sys.exc_info()
        traceback.print_tb(tb)  # Fixed format
        tb_info = traceback.extract_tb(tb)
        filename, line, func, text = tb_info[-1]
        print("An error occurred on line {} in statement {}".format(line, text))
        print(y)
        print("drop " + str(j))

    if nonoutlayer:
        return [x[idxStart:idxThres], y[idxStart:idxThres, :], True]
    else:
        print("outlayer")
        return [[], [], False]


def dataUpAndDown(df, j, channel, seuilDetection):
    Nb = ""
    if channel == "RFP":
        Nb = "2"
    elif channel == "GFP":
        Nb = "3"
    elif channel == "PVD":
        Nb = "1"
    else:
        print("Not a simple fluo channel")

    if Nb != "":
        p = df[j]
        # idxThresDetect=bisect.bisect(p['fluo_'+Nb+'_median'],seuilDetection)#find threshold of data above detection 5e-3V

        # y=np.array( np.log(p['fluo_2_area']*p['speed']/p['size']) - np.log(p['fluo_2_area'][1]*p['speed'][1]/p['size'][1]) )

        x = np.array(p.time / 3600)
        y = np.array(p["fluo_" + Nb + "_area"] * p["speed"] / p["size"])

        # remove nan value
        threshArea = seuilDetection  # threshold of detection for the area
        idxNan = [i for i, j in enumerate(np.isnan(y)) if j]
        x = np.delete(x, idxNan)
        y = np.delete(y, idxNan)
        fluo_cov = np.delete(p["fluo_" + Nb + "_cov"], idxNan)
        yUp = np.array(
            y * (1 + fluo_cov)
        )  # the data cov are the std/mean for the channel
        yDown = np.array(y * (1 - fluo_cov))
        yDown[yDown < 0] = threshArea

        # test inès
        og_order = np.arange(len(y))
        idx_sorted = np.argsort(y)  # tri selon x
        x = x[idx_sorted]
        y = y[idx_sorted]
        sorted_order = og_order[idx_sorted]  # ordre trié
        # fin test inès

        idxThresDetect = bisect.bisect(
            y, seuilDetection
        )  # find threshold of data above detection 5e-3V


        # to have a nice fit, set to the value that are below the detection threshold to the detection threshold


        if idxThresDetect > 0:
            print("idxThresDetect: " + str(idxThresDetect))
            print("yUp" + str(yUp))
            if idxThresDetect >= len(y):
                print("idxThresDetect is too high, setting to last value")
                idxThresDetect = len(y) - 1
            indices = yUp <= y[idxThresDetect]
            yUp[indices] = threshArea * 2  # replacing 0s with seuilDetection
            indices = yDown <= y[idxThresDetect]
            yDown[indices] = threshArea / 2  # replacing 0s with seuilDetection
            indices = y <= y[idxThresDetect]
            y[indices] = threshArea  # replacing 0s with seuilDetection

        try:
            print("calculate the values")
            print("drop " + str(j))
            logy = np.log(y)
            logyUp = np.log(yUp)
            logyDown = np.log(yDown)

        except Exception:
            _, _, tb = sys.exc_info()
            traceback.print_tb(tb)  # Fixed format
            tb_info = traceback.extract_tb(tb)
            filename, line, func, text = tb_info[-1]
            print("An error occurred on line {} in statement {}".format(line, text))
            print('y:', y)
            print('yUp:', yUp)
            print('yDown:', yDown)
            return [x, np.column_stack(([], [], [])), idxThresDetect]

        # restoring original order (before bisect) to keep the right data for each label
        x_og = np.empty_like(x)
        x_og[sorted_order] = x
        logy_og = np.empty_like(logy)
        logy_og[sorted_order] = logy
        logyUp_og = np.empty_like(logyUp)
        logyUp_og[sorted_order] = logyUp
        logyDown_og = np.empty_like(logyDown)
        logyDown_og[sorted_order] = logyDown

        return [x_og, np.column_stack((logy_og, logyUp_og, logyDown_og)), idxThresDetect]
    else:
        print("drop does not exist " + str(j))
        return [x, np.column_stack(([], [], [])), idxThresDetect]


def calibrationFluo(df, dropMap, label, channel, runToConsider):
    calib = pd.DataFrame()
    for k, dataFile in enumerate(label):
        fluo = []

        for j, content in enumerate(dropMap[:, 1]):
            y = []
            x = []

            if content == dataFile:
                [x, y, flag] = getfitData(df, j, channel, 0, 20, -1, -15, -2)
                
                print("y:", y)
                yy = y[:, 0]
                #print(yy[runToConsider])
                fluo.append(yy[runToConsider])

        #print(dataFile)
        #print(fluo)
        calib[dataFile] = pd.Series(fluo)

    return calib


def calib(f, calibParam):
    c = np.exp(f * calibParam[0] + calibParam[1])
    return c


def apply_calib_to_df(df, calibParam):
    threshN = 8

    for key, data in df.items():
        if "fluo_3_area" in data:
            # Apply this function to row fluo_area_3 and overwrite values -> now c = cells/mL (get exponential function from area)
            data["fluo_3_calib"] = calib(data["fluo_3_area"], calibParam)
            # save
            data.to_csv(f"path_to_save/{key}.csv", index=False)

            if threshN in data["fluo_3_area"].values:
                y = data["fluo_3_area"](threshN)

            if threshN not in data["fluo_3_area"].values:
                # interpolation
                y = data["fluo_3_area"].interpolate(method="linear")(threshN)


def checkfitData(df, dropMap, label, channel, folder, parameters):
    parameters = parameters.loc[channel]

    startTime = parameters["startTime"]
    timeThresh = parameters["timeThresh"]
    threshOutlayer = parameters["threshOutlayer"]
    incCarte = parameters["incCarte"]
    ymin = parameters["ymin"]
    ymax = parameters["ymax"]
    noGrowthThresh = parameters["noGrowthThresh"]

    for k, dataFile in enumerate(label):
        fig, ax = plt.subplots()
        date = folder.split("/")
        plt.title(date[-3] + " " + dataFile)
        # plt.ylim(np.log([1e-4,5]))
        ax.set_xlim([0, timeThresh])

        err = np.array([])
        for j, content in enumerate(dropMap[:, 1]):
            y = []
            x = []

            if content == dataFile:
                [x, y, flag] = getfitData(
                    df,
                    j,
                    channel,
                    startTime,
                    timeThresh,
                    threshOutlayer,
                    incCarte,
                    noGrowthThresh,
                )

                if flag == True:
                    if channel == "RFP":
                        c = "r"
                    if channel == "GFP":
                        c = "b"
                        if not ymin and not ymax:
                            ax.set_ylim([np.min(y[:, 0]) / 2, 2 * np.max(y[:, 0])])
                        else:
                            ax.set_ylim([ymin, ymax])

                    if channel == "PVD":
                        c = "g"
                    if channel == "PVDoverGFP":
                        c = "k"
                        plt.ylim([0, 3])
                    if channel == "SpeedSize":
                        c = "m"
                        a = 100 * np.std(y) / np.mean(y)
                        err = np.concatenate([err, [a]])
                        ax.title(
                            date[-3]
                            + " "
                            + dataFile
                            + " std/m="
                            + str("{0:.2g}".format(np.mean(err)))
                            + "%"
                        )
                        ax.set_ylim([np.min(y) / 2, 2 * np.max(y)])
                        ax.grid(which="both")

                    ax.scatter(x, y[:, 0], color=c)
                    # ax.errorbar(x, y, ystd,  fmt='.r')
                else:
                    print("data removed by thresh drop=" + str(j))

        # ax.grid(b=True, which='major', color='b', linestyle='-')
        # ax.grid(b=True, which='minor', color='r', linestyle='--')
        # plt.draw()
        print(np.mean(np.mean(err)))
        print('folder:', folder)
        fig.savefig(folder + '/'+ dataFile + "checkDataFit" + channel + ".png")


def checkfitDataCalib(df, dropMap, label, channel, folder, parameters):
    parameters = parameters.loc[channel]

    startTime = parameters["startTime"]
    timeThresh = parameters["timeThresh"]
    threshOutlayer = parameters["threshOutlayer"]
    incCarte = parameters["incCarte"]
    ymin = parameters["ymin"]
    ymax = parameters["ymax"]
    noGrowthThresh = parameters["noGrowthThresh"]
    aa = parameters["calib"]
    print("calib " + channel)
    print(aa)
    calibParam = [float(aa[0]), float(aa[1])]

    for k, dataFile in enumerate(label):
        fig, ax = plt.subplots()
        date = folder.split("/")
        plt.title(date[-3] + " " + dataFile)
        # plt.ylim(np.log([1e-4,5]))
        ax.set_xlim([0, timeThresh])

        err = np.array([])
        for j, content in enumerate(dropMap[:, 1]):
            y = []
            x = []

            if content == dataFile:
                [x, y, flag] = getfitData(
                    df,
                    j,
                    channel,
                    startTime,
                    timeThresh,
                    threshOutlayer,
                    incCarte,
                    noGrowthThresh,
                )

                if flag == True:
                    y = (y * calibParam[0] + calibParam[1]) / np.log(
                        10
                    )  # get the log of concentration in log10

                    if not ymin and not ymax:
                        ax.set_ylim([np.min(y[:, 0]) * 0.95, 1.05 * np.max(y[:, 0])])
                    else:
                        ax.set_ylim([ymin, ymax])
                        ax.set_yticks(np.arange(ymin, ymax))

                    if channel == "RFP":
                        c = "r"

                    if channel == "GFP":
                        c = "b"

                    if channel == "PVD":
                        c = "g"

                    if channel == "PVDoverGFP":
                        c = "k"
                        plt.ylim([0, 3])
                    if channel == "SpeedSize":
                        c = "m"
                        a = 100 * np.std(y) / np.mean(y)
                        err = np.concatenate([err, [a]])
                        ax.title(
                            date[-3]
                            + " "
                            + dataFile
                            + " std/m="
                            + str("{0:.2g}".format(np.mean(err)))
                            + "%"
                        )
                        ax.set_ylim([np.min(y) / 2, 2 * np.max(y)])
                        ax.grid(which="both")

                    ax.plot(x, y[:, 0], color=c, linewidth=0.5)
                    # ax.errorbar(x, y, ystd,  fmt='.r')
                else:
                    print("data removed by thresh drop=" + str(j))

        ax.grid(which="major", color="k", linestyle="-")
        plt.xlabel("Time (h)")
        plt.ylabel("log10(cells/mL)")
        # ax.grid(b=True, which='minor', color='r', linestyle='--')
        # plt.draw()
        print(np.mean(np.mean(err)))
        fig.savefig(folder + dataFile + "checkDataFit" + channel + ".png")


# def fitDataIndiv(df, dropMap, label, channel, folder, parameters):


#     parameters = parameters.loc[channel]

#     startTime =  parameters['startTime']
#     timeThresh =  parameters['timeThresh']
#     threshOutlayer = parameters['threshOutlayer']
#     incCarte = parameters['incCarte']
#     display = parameters['display']
#     deletionData = parameters['deletionData']
#     ymin = parameters['ymin']
#     ymax = parameters['ymax']
#     noGrowthThresh= parameters['noGrowthThresh']
#     fromLag = parameters['fromLag']
#     calibParam = parameters['calib']
#     if not calibParam:
#         calibParam = [1, 0]

#     #print('RESET TIME!')
#     pathResults = folder+'resultIndiv/'
#     if not os.path.exists(pathResults):
#         os.makedirs(pathResults)
#     else :
#         if deletionData==True:
#             shutil.rmtree(pathResults)
#             os.makedirs(pathResults)

#     if os.path.exists(folder + 'outliersbyAlgo_'+channel+'.csv') and not fromLag:
#         os.remove(folder + 'outliersbyAlgo_'+channel+'.csv')
#     elif not os.path.exists(folder + 'outliersbyAlgo_'+channel+'.csv') and not fromLag:
#         with open(folder + 'outliersbyAlgo_'+channel+'.csv', 'a') as f:
#             writer = csv.writer(f)
#             writer.writerow([-1])


#     for j,dataFile in enumerate(label):


#         pathFitPlot=folder+'resultIndiv/'+dataFile+'/'
#         combined_csv = pd.DataFrame()

#         if not os.path.exists(pathFitPlot):
#             os.makedirs(pathFitPlot)
#         else :
#             if deletionData==True:
#                 shutil.rmtree(pathFitPlot)
#                 os.makedirs(pathFitPlot)

#         if fromLag:
#             try :
#                 print('fit the data from the lag time found in the previous treatment')
#                 filePreTreatment = pd.read_csv(folder+'resultIndiv/'+dataFile+ 'resultIndiv_Drop_'+channel+'.csv')
#                 filePreTreatment = filePreTreatment.set_index('drop')
#                 fileOutliersAlgo = pd.read_csv(folder+'outliersbyAlgo_'+channel+'.csv', header=None)
#             except:
#                 print(folder+'resultIndiv/'+dataFile+ 'resultIndiv_Drop_'+channel+'.csv')
#                 _, _, tb = sys.exc_info()
#                 traceback.print_tb(tb) # Fixed format
#                 tb_info = traceback.extract_tb(tb)
#                 filename, line, func, text = tb_info[-1]
#                 print('An error occurred on line {} in statement {}'.format(line, text))

#         for j,content in enumerate(dropMap[:,1]):
#             y=[]
#             x=[]

#             if content==dataFile :


#                 if fromLag:
#                     if j in fileOutliersAlgo.values:
#                         print('flagged drop by algo do not fit ' +str(j) )
#                         startTime =  parameters['startTime']
#                         flag = False
#                     else:
#                         startTime = filePreTreatment['lag time'][j]-1
#                         [x,y,flag] = getfitData(df, j, channel,startTime,timeThresh,threshOutlayer, incCarte,noGrowthThresh)
#                         [xPVD,yPVD,flagPVD] = getfitData(df, j, 'PVD',startTime,timeThresh,threshOutlayer, incCarte,noGrowthThresh)

#                 else :
#                     [x,y,flag] = getfitData(df, j, channel,startTime,timeThresh,threshOutlayer, incCarte,noGrowthThresh)
#                     [xPVD,yPVD,flagPVD] = getfitData(df, j, 'PVD',startTime,timeThresh,threshOutlayer, incCarte,noGrowthThresh)


#                 if flag==True:
#                     idx=np.where(np.isinf(y))
#                     y=np.delete(y,idx[0],0)
#                     x=np.delete(x,idx[0],0)

#                     try:
#                         # redefine bounds and run inference
#                         #hyperparameters
#                         #0: amplitude
#                         #1: flexibility
#                         #2: measurement error
#                         #default sets the boundaries for the first hyperparameter to be 10^-1 and 10^4 and the boundaries for the third hyperparameter to be 10^2 and 10^6
#                         #b= {0: [-1, 4], 1: [-1,2], 2: [-2, 6]}
#                         #initial
#                         #b= {0: [-1,8], 1: [-5,50], 2: [-3,2]}
#                         b= {0: [-1, 8], 1: [-2,2], 2: [-3,2]}#paper


#                         #apply on the range from start to max time

#                         q= fitderiv(x, y, bd= b, logs= False, showstaterrors=False)

#                         if display==True:
#                             # plot results
#                             plt.ioff()
#                             plotDerivCurve(j, [np.column_stack((x,x,x)), xPVD ], [y, yPVD[:,0]], [2*(y[:,1]-y[:,0]), 2*(yPVD[:,1]-yPVD[:,0])] , q, startTime, timeThresh, pathFitPlot, dataFile, channel, ymin,ymax, calibParam, incCarte, content, fromLag, display)


#                         # export results
#                         statdict=q.printstats()
#                         # with open(folder+'resultIndiv/'+dataFile+'resultIndiv_Drop'+str(j)+channel+'.csv','w') as csv_file:
#                         #     writer = csv.writer(csv_file)
#                         #     for key, value in statdict.items():
#                         #         writer.writerow([key, value])

#                         statdict['drop'] = j

#                         combined_csv = pd.concat( [combined_csv, pd.DataFrame(statdict,index=[0])])
#                         #print(combined_csv)

#                     except Exception as inst:
#                         print(folder)
#                         print('bug od for drop: '+str(j))
#                         print(type(inst))     # the exception instance
#                         print(inst.args)     # arguments stored in .args
#                         _, _, tb = sys.exc_info()
#                         traceback.print_tb(tb) # Fixed format
#                         tb_info = traceback.extract_tb(tb)
#                         filename, line, func, text = tb_info[-1]
#                         print('An error occurred on line {} in statement {}'.format(line, text))

#                 else:
#                     print("outlayer not fitted drop="+str(j)+" flag="+str(flag)+" label= "+content)
#                     combined_csv = pd.concat( [combined_csv, pd.DataFrame(
#                     {'max df': np.nan,
#                      'max df std': np.nan,
#                      'max df stderr': np.nan,
#                      'time of max df': np.nan,
#                      'time of max df std': np.nan,
#                      'time of max df stderr': np.nan,
#                      'inverse max df': np.nan,
#                      'inverse max df std': np.nan,
#                      'inverse max df stderr': np.nan,
#                      'max y': np.nan,
#                      'max y std': np.nan,
#                      'max y stderr': np.nan,
#                      'lag time': np.nan,
#                      'lag time std': np.nan,
#                      'lag time stderr': np.nan,
#                      'drop' : j},index=[0])])
#                     if not fromLag:
#                         with open(folder + 'outliersbyAlgo_'+channel+'.csv', 'a') as f:
#                             writer = csv.writer(f)
#                             writer.writerow([str(j)])

#         if fromLag:
#             combined_csv.to_csv(folder+'resultIndiv/'+dataFile+'resultIndivFitFromDetectionTime_Drop_'+channel+'.csv')
#         else:
#             combined_csv.to_csv(folder+'resultIndiv/'+dataFile+'resultIndiv_Drop_'+channel+'.csv')


def fitDataIndivParallel(df, dropMap, label, channel, folder, parameters):
    if (
        os.path.exists(folder + "outliersbyAlgo_" + channel + ".csv")
        and not parameters["fromLag"][channel]
    ):
        os.remove(folder + "outliersbyAlgo_" + channel + ".csv")
    elif (
        not os.path.exists(folder + "outliersbyAlgo_" + channel + ".csv")
        and not parameters["fromLag"][channel]
    ):
        with open(folder + "outliersbyAlgo_" + channel + ".csv", "a") as f:
            writer = csv.writer(f)
            writer.writerow([-1])

    num_cores = multiprocessing.cpu_count()
    # processed_list = Parallel(n_jobs=num_cores)(delayed(wrapFitData)(df, folder,dataFile,dropMap,parameters,channel) for j,dataFile in enumerate(label))
    # print(processed_list)
    for j, dataFile in enumerate(label):
        wrapFitData(df, folder, dataFile, dropMap, parameters, channel)


# function used in gui
def wrapFitData(df, folder, dataFile, dropMap, parameters, channel):
    parameters = parameters.loc[channel]

    deletionData = parameters["deletionData"]
    startTime = parameters["startTime"]
    timeThresh = parameters["timeThresh"]
    threshOutlayer = parameters["threshOutlayer"]
    incCarte = parameters["incCarte"]
    display = parameters["display"]
    fromLag = parameters["fromLag"]
    ymin = parameters["ymin"]
    ymax = parameters["ymax"]
    noGrowthThresh = parameters["noGrowthThresh"]
    calibParam = parameters["calib"]
    if not calibParam:
        calibParam = [1, 0]

    frameinfo = getframeinfo(currentframe())
    print(frameinfo.filename, frameinfo.lineno)
    print(calibParam)

    pathFitPlot = folder + "resultIndiv/" + dataFile + "/"
    combined_csv = pd.DataFrame()

    pathResults = folder + "resultIndiv/"
    if not os.path.exists(pathFitPlot):
        os.makedirs(pathFitPlot)
    else:
        if deletionData == True:
            shutil.rmtree(pathFitPlot)
            os.makedirs(pathFitPlot)

    if fromLag:
        try:
            print("fit the data from the lag time found in the previous treatment")
            filePreTreatment = pd.read_csv(
                folder
                + "resultIndiv/"
                + dataFile
                + "resultIndiv_Drop_"
                + channel
                + ".csv"
            )
            filePreTreatment = filePreTreatment.set_index("drop")
            fileOutliersAlgo = pd.DataFrame([])
            if os.path.exists(folder + "outliersbyAlgo_" + channel + ".csv"):
                fileOutliersAlgo = pd.read_csv(
                    folder + "outliersbyAlgo_" + channel + ".csv", header=None
                )

        except:
            print(
                folder
                + "resultIndiv/"
                + dataFile
                + "resultIndiv_Drop_"
                + channel
                + ".csv"
            )
            _, _, tb = sys.exc_info()
            traceback.print_tb(tb)  # Fixed format
            tb_info = traceback.extract_tb(tb)
            filename, line, func, text = tb_info[-1]
            print("An error occurred on line {} in statement {}".format(line, text))

    for j, content in enumerate(dropMap[:, 1]):
        y = []
        x = []

        if content == dataFile:
            if fromLag:
                if j in fileOutliersAlgo.values:
                    print("flagged drop by algo do not fit " + str(j))
                    startTime = parameters["startTime"]
                    flag = False
                else:
                    startTime = filePreTreatment["lag time"][j] - 2
                    [x, y, flag] = getfitData(
                        df,
                        j,
                        channel,
                        startTime,
                        timeThresh,
                        threshOutlayer,
                        incCarte,
                        noGrowthThresh,
                    )
                    # [xPVD,yPVD,flagPVD] = getfitData(df, j, 'PVD',startTime,timeThresh,threshOutlayer, incCarte,noGrowthThresh)

            else:
                [x, y, flag] = getfitData(
                    df,
                    j,
                    channel,
                    startTime,
                    timeThresh,
                    threshOutlayer,
                    incCarte,
                    noGrowthThresh,
                )
                # [xPVD,yPVD,flagPVD] = getfitData(df, j, 'PVD',startTime,timeThresh,threshOutlayer, incCarte,noGrowthThresh)

            if flag == True:
                idx = np.where(np.isinf(y))
                y = np.delete(y, idx[0], 0)
                x = np.delete(x, idx[0], 0)

                try:
                    # redefine bounds and run inference
                    # hyperparameters
                    # 0: amplitude
                    # 1: flexibility
                    # 2: measurement error
                    # default sets the boundaries for the first hyperparameter to be 10^-1 and 10^4 and the boundaries for the third hyperparameter to be 10^2 and 10^6
                    # b= {0: [-1, 4], 1: [-1,2], 2: [-2, 6]}
                    # initial
                    b = {0: [-1, 8], 1: [-2, 5], 2: [-3, 2]}  # paper
                    # b= {0: [-1, 8], 1: [-2,2], 2: [-3,2]}

                    # apply on the range from start to max time

                    q = fitderiv(x, y, bd=b, logs=False, showstaterrors=False)

                    if display == True:
                        plt.ioff()
                        # plot results
                        # plotDerivCurve(j, [np.column_stack((x,x,x)), xPVD ], [y, yPVD[:,0]], [2*(y[:,1]-y[:,0]), 2*(yPVD[:,1]-yPVD[:,0])] , q, startTime, timeThresh, pathFitPlot, dataFile, channel, ymin,ymax, incCarte, content, fromLag, display)
                        plotDerivCurve(
                            j,
                            np.column_stack((x, x, x)),
                            y,
                            [2 * (y[:, 1] - y[:, 0])],
                            q,
                            startTime,
                            timeThresh,
                            pathFitPlot,
                            dataFile,
                            channel,
                            ymin,
                            ymax,
                            calibParam,
                            incCarte,
                            content,
                            fromLag,
                            display,
                        )

                    # export results
                    statdict = q.printstats()
                    # with open(folder+'resultIndiv/'+dataFile+'resultIndiv_Drop'+str(j)+channel+'.csv','w') as csv_file:
                    #     writer = csv.writer(csv_file)
                    #     for key, value in statdict.items():
                    #         writer.writerow([key, value])

                    statdict["drop"] = j

                    combined_csv = pd.concat([combined_csv, pd.DataFrame(statdict, index=[0])])
                    # print(combined_csv)

                except Exception as inst:
                    print(folder)
                    print("bug od for drop: " + str(j))
                    print(type(inst))  # the exception instance
                    print(inst.args)  # arguments stored in .args
                    _, _, tb = sys.exc_info()
                    traceback.print_tb(tb)  # Fixed format
                    tb_info = traceback.extract_tb(tb)
                    filename, line, func, text = tb_info[-1]
                    print(
                        "An error occurred on line {} in statement {}".format(
                            line, text
                        )
                    )

            else:
                print(
                    "outlayer not fitted drop="
                    + str(j)
                    + " flag="
                    + str(flag)
                    + " label= "
                    + content
                )
                combined_csv = pd.concat(
                    [
                        combined_csv,
                        pd.DataFrame(
                            {
                                "max df": np.nan,
                                "max df std": np.nan,
                                "max df stderr": np.nan,
                                "time of max df": np.nan,
                                "time of max df std": np.nan,
                                "time of max df stderr": np.nan,
                                "inverse max df": np.nan,
                                "inverse max df std": np.nan,
                                "inverse max df stderr": np.nan,
                                "max y": np.nan,
                                "max y std": np.nan,
                                "max y stderr": np.nan,
                                "lag time": np.nan,
                                "lag time std": np.nan,
                                "lag time stderr": np.nan,
                                "drop": j,
                            },
                            index=[0],
                        ),
                    ]
                )
                if not fromLag:
                    with open(folder + "outliersbyAlgo_" + channel + ".csv", "a") as f:
                        writer = csv.writer(f)
                        writer.writerow([str(j)])

    if fromLag:
        combined_csv.to_csv(
            folder
            + "resultIndiv/"
            + dataFile
            + "resultIndivFitFromDetectionTime_Drop_"
            + channel
            + ".csv"
        )
    else:
        combined_csv.to_csv(
            folder + "resultIndiv/" + dataFile + "resultIndiv_Drop_" + channel + ".csv"
        )


def plotDerivCurve(
    j,
    x,
    y,
    ystd,
    q,
    startTime,
    timeThresh,
    pathFitPlot,
    dataFile,
    channel,
    ymin,
    ymax,
    calibParam,
    incCarte,
    content,
    fromLag,
    display=False,
):
    # x = xlist[0]
    # y = ylist[0]
    # ystd = ystdlist[0]

    # xPVD = xlist[1]
    # yPVD = ylist[1]
    # ystdPVD = ystdlist[1]

    if channel == "RFP":
        c = "r"
    if channel == "GFP":
        c = "b"
    if channel == "PVD":
        c = "g"
    if channel == "PVDoverGFP":
        c = "b"
    if channel == "SpeedSize":
        c = "m"

    mkr = ["o", "+", "x"]
    fig, ax = plt.subplots()
    # ax.errorbar(x,y,yerr=ystd,fmt='+'+c)
    calibParam = [float(calibParam[0]), float(calibParam[1])]
    y = (y * calibParam[0] + calibParam[1]) / np.log(10)

    ax.plot(x[:, 0], y[:, 0], ".-", color=c)
    ax.fill_between(x[:, 0], y[:, 1], y[:, 2])
    # for i in range(3):
    #    ax.plot(x[:,i],y[:,i], marker = mkr[i], color = c, linestyle = '')
    # ax.errorbar(xPVD,yPVD,yerr=ystdPVD,fmt='+g')
    # q.plotfit('f', color = 'r')
    qf = (q.f * calibParam[0] + calibParam[1]) / np.log(10)
    # ax.plot(q.t, qf,'-r')

    ax.set_title(
        dataFile
        + "\n start fit after detection threshold= "
        + str(incCarte)
        + "\n drop="
        + str(j)
    )
    ax.set_xlim(0, timeThresh)
    ax.set_xticks(np.arange(0, timeThresh, 10))
    # ax.set_xticks(np.arange(0,timeThresh,1), minor=True)
    # ax.set_yticks(np.arange(-1,10,1), minor=True)
    ax.set_yticks(np.arange(ymin, ymax))
    ax.grid(color="k")
    # ax.grid(which='minor', alpha=1, linestyle='--')
    ax.set_ylabel("log(raw fluo)")
    if not np.isnan(ymin) and not np.isnan(ymax):
        ax.set_ylim(ymin, ymax)

    # plot a vertical line on the lag time
    l = q.printstats()
    ax.plot([l["lag time"], l["lag time"]], [-2, 5], "r")

    # fig.savefig(pathFitPlot+dataFile+'fitIndiv_Drop'+str(j)+channel+'.png')
    # plt.subplot(2,1,2)
    # fig, ax=plt.subplots()
    ax2 = ax.twinx()
    q.plotfit("df", color="r")
    # ax2.errorbar(q.t, q.df, q.dfvar,fmt='-o'+c)
    ax2.set_ylabel("derivative")
    ax2.set_ylim([-0.4, 1])
    ax2.set_yticks(np.arange(-0.4, 1.5, 0.1), minor=True)
    # ax.set_title(dataFile + ' '+ j )
    # ax.set_xticks(np.arange(startTime,timeThresh,1), minor=True)
    # ax.set_xlim(startTime,timeThresh)
    ax2.grid(color="r", linestyle="--")

    if display == True:
        if fromLag:
            fig.savefig(
                pathFitPlot
                + dataFile
                + "fitAndDerivativeIndivFromDetectionTime_Drop"
                + str(j)
                + channel
                + ".png",
                bbox_inches="tight",
            )
        else:
            fig.savefig(
                pathFitPlot
                + dataFile
                + "fitAndDerivativeIndiv_Drop"
                + str(j)
                + channel
                + ".png",
                bbox_inches="tight",
            )


def barPlot(labelList, y, stdy, ymin, ymax, title, ylabel, path):
    # Build the plot
    fig, ax = plt.subplots()
    x_pos = np.arange(len(labelList))
    ax.bar(x_pos, y, yerr=stdy, align="center", alpha=0.5, ecolor="black", capsize=10)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x_pos)
    plt.xticks(rotation=90)
    ax.set_xticklabels(labelList)
    ax.set_ylim([ymin, ymax])
    ax.set_title(title)
    ax.yaxis.grid(True)

    # Save the figure and show
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight", format="png")
    plt.show()


def boxPlot(labelList, y, ymin, ymax, title, ylabel, path):
    # Build the plot
    fig, ax = plt.subplots()
    x_pos = np.arange(len(labelList))
    plt.boxplot(y)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x_pos)
    plt.xticks(rotation=90)
    ax.set_xticklabels(labelList)
    ax.set_ylim([ymin, ymax])
    ax.set_title(title)
    ax.yaxis.grid(True)

    # Save the figure and show
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight", format="png")
    plt.show()


def findIdx(your_list, item):
    lst = pd.DataFrame(your_list)
    lst2 = lst[0].values
    lst3 = lst2.tolist()
    return lst3.index(item)


def wrapFitLoop(parameters, rootPath, path, channel):
    print(path)
    parameters.to_csv(path + "parametersAnalysis.csv")
    [dropMap, df, label] = loadData(path)


    # if 'checkData' in parameters.columns and 'GFP' in parameters['checkData'].iloc[0]:
    print('path', path, 'path.split("/")', path.split("/"))
    checkfitData(
        df, dropMap, label, channel, rootPath + "/" + path.split("/")[-3], parameters
    )

    fitDataIndivParallel(df, dropMap, label, channel, path, parameters)


def getDataHistogram(label, path, channel, fromLag):
    folder = path + "resultIndiv/"
    [dropMap, df, nn] = loadData(path)

    stdgRate = pd.DataFrame()
    gRate = pd.DataFrame()
    lag = pd.DataFrame()
    stdlag = pd.DataFrame()
    yld = pd.DataFrame()
    stdyield = pd.DataFrame()
    halfTime = pd.DataFrame()

    with open(path + "outliersbyAlgo_" + channel + ".csv", "r") as f:
        reader = csv.reader(f)
        outlayerAlgo = np.squeeze(list(reader))
    if outlayerAlgo.size > 1:
        outlayerAlgo = [int(o) for o in outlayerAlgo]
    else:
        outlayerAlgo = [int(outlayerAlgo)]

    with open(path + "outliersbyEye.csv", "r") as f:
        reader = csv.reader(f)
        outlayer = np.squeeze(list(reader))
    if outlayer.size > 1:
        outlayer = [int(o) for o in outlayer]
    else:
        outlayer = [int(outlayer)]

    outlayer.extend(outlayerAlgo)

    for n, labelName in enumerate(label):
        if fromLag:
            df = pd.read_csv(
                folder
                + labelName
                + "resultIndivFitFromDetectionTime_Drop_"
                + channel
                + ".csv"
            )
            ht = pd.read_csv(folder + labelName + "halftime" + channel + ".csv")
            ht.columns = ["drop", labelName + "_ht"]
            ht = ht.set_index("drop")
        else:
            df = pd.read_csv(
                folder + labelName + "resultIndiv_Drop_" + channel + ".csv"
            )
            ht = pd.DataFrame()

        idx = df.index[df["drop"].isin(outlayer) == True].tolist()
        if idx:
            print(labelName + " \n->outlayer drop index: " + str(idx))
            df = df.drop(idx)

        df = df.set_index("drop")

        tmp = pd.DataFrame({labelName: df["max df"]})
        gRate = pd.concat([gRate, tmp], axis=1)
        halfTime = pd.concat([halfTime, ht], axis=1)
        halfTime = pd.concat([halfTime, tmp], axis=1)
        tmp = pd.DataFrame({labelName: df["lag time"]})
        lag = pd.concat([lag, tmp], axis=1)
        tmp = pd.DataFrame({labelName: df["max y"]})
        yld = pd.concat([yld, tmp], axis=1)
        tmp = pd.DataFrame({labelName: df["max df std"]})
        stdgRate = pd.concat([stdgRate, tmp], axis=1)
        tmp = pd.DataFrame({labelName: df["lag time std"]})
        stdlag = pd.concat([stdlag, tmp], axis=1)
        tmp = pd.DataFrame({labelName: df["max y std"]})
        stdyield = pd.concat([stdyield, tmp], axis=1)

    return [stdgRate, gRate, lag, stdlag, yld, stdyield, halfTime]


def poolDataInterpolate(dropMap, df, label, channel, path, parameters, savefig):
    parameters = parameters.loc[channel]

    deletionData = parameters["deletionData"]
    startTime = parameters["startTime"]
    timeThresh = parameters["timeThresh"]
    threshOutlayer = parameters["threshOutlayer"]
    incCarte = parameters["incCarte"]
    display = parameters["display"]
    fromLag = parameters["fromLag"]
    ymin = parameters["ymin"]
    ymax = parameters["ymax"]
    noGrowthThresh = parameters["noGrowthThresh"]
    nbReps = parameters["nbReps"]

    for j, tmp in enumerate(label):
        data = pd.DataFrame()

        print(tmp)
        nbAdded = 0
        p = df[0]
        x0 = np.array(p.time / 3600, dtype=np.float64)

        p = df[
            len(df) - 4
        ]  # antepenultien run pour creer le range de temps sur lequel interpoler
        xLast = np.array(p.time / 3600, dtype=np.float64)
        # x2=np.arange(.5, min(max(x0)-min(x0),max(xLast)-min(xLast)), .5)
        xinterp = np.arange(0, timeThresh, 0.5)

        data["time"] = pd.Series(xinterp).values
        data.set_index("time")

        for j, content in enumerate(dropMap[:, 1]):
            y = []
            x = []

            if content == tmp:
                [x, yyy, flag] = getfitData(
                    df,
                    j,
                    channel,
                    startTime,
                    timeThresh,
                    threshOutlayer,
                    incCarte,
                    noGrowthThresh,
                )

                if flag == True:
                    try:
                        y = yyy[:, 0]
                        f = si.interp1d(x, y)
                        x2 = np.arange(np.ceil(x[0]), np.floor(x[-1]), 0.1)
                        y2 = f(x2)
                        dtmp = pd.DataFrame()
                        dtmp["time"] = pd.Series(x2).values
                        dtmp.set_index("time")
                        dtmp[str(j)] = pd.Series(y2).values

                        data = pd.merge(data, dtmp, how="outer", on="time")

                        nbAdded += 1

                    except Exception as inst:
                        print("poolData() remove this data")
                        _, _, tb = sys.exc_info()
                        traceback.print_tb(tb)  # Fixed format
                        tb_info = traceback.extract_tb(tb)
                        filename, line, func, text = tb_info[-1]
                        print(
                            "An error occurred on line {} in statement {}".format(
                                line, text
                            )
                        )

                    # plot the growth curve with the fit and the derivative
                    if savefig:
                        pathgrowth = path + "plotInterp/"
                        if not os.path.exists(pathgrowth):
                            os.makedirs(pathgrowth)

                        # plt.plot(x2,y2,marker='o', linestyle='-',color='blue')
                        plt.plot(x, y, marker="o", linestyle="-", color="blue")

                        # ax1.set_ylim([np.log(low), np.log(high)])
                        # ax1.set_yticks(range(np.int(np.log(low)), np.int(np.log(high))))
                        plt.ylabel("log(fluo)")
                        plt.xlabel("time (h)")

                        plt.title(tmp + " drp" + str(j) + " area")
                        plt.savefig(
                            pathgrowth + "drp" + str(j) + "_deriveAll_nbpt" + "_" + tmp,
                            format="png",
                            bbox_inches="tight",
                        )
                        plt.close()
                else:
                    print("drop " + str(j) + "is an outliers")

        # print(data)
        print(path + tmp + channel + "Interp.csv")
        data.to_csv(path + tmp + channel + "Interp.csv", index=False)


def checkDataInterpolated(folder, label, channel):
    for dataFile in label:
        file = folder + dataFile + channel + "Interp.csv"
        data = pd.read_csv(file)
        # load data
        t = np.array(data["time"])
        fig = plt.figure()

        if channel == "GFP":
            c = "b"
        if channel == "RFP":
            c = "r"
        if channel == "PVD":
            c = "g"
        if channel == "PVDoverGFP":
            c = "k"

        for column in data:
            if column != "time":
                od = np.array(data[column])
                plt.scatter(t, od, color=c)
                plt.title(dataFile)
        #        plt.ylim([-10,-2])
        plt.xlabel("time (h)")
        plt.ylabel("log(fluo)")
        fig.savefig(
            folder + dataFile + "checkDataInterp" + channel + ".png",
            bbox_inches="tight",
        )
        plt.close()


def plot_distribution(ax, key, base, val, m_path):
    if key == 5 or key == 4.5 or key == 4 or key == 14:
        if "2018-02-21_dilution_CAA_WT-stock2" in m_path:
            ax.plot(base[:-1], val, "-*k", markerSize=10)
        else:
            ax.plot(base[:-1], val, "-k")
    elif key == 50 or key == 25 or key == 6.5 or key == 6:
        ax.plot(base[:-1], val, "-b")
    elif key == 500 or key == 100 or key == 24 or key == 19:
        ax.plot(base[:-1], val, "-r")
    elif key == 3500 or key == 1000 or key == 52 or key == 43:
        ax.plot(base[:-1], val, "-g")
    else:
        print("not printed:" + str(key))


def my_legend(rootPath):
    if "Exp_Incubation_4h30_6h30_24h_43h13_1bactPerDrop" in rootPath:
        lgd = "black=4h30 \nblue=6h30 \nred=24h \ngreen=43h"

    elif "Exp_Incub_4h30-6h30-24h-52h" in rootPath:
        lgd = "black=4h30 \nblue=6h30 \nred=24h \ngreen=52h"

    elif "Exp_DilutionWTstock2" in rootPath:
        lgd = "blue=25 \nred=100 \ngreen=1000"

    elif "Exp_Dilution_5_50_500_3500" in rootPath:
        lgd = "black=5 \nblue=50 \nred=500 \ngreen=3500"

    elif "Exp_Incubation_14h-19h-43h" in rootPath:
        lgd = "black=14h \nred=19h \ngreen=43h"

    elif "Exp_Dilution_M9glycerol" in rootPath:
        lgd = "black=1 \nblue=10 \nred=100 \ngreen=1000"

    elif "Exp_Dilution_PvdS" in rootPath:
        lgd = "black=5 \nblue=50 \nred=500 \ngreen=3500"

    else:
        lgd = "No legend known"

    return lgd


def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction="out")
    ax.xaxis.set_ticks_position("bottom")
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel("Sample name")


def getmarker(l):
    if l == "CAA-WT-5" or l == "SBW25-WT_CAA_1" or l == "SBW25-WT-M9gly-1":
        return 0
    if l == "CAA-WT-50" or l == "SBW25-WT_CAA_4" or l == "SBW25-WT-M9gly-4":
        return 1
    if l == "CAA-WT-500" or l == "SBW25-WT_CAA_16" or l == "SBW25-WT-M9gly-16":
        return 2
    if l == "CAA-WT-3500" or l == "SBW25-WT_CAA_64" or l == "SBW25-WT-M9gly-64":
        return 3
    if l == "CAA-WT-0.5" or l == "SBW25-WT_CAA_256" or l == "SBW25-WT-M9gly-256":
        return 4
    if l == "CAA-WT-350" or l == "SBW25-WT_CAA_1024" or l == "SBW25-WT-M9gly-1024":
        return 5


def getThresHalfTime(rootPath, source, label):
    if (
        rootPath == source + "Exp_Dilution_1_4_16_64_256_1024/CAA/"
        or rootPath == source + "Exp_Dilution_1_4_16_64_256_1024/CAAnoGly/"
        or rootPath == source + "Exp_Dilution_1_4_16_64_256_1024/CAA1/"
        or rootPath == source + "Exp_Dilution_1_4_16_64_256_1024/AyshaCAA1/"
        or rootPath == source + "Exp_Dilution_1_4_16_64_256_1024/CAAspp/"
    ):
        label.sort(key=lambda x: float(x.rsplit("_")[-1]))
        getValN0 = lambda l: l.split("_")[-1]
        thresh = 3
    if rootPath == source + "Exp_Dilution_1_4_16_64_256_1024/CAAPvdS/":
        label.sort(key=lambda x: float(x.rsplit("_")[-1]))
        getValN0 = lambda l: l.split("_")[-1]
        thresh = (np.log(64000 / 0.4e-3) - 16.83) / 0.81
    if (
        rootPath == source + "Exp_Dilution_1_4_16_64_256_1024/M9gly"
        or rootPath == source + "Exp_Dilution_1_4_16_64_256_1024/M9glyFromCAA/"
    ):
        label.sort(key=lambda x: float(x.rsplit("-")[-1]))
        getValN0 = lambda l: l.split("-")[-1]
        thresh = 4.5
    if rootPath == source + "Exp_Dilution_1_4_16_64_256_1024/SMM/":
        label.sort(key=lambda x: float(x.rsplit("_")[-1]))
        getValN0 = lambda l: l.split("-")[-1]
        thresh = 4.5
    if rootPath == source + "Exp_Dilution_5_50_500_3500/":
        label.sort(key=lambda x: float(x.rsplit("-")[-1]))
        getValN0 = lambda l: l.split("-")[-1]

    if (
        rootPath == source + "Exp_BipyCAA_WT/"
        or rootPath == source + "Exp_BipyCAA_PvdS/"
    ):
        thresh = 3

    if rootPath == source + "Exp_GreenRed/Maelle/Azur3":
        thresh = 3

    else:
        thresh = 3

    return thresh


def getValueFromString(av, ap, x):
    xxx = []
    if av != "":
        xx = x.split(av)
    else:
        xx = x

    if ap != "":
        xxx = xx[-1].split(ap)
        xxx = xxx[0]
    else:
        xxx = xx[-1]
    return xxx


def _1Lorentzian(x, amp1, cen1, wid1):
    return amp1 * wid1**2 / ((x - cen1) ** 2 + wid1**2)

def func(p, x):
    a, b = p
    return a * x + b

def calcLag(h, g, c0, threshN):
    lag = h - np.log(threshN / c0) / g
    return lag

