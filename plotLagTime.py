"""
Created on Thu Jun 5 2025

@author: maxime ardré, inès gabert
This script is used to plot the lag time as a function of the inoculum size.
"""

# In[]:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ast
import matplotlib
import seaborn as sns
import os
from natsort import natsorted
from os import listdir
from os.path import join

from functionsCleanPipeline import (
    loadData,
    getValueInLabel,
    getDataHistogram,
    calib,
    poolDataInterpolate,
    checkDataInterpolated,
    getThresHalfTime,
    calcLag
)

# In[]:

# ---------------- Entries ------------------
rootPathList = ["/Users/inesgabert/Documents/LBE/experiences/RFP_inoculum/"]
source = "/Users/inesgabert/Documents/LBE/experiences/"
# --------------- End of entries ------------
plt.close("all")
slope = []
inc = []

# get the halftime of the logistic growth to measure the lag
for rootPath in rootPathList:
    folder = sorted(
        [
            join(join(rootPath, o), "analysis/")
            for o in listdir(rootPath)
            if (o[:3] == "201" or o[:3] == "EXP")
            and os.path.isdir(os.path.join(rootPath, o))
        ]
    )
    channelList = ["RFP"]
    colorIndex = -1
    calculate = True

    for channel in channelList:
        for path in folder:
            [dropMap, df, labelList] = loadData(path)
            parameters = pd.read_csv(path + "parametersAnalysis.csv")
            parameters = parameters.set_index("channelList")
            parameters["display"] = False

            if calculate:
                print(
                    "calculate the interpolation of data in a single time vector for all drop"
                )
                poolDataInterpolate(
                    dropMap, df, labelList, channel, path, parameters, savefig=False
                )  # used after to measure halftime
                checkDataInterpolated(path, labelList, channel)
            else:
                print("no calculation of interpolation")

            print("plotting the results")
            threshold = getThresHalfTime(rootPath, source, labelList)

            for label in labelList:
                filePath = path + label + channel
                print(filePath)
                df = pd.read_csv(filePath + "Interp.csv")

                time = df["time"].to_list()
                df = df.drop("time", axis=1)
                indexHalf = df[df > threshold].apply(pd.Series.first_valid_index)
                indexHalf = indexHalf.dropna()
                timeLag = [time[int(i)] for i in indexHalf.values]
                dfTime = pd.DataFrame({label: timeLag}, index=indexHalf.index)
                dfTime.to_csv(
                    path + "resultIndiv/" + label + "halftime" + channel + ".csv"
                )

# In[]:

# plot halftime as a function of N0 to extract lag time t_lag = t_theta - (ln theta  – ln N_0)/Lambda
channelList = ["RFP"]
dropVolume = 4e-4  # ml



for rootPath in rootPathList:
    folder = sorted(
        [
            join(join(rootPath, o), "analysis/")
            for o in listdir(rootPath)
            if (o[:3] == "201" or o[:3] == "EXP")
            and os.path.isdir(os.path.join(rootPath, o))
        ]
    )
    print("folder:", folder)
    fig, ax = plt.subplots()
    colorMap = matplotlib.cm.Dark2(np.linspace(0, 1, len(folder)))

    for channel in channelList:
        jitter = -0.1
        colorIdx = 0
        for path in folder:
            lagTable = pd.DataFrame()
            initialCounts = []
            [dropMap, df, labelList] = loadData(path)
            [stdgRate, gRate, lag, stdLag, yld, stdYield, halfTime] = getDataHistogram(
                labelList, path, channel, True
            )
            labelList = sorted(labelList)
            for label in labelList:
                h = halfTime[[label + "_ht", label]].dropna()
                inoculum = getValueInLabel(label, path)
                if inoculum < 1:
                    inoculum = 1

                c0 = inoculum / dropVolume
                initialCounts.append(getValueInLabel(label, path))
                thresholdN = calib(
                    getThresHalfTime(rootPath, source, labelList),
                    list(map(float, ast.literal_eval(parameters.calib[channel]))),
                )
                t = pd.DataFrame()
                t[label] = calcLag(h[label + "_ht"], h[label], c0, thresholdN)
                # print(
                #     f"tlag: {t[label]}",
                #     f"ht: {h[label + '_ht']}",
                #     f"h: {h[label]}",
                #     f"C0: {c0}",
                #     f"threshN:{thresholdN}",
                # )

                lagTable = pd.concat([lagTable, t[label]], axis=1)

            pointList = []
            for columnName, columnData in lagTable.items():
                points = columnData.values
                pointList.append(points[~(pd.isnull(points))])

            flierProps = dict(
                marker="o",
                markerfacecolor=colorMap[colorIdx],
                markersize=2,
                markeredgecolor=colorMap[colorIdx],
                alpha=0.5,
            )

            box = plt.boxplot(
                pointList,
                positions=np.log(initialCounts) + jitter,
                showfliers=True,
                widths=0.2,
                whis=0.2,
                flierprops=flierProps,
            )
            for item in ["boxes", "whiskers", "fliers", "medians", "caps"]:
                plt.setp(box[item], color=colorMap[colorIdx])
                plt.setp(box[item], alpha=0.8)
                plt.setp(box[item], lw=2)
            jitter += 0.1
            colorIdx += 1

            plotPath = path + "plotIndiv/"
            if not os.path.exists(plotPath):
                os.makedirs(plotPath)
            lagTable.to_csv(plotPath + "result_tlagGood_" + channel + ".csv")


    plt.ylabel("lag (h)")
    yMin = 0
    yMax = 9
    plt.ylim([yMin, yMax])
    gridYMinor = np.arange(yMin, yMax, 1)
    gridYMajor = np.arange(yMin, yMax, 2)
    plt.xticks(np.log(initialCounts))
    ax.set_xticklabels(list(initialCounts), rotation=45, ha="right")
    plt.xlabel("inoculum on logscale")
    ax.set_yticks(gridYMajor)
    ax.set_yticks(gridYMinor, minor=True)
    ax.grid(which="both")
    date = rootPath.split("/")
    plt.title(date[-2] + " lag " + channel)
    fig.savefig(
        rootPath + "t_lagVSN0_" + channel + ".pdf",
        format="pdf",
        dpi=500,
        bbox_inches="tight",
    )
    plt.close()
# In[]:
# plot le boxplot du lag, gRate, yield en allant chercher les bons label directement et en les ordonnant
# sauve le plot qui montre tout les replicat dans rootpath
# remove the outlayers

for rootPath in rootPathList:
    folder = sorted(
        [
            join(join(rootPath, o), "analysis/")
            for o in listdir(rootPath)
            if (o[:3] == "201" or o[:3] == "EXP")
            and os.path.isdir(os.path.join(rootPath, o))
        ]
    )

    paramList = ["lag"]  # ['gRate', 'lag', 'yld']
    channel = "RFP"  # GFP RFP PVD
    logX = True
    xLabel = "inoculum " + channel + " (cell/drop)"

    dfDataAllDate = pd.DataFrame()

    for dataType in paramList:
        for path in folder:
            print(path)
            dfData = pd.DataFrame()
            pathSplit = path.split("/")
            date = pathSplit[-3]
            rawData = []
            cleanedData = []
            listData = []

            [_, _, labelList] = loadData(path)
            [
                stdGrowthRate,
                growthRate,
                lagTime,
                stdLagTime,
                yieldValue,
                stdYield,
                halfTime,
            ] = getDataHistogram(labelList, path, channel, fromLag=True)
            for label in labelList:
                # définir le yrange
                if dataType == "gRate":
                    rawData = growthRate
                    if (
                        "/CAA/" in path
                        or "/CAAnoGly/" in path
                        or "/CAAPvdS/" in path
                        or "/CAA1/" in path
                    ):
                        yMin = 0.5
                        yMax = 2
                    if "M9gly" in path:
                        yMin = 0
                        yMax = 0.5
                    if "BipyCAA" in path:
                        yMin = 0.2
                        yMax = 1
                    gridYTicksMinor = np.arange(yMin, yMax, 0.01)
                    gridYTicksMajor = np.arange(yMin, yMax, 0.1)

                elif dataType == "lag":
                    rawData = lagTime
                    yMin = 0
                    yMax = 23
                    yLabel = "lag time (h)"
                    if "BipyCAA" in path:
                        yMin = 0
                        yMax = 35
                    gridYTicksMinor = np.arange(yMin, yMax, 0.5)
                    gridYTicksMajor = np.arange(yMin, yMax, 1)

                elif dataType == "yld":
                    rawData = yieldValue
                    yLabel = "maximal concentration (cell/ml)"
                    if (
                        "/CAA/" in path
                        or "/CAAnoGly/" in path
                        or "/CAAPvdS/" in path
                        or "/CAA1/" in path
                        or "/Maelle/" in path
                    ):
                        scale = 1e9
                        yMin = -0.1 * scale
                        yMax = 6 * scale
                        rawData = calib(
                            yieldValue,
                            list(
                                map(float, ast.literal_eval(parameters.calib[channel]))
                            ),
                        )
                        gridYTicksMinor = np.arange(yMin, yMax, scale / 2)
                        gridYTicksMajor = np.arange(yMin, yMax, scale)
                    if "M9gly" in path:
                        yMin = 5
                        yMax = 7
                    if "BipyCAA" in path:
                        yMin = 2
                        yMax = 7
                        gridYTicksMinor = np.arange(yMin, yMax, 0.1)
                        gridYTicksMajor = np.arange(yMin, yMax, 0.2)

                elif dataType == "stdgRate":
                    rawData = stdGrowthRate
                    if (
                        "/CAA/" in path
                        or "/CAAnoGly/" in path
                        or "/CAAPvdS/" in path
                        or "/CAA1/" in path
                    ):
                        yMin = 0
                        yMax = 0.06
                    if "M9gly" in path:
                        yMin = 0
                        yMax = 0.02
                    if "BipyCAA" in path:
                        yMin = 0
                        yMax = 0.06
                    gridYTicksMinor = np.arange(yMin, yMax, 0.001)
                    gridYTicksMajor = np.arange(yMin, yMax, 0.01)
                else:
                    print("do not get kind of data")

                # get the non nan data for every label
                cleanedData = [x for x in rawData[label] if str(x) != "nan"]
                listData.append(cleanedData)
                dfData[getValueInLabel(label, path)] = pd.Series(cleanedData)

            dfData["date"] = pd.Series(dfData.shape[0] * [date])
            dfData = dfData.reindex(natsorted(dfData.columns), axis=1)

            dfDataAllDate = pd.concat([dfDataAllDate, dfData])
            dfPlot = dfDataAllDate.melt(id_vars="date", value_vars=dfData.columns)

        dfDataAllDate.to_csv(rootPath + "_df_" + dataType + "_" + channel + ".csv")

        fig, ax = plt.subplots()
        sns.boxplot(x="variable", y="value", hue="date", data=dfPlot, showfliers=False)

        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        ax.set_xlabel(xLabel)
        ax.set_ylim([yMin, yMax])
        ax.set_ylabel(yLabel)
        ax.set_yticks(gridYTicksMajor)
        ax.set_yticks(gridYTicksMinor, minor=True)
        ax.grid(which="both")

        fig.savefig(
            rootPath + "_boxplot_" + dataType + "_" + channel + ".pdf",
            bbox_inches="tight",
            format="pdf",
        )

        dfPlot2 = pd.DataFrame()
        dfDataMean = pd.DataFrame()
        dfDataStd = pd.DataFrame()

        dfDataMean = dfDataAllDate.groupby("date").mean()
        dfDataMean["date"] = dfDataMean.index
        dfPlot2 = dfDataMean.melt(
            id_vars="date", value_vars=dfDataMean.columns, value_name="mean"
        )

        dfDataStd = dfDataAllDate.groupby("date").std()
        dfDataStd["date"] = dfDataStd.index
        dfDataStd = dfDataStd.melt(
            id_vars="date", value_vars=dfDataStd.columns, value_name="std"
        )
        dfPlot2["std"] = dfDataStd["std"]

        fig, ax = plt.subplots()

        uniqueDates = dfPlot2.date.unique()
        for dateValue in uniqueDates:
            ax.errorbar(
                x=dfPlot2[dfPlot2["date"] == dateValue]["variable"],
                y=dfPlot2[dfPlot2["date"] == dateValue]["mean"],
                yerr=dfPlot2[dfPlot2["date"] == dateValue]["std"],
                label=dateValue,
                linestyle="",
                fmt="o",
            )

        ax.legend()
        ax.set_xlabel(xLabel)
        ax.set_ylim([yMin, yMax])
        ax.set_ylabel(yLabel)

        fig.savefig(
            rootPath + "_errorbar_" + dataType + "_" + channel + ".pdf",
            bbox_inches="tight",
            format="pdf",
        )

# %%
