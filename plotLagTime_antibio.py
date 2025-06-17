"""
Created on Thu Jun 5 2025

@author: maxime ardré, inès gabert
This script is used to plot the lag time as a function of the inoculum size.
"""

# ----------Imports----------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.colors as mcolors
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

# ---------------- Entries ------------------
rootPathList = ["/Users/inesgabert/Documents/LBE/experiences/GFP2_antib_amp/"]
source = "/Users/inesgabert/Documents/LBE/experiences/"
channelList = ["GFP"]
channel = "GFP"
paramList = ["lag"]  # ['gRate', 'lag', 'yld']
antibio = True
# --------------- End of entries ------------
plt.close("all")

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


# plot halftime as a function of N0 to extract lag time t_lag = t_theta - (ln theta  – ln N_0)/Lambda
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
            inoculumList = []
            if antibio:
                antibioList = []
            [dropMap, df, labelList] = loadData(path)
            [stdgRate, gRate, timeToDetection, stdLag, yld, stdYield, halfTime] = getDataHistogram(
                labelList, path, channel, True
            )
            labelList = sorted(labelList)
            for label in labelList:
                h = halfTime[[label + "_ht", label]].dropna()

                if antibio:
                    inoculum, antibioC = getValueInLabel(label, path)
                    if antibioC not in antibioList:
                        antibioList.append(antibioC)

                else:
                    inoculum = getValueInLabel(label, path)

                if inoculum < 1:
                    inoculum = 1
                c0 = inoculum / dropVolume
                if inoculum not in inoculumList:
                    inoculumList.append(inoculum)
                thresholdN = calib(
                    getThresHalfTime(rootPath, source, labelList),
                    list(map(float, ast.literal_eval(parameters.calib[channel]))),
                )
                t = pd.DataFrame()
                t[label] = calcLag(h[label + "_ht"], h[label], c0, thresholdN)

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
            if antibio:
                positions = []
                width = width = fig.get_figwidth() / (
                    len(inoculumList) * len(antibioList)
                )
                for i, inoc in enumerate(inoculumList):
                    for j, ab in enumerate(antibioList):
                        positions.append(
                            np.log(inoc) + (j - len(antibioList) / 2) * width + jitter
                        )

            else:
                positions = np.log(inoculumList) + jitter
                width = 0.2

            box = plt.boxplot(
                pointList,
                positions=positions,
                showfliers=True,
                widths=width,
                whis=0.2,
                flierprops=flierProps,
                patch_artist=True,
            )
            if antibio:
                colors = plt.cm.Blues(np.linspace(0, 1, len(antibioList)))
                for patch, pos in zip(box["boxes"], positions):
                    ab_idx = int(round((pos - np.log(inoculumList[0])) / width)) % len(
                        antibioList
                    )
                    patch.set_facecolor(colors[ab_idx])

            for item in ["boxes", "whiskers", "fliers", "medians", "caps"]:
                if not antibio:
                    plt.setp(box[item], color=colorMap[colorIdx])
                    plt.setp(box[item], alpha=0.8)
                    plt.setp(box[item], lw=2)
            jitter += 0.1
            colorIdx += 1

            plotPath = path + "plotIndiv/"
            if not os.path.exists(plotPath):
                os.makedirs(plotPath)
            lagTable.to_csv(plotPath + "result_tlagGood_" + channel + ".csv")

    # plot
    if antibio:
        legend_handles = [
            Patch(facecolor=colors[i], label=f"{antibioList[i]}")
            for i in range(len(antibioList))
        ]
        plt.legend(
            handles=legend_handles,
            loc="best",
            title="antibioC (µg/mL)",
            fontsize="small",
        )
    plt.ylabel("lag (h)")
    yMin = -9
    yMax = 9
    plt.ylim([yMin, yMax])
    gridYMinor = np.arange(yMin, yMax, 1)
    gridYMajor = np.arange(yMin, yMax, 2)
    plt.xticks(np.log(inoculumList))
    ax.set_xticklabels(list(inoculumList), rotation=45, ha="right")
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
                timeToDetection,
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
                    rawData = timeToDetection
                    yMin = 0
                    yMax = 23
                    yLabel = "time to detection (h)"
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
            if antibio:
                tuple_cols = [col for col in dfData.columns if isinstance(col, tuple)]
                str_cols = [col for col in dfData.columns if isinstance(col, str)]

                # Trie les tuples, puis concatène avec les colonnes str (ex: 'date')
                sorted_cols = natsorted(tuple_cols) + str_cols
                dfData = dfData.reindex(sorted_cols, axis=1)
                antibio_c_list = sorted(antibioList)

            else:
                dfData = dfData.reindex(natsorted(dfData.columns), axis=1)


            dfData.columns = [
                str(col) if isinstance(col, tuple) else col for col in dfData.columns
            ]
            dfDataAllDate = pd.concat([dfDataAllDate, dfData])
            dfPlot = dfDataAllDate.melt(id_vars="date", value_vars=dfData.columns)
            if antibio:
                dfPlot["antibio_c"] = dfPlot["variable"].apply(
                    lambda x: str(x).split(", ")[1].replace(")", "").replace("'", "")
                )
                date_list = sorted(dfPlot["date"].unique())
                base_colors = sns.color_palette("tab10", n_colors=len(date_list))

        dfDataAllDate.to_csv(rootPath + "_df_" + dataType + "_" + channel + ".csv")

        fig, ax = plt.subplots()
        if not antibio:
            sns.boxplot(
                x="variable", y="value", hue="date", data=dfPlot, showfliers=False
            )
            ax.set_xticks(np.arange(len(dfPlot["variable"].unique())))
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        else:
            palette = {}
            for i, date in enumerate(date_list):
                base_color = base_colors[i]
                cmap = mcolors.LinearSegmentedColormap.from_list(
                    f"date_{date}", ["white", base_color]
                )
                for j, ab in enumerate(antibio_c_list):
                    # Position dans le dégradé selon la concentration
                    color = (
                        cmap(j / (len(antibio_c_list) - 1))
                        if len(antibio_c_list) > 1
                        else base_color
                    )
                    palette[f"{date}_{ab}"] = color
            dfPlot["date_antibio"] = (
                dfPlot["date"].astype(str) + "_" + dfPlot["antibio_c"].astype(str)
            )
            sns.boxplot(
                x="variable",
                y="value",
                hue="date_antibio",
                data=dfPlot,
                palette=palette,
                showfliers=False,
            )
            xtick_labels = [
                str(v).split(",")[0].replace("(", "").replace("'", "")
                for v in dfPlot["variable"].unique()
            ]
            unique_inoc = []
            final_labels = []
            for label in xtick_labels:
                if label not in unique_inoc:
                    unique_inoc.append(label)
                    final_labels.append(label)
                else:
                    final_labels.append("")
            ax.set_xticks(np.arange(len(final_labels)))
            ax.set_xticklabels(final_labels, rotation=45, ha="right")
            from matplotlib.patches import Patch

            # Crée une légende personnalisée
            legend_handles = []
            for date in date_list:
                for ab in antibio_c_list:
                    key = f"{date}_{ab}"
                    color = palette[key]
                    date2 = date.split("_")[0]  # Extract date part
                    label = f"{date2}_{ab}"
                    legend_handles.append(Patch(facecolor=color, label=label))
            ax.legend(
                handles=legend_handles,
                loc="best",
                fontsize="small",
                title="date_antibioC (µg/mL)",
            )

        ax.set_xlabel(xLabel)
        ax.set_ylim([yMin, yMax])
        if yLabel is not None:
            ax.set_ylabel(yLabel)
        ax.set_yticks(gridYTicksMajor)
        ax.set_yticks(gridYTicksMinor, minor=True)
        ax.grid(which="both")

        fig.savefig(
            rootPath + "_boxplot_" + dataType + "_" + channel + ".pdf",
            bbox_inches="tight",
            format="pdf",
        )

        ## error bar plot
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
        if not antibio:
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
        else:
            # Ajoute la colonne antibio_c et inoculum à dfPlot2
            dfPlot2["antibio_c"] = dfPlot2["variable"].apply(
                lambda x: str(x).split(", ")[1].replace(")", "").replace("'", "")
                if "," in str(x)
                else None
            )
            dfPlot2["inoculum"] = dfPlot2["variable"].apply(
                lambda x: float(str(x).split(",")[0].replace("(", "").replace("'", ""))
                if "," in str(x)
                else np.nan
            )
            log_inoc = np.log(dfPlot2["inoculum"].unique())
            inoculumList = sorted(inoculumList)
            width = fig.get_figwidth() / (len(inoculumList) * len(antibio_c_list))
            # Mapping antibio_c to index
            ab_idx_map = {ab: j for j, ab in enumerate(antibio_c_list)}
            # Pour la légende
            from matplotlib.patches import Patch

            # Plot groupé
            for i, inoc in enumerate(inoculumList):
                for ab in antibio_c_list:
                    for date in date_list:
                        mask = (
                            (dfPlot2["inoculum"] == inoc)
                            & (dfPlot2["antibio_c"] == ab)
                            & (dfPlot2["date"] == date)
                        )
                        if mask.any():
                            key = f"{date}_{ab}"
                            color = palette[key]
                            # Position groupée
                            pos = (
                                np.log(inoc)
                                + (ab_idx_map[ab] - len(antibio_c_list) / 2) * width
                            )
                            ax.errorbar(
                                x=[pos],
                                y=dfPlot2.loc[mask, "mean"],
                                yerr=dfPlot2.loc[mask, "std"],
                                fmt="o",
                                color=color,
                                label=None,
                            )
            # Xticks groupés
            ax.set_xticks([np.log(inoc) for inoc in inoculumList])
            ax.set_xticklabels(
                [str(int(inoc)) for inoc in inoculumList], rotation=45, ha="right"
            )
            ax.legend(
                handles=legend_handles,
                loc="best",
                fontsize="small",
                title="date_antibioC (µg/mL)",
            )

        ax.set_xlabel(xLabel)
        ax.set_ylim([yMin, yMax])
        ax.set_ylabel(yLabel)

        fig.savefig(
            rootPath + "_errorbar_" + dataType + "_" + channel + ".pdf",
            bbox_inches="tight",
            format="pdf",
        )
