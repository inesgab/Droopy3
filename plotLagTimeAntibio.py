"""
Created on Fri Jun 13 2025

@author: maxime ardré, inès gabert
This script is used to box plot the lag time as a function of the inoculum size.
"""

# In[]
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
    calcLag,
)

# ---------------- Entries ------------------
rootPathList = [
    "/Users/inesgabert/Documents/LBE/experiences/GFPmut3b_antib_nal_1machine_41/"
]
source = "/Users/inesgabert/Documents/LBE/experiences/"
channelList = ["GFP"]
channel = "GFP"
plot_grate = True  # plot the growth rate
# --------------- End of entries ------------
plt.close("all")
dropVolume = 4e-4  # ml

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

    for channel in channelList:
        jitter = -0.1
        colorIdx = 0
        for path in folder:
            lagTable = pd.DataFrame()
            inoculumList = []
            antibioList = []
            [dropMap, df, labelList] = loadData(path)
            [_, gRate, timeToDetection, stdLag, yld, _, halfTime] = getDataHistogram(
                labelList, path, channel, True
            )
            labelList = sorted(labelList)
            for label in labelList:
                h = halfTime[[label + "_ht", label]].dropna()

                inoculum, antibioD = getValueInLabel(label, path)
                if antibioD not in antibioList:
                    antibioList.append(antibioD)

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

            # Préparation des données pour le boxplot
            # On veut : x = antibio, y = lag, couleur = inoculum

            # 1. Récupère toutes les valeurs d'inoculum et d'antibio
            inoculumList = sorted(
                list({getValueInLabel(col, path)[0] for col in lagTable.columns})
            )
            antibioList = sorted(
                list({getValueInLabel(col, path)[1] for col in lagTable.columns})
            )

            # 2. Prépare les données au format "long"
            dfPlots = pd.DataFrame()
            for col in lagTable.columns:
                inoc, antibioD = getValueInLabel(col, path)
                temp = pd.DataFrame(
                    {
                        "lag": lagTable[col].dropna(),
                        "gRate": gRate[col].dropna(),
                        "yld": yld[col].dropna(),
                        "antibio_d": antibioD,
                        "inoculum": inoc,
                        "dropId": lagTable[col]
                        .dropna()
                        .index,  # index of the data in this col
                    }
                )
                dfPlots = pd.concat([dfPlots, temp], ignore_index=True)

            c_init = 1000  # µg/mL
            dfPlots["antibio_c"] = c_init / dfPlots["antibio_d"]

            dfPlots["antibio_c"].replace([np.inf, -np.inf], 0, inplace=True)
            # remove lines with extreme values
            dfPlots = dfPlots[(dfPlots["yld"] > 5.9) & (dfPlots["gRate"] < 1)]
            #add one line to represent the 1E2 antibio_d
            in2row = pd.DataFrame(
                {
                    "lag": None,
                    "gRate": None,
                    "yld": 0,
                    "antibio_d": 1E2,
                    "inoculum": 2,
                    "dropId": None,
                    "antibio_c": c_init / 1E2,
                },
                index=[0]
            )
            in1024row = pd.DataFrame(
                {
                    "lag": None,
                    "gRate": None,
                    "yld": 0,
                    "antibio_d": 1E2,
                    "inoculum": 1024,
                    "dropId": None,
                    "antibio_c": c_init / 1E2,
                },
                index=[0]
            )
            dfPlots = pd.concat([dfPlots, in2row, in1024row], ignore_index=True, sort=False)

            # 3. Palette de couleurs pour les inoculums
            palette = dict(
                zip(
                    sorted(dfPlots["inoculum"].unique()),
                    sns.color_palette(
                        "tab10", n_colors=len(dfPlots["inoculum"].unique())
                    ),
                )
            )

            # 4. Boxplot seaborn
            fig, ax = plt.subplots(figsize=(10, 6))
            sns.boxplot(
                x="antibio_c",
                y="lag",
                hue="inoculum",
                data=dfPlots,
                palette=palette,
                showfliers=True,
                ax=ax,
            )
            print("rootPath:", rootPath)
            ax.set_title(
                f"Lag time in function of antibiotic concentration: {rootPath.split('/')[-2]}"
            )
            ax.set_xlabel("antibio concentration (µg/mL)")
            ax.set_ylabel("lag time (h)")
            ax.legend(title="inoculum (cell/drop)")
            for tick, label in zip(ax.get_xticks(), ax.get_xticklabels()):
                try:
                    value = float(label.get_text())
                    if value == 10:
                        # Place un symbole au-dessus de la boîte (yMax à ajuster selon tes données)
                        yMax = ax.get_ylim()[1]
                        yMin = ax.get_ylim()[0]
                        ax.text(
                            tick, (yMax-yMin)*0.5 + yMin, "✖",  # ou "*" ou "?", selon ton choix
                            color="red", fontsize=20, ha="center", va="bottom"
                        )
                except ValueError:
                    continue
            plt.tight_layout()
            fig.savefig(
                rootPath + "/boxplotLagAntibio_" + channel + ".pdf",
                format="pdf",
                bbox_inches="tight"
            )

            if plot_grate:
                # 3. Boxplot seaborn pour gRate
                fig, ax = plt.subplots(figsize=(10, 6))
                sns.boxplot(
                    x="antibio_c",
                    y="gRate",
                    hue="inoculum",
                    data=dfPlots,
                    palette=palette,
                    showfliers=True,
                    ax=ax,
                )
                ax.set_title(
                    f"gRate in function of antibiotic concentration: {rootPath.split('/')[-2]}"
                )
                ax.set_xlabel("antibio concentration (µg/mL)")
                ax.set_ylabel("gRate")
                ax.legend(title="inoculum (cell/drop)")
                for tick, label in zip(ax.get_xticks(), ax.get_xticklabels()):
                    try:
                        value = float(label.get_text())
                        if value == 10:
                            # Place un symbole au-dessus de la boîte (yMax à ajuster selon tes données)
                            yMax = ax.get_ylim()[1]
                            yMin = ax.get_ylim()[0]
                            ax.text(
                                tick, (yMax-yMin)*0.5 + yMin, "✖",  # ou "*" ou "?", selon ton choix
                                color="red", fontsize=20, ha="center", va="bottom"
                            )
                    except ValueError:
                        continue
                plt.tight_layout()
                fig.savefig(
                    rootPath + "/boxplotgRateAntibio_" + channel + ".pdf",
                    format="pdf",
                    bbox_inches="tight"
                )
# %%
