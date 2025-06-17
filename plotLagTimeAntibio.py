"""
Created on Fri Jun 13 2025

@author: maxime ardré, inès gabert
This script is used to box plot the lag time as a function of the inoculum size.
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
rootPathList = ["/Users/inesgabert/Documents/LBE/experiences/GFPmut3b_antib_nal_1machine_41/"]
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
            [_, gRate, timeToDetection, stdLag, _, _, halfTime] = getDataHistogram(
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
            inoculumList = sorted(list({getValueInLabel(col, path)[0] for col in lagTable.columns}))
            antibioList = sorted(list({getValueInLabel(col, path)[1] for col in lagTable.columns}))
            
            # 2. Prépare les données au format "long"
            df_long = pd.DataFrame()
            for col in lagTable.columns:
                inoc, antibioD = getValueInLabel(col, path)
                temp = pd.DataFrame({
                    "lag": lagTable[col].dropna(),
                    "antibio_d": antibioD,
                    "inoculum": inoc
                })
                df_long = pd.concat([df_long, temp], ignore_index=True)

            c_init = 1000  # µg/mL
            df_long["antibio_c"] = c_init / df_long["antibio_d"]

            df_long["antibio_c"] = df_long["antibio_c"].round(2)
            df_long["antibio_c"].replace([np.inf, -np.inf], 0, inplace=True)


            # 3. Palette de couleurs pour les inoculums
            palette = dict(zip(
                sorted(df_long["inoculum"].unique()),
                sns.color_palette("tab10", n_colors=len(df_long["inoculum"].unique()))
            ))

            # 4. Boxplot seaborn
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.set_ylim([-5, 25])
            sns.boxplot(
                x="antibio_c",
                y="lag",
                hue="inoculum",
                data=df_long,
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
            plt.tight_layout()
            fig.savefig(
                rootPath + "/boxplotLagAntibio_" + channel + ".pdf",
                format="pdf",
                bbox_inches="tight"
            )

if plot_grate:
    # ...après le plot du lag time...

    # ----------- Plot gRate en fonction de antibio et inoculum -----------
    # 1. Prépare les données au format "long" pour gRate
    df_gRate = pd.DataFrame()
    for col in gRate.columns:
        inoc, antibioD = getValueInLabel(col, path)
        temp = pd.DataFrame({
            "gRate": gRate[col].dropna(),
            "antibio_d": antibioD,
            "inoculum": inoc
        })
        df_gRate = pd.concat([df_gRate, temp], ignore_index=True)

    df_gRate["antibio_c"] = c_init / df_gRate["antibio_d"]
    df_gRate["antibio_c"] = df_gRate["antibio_c"].round(2)
    df_gRate["antibio_c"].replace([np.inf, -np.inf], 0, inplace=True)

    # 2. Palette de couleurs pour les inoculums (même que pour lag)
    palette_g = dict(zip(
        sorted(df_gRate["inoculum"].unique()),
        sns.color_palette("tab10", n_colors=len(df_gRate["inoculum"].unique()))
    ))

    # 3. Boxplot seaborn pour gRate
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(
        x="antibio_c",
        y="gRate",
        hue="inoculum",
        data=df_gRate,
        palette=palette_g,
        showfliers=True,
        ax=ax,
    )
    ax.set_title(
        f"gRate in function of antibiotic concentration: {rootPath.split('/')[-2]}"
    )
    ax.set_xlabel("antibio concentration (µg/mL)")
    ax.set_ylabel("gRate")
    ax.legend(title="inoculum (cell/drop)")
    plt.tight_layout()
    fig.savefig(
        rootPath + "/boxplotgRateAntibio_" + channel + ".pdf",
        format="pdf",
        bbox_inches="tight"
    )