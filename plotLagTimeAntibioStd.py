"""
Created on Fri Jun 26 2025

@author: maxime ardré, inès gabert
This script is used to plot the lag time and growth rateas a function of the inoculum size.
Several experiments can be plotted together: rootpathlist
"""

# In[]
# ----------Imports----------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from os import listdir
from os.path import join
from scipy.stats import mannwhitneyu, kruskal
from itertools import combinations

from functionsCleanPipeline import (
    loadData,
    getValueInLabel,
    getDataHistogram,
    poolDataInterpolate,
    checkDataInterpolated,
    getThresHalfTime,
    calcLag,
)

# ---------------- Entries ------------------
rootPathList = [
    "/Users/inesgabert/Documents/LBE/experiences/GFPmut3b_antib_nal_1machine_39",
    "/Users/inesgabert/Documents/LBE/experiences/GFPmut3b_antib_nal_1machine_39_2",
    "/Users/inesgabert/Documents/LBE/experiences/GFPmut3b_antib_nal_2machine_35"
]
source = "/Users/inesgabert/Documents/LBE/experiences/"
channelList = ["GFP"]
channel = "GFP"
missingAntibio = [
    1e2
]  # à remplir si des concentrations d'antibiotiques sont manquantes dans les plots
# =[] sinon
cInitAntib = 1000  # µg/mL
inoculum_symbols = {2: "^", 1024: "s"}
# --------------- End of entries ------------

dropVolume = 4e-4  # ml
allPlots = []
paletteLag = dict(zip(rootPathList,
        sns.color_palette("Blues", n_colors=len(rootPathList) ) ) )
paletteGRate = dict(zip(rootPathList,
        sns.color_palette("Greens", n_colors=len(rootPathList) ) ) )# get the halftime of the logistic growth to measure the lag

parameterPalettes = {
    "lag": paletteLag,
    "gRate": paletteGRate,
}

parameterYLims = {
    "lag": (0, 15),
    "gRate": (0.1, 1.1),
}


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
            print("for path in folder: path:", path)
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
        dfPlots = pd.DataFrame()
        for path in folder:
            lagTable = pd.DataFrame()
            inoculumList = []
            antibioDList = []
            [dropMap, df, labelList] = loadData(path)
       
            labelList = sorted(labelList)
            for label in labelList:
          
                [_, gRate, timeToDetection, stdLag, yld, _, halfTime, conc_halfTime] = getDataHistogram(
                labelList, path, channel, True)

                inoculum, antibioD = getValueInLabel(label, path)
                if antibioD not in antibioDList:
                    antibioDList.append(antibioD)

                c0 = inoculum / dropVolume
                if inoculum not in inoculumList:
                    inoculumList.append(inoculum)
                t = pd.DataFrame()
                t[label] = calcLag(halfTime[label], gRate[label], c0, 10**conc_halfTime[label])

                lagTable = pd.concat([lagTable, t[label]], axis=1)

            # Préparation des données pour le boxplot

            # Récupérer toutes les valeurs d'inoculum et d'antibio
            inoculumList = sorted(
                list({getValueInLabel(col, path)[0] for col in lagTable.columns})
            )
            antibioDList = sorted(
                list({getValueInLabel(col, path)[1] for col in lagTable.columns})
            )

            # Préparer les données

            for col in lagTable.columns:
                inoculum, antibioD = getValueInLabel(col, path)
                common_idx = (
                    lagTable[col]
                    .dropna()
                    .index.intersection(gRate[col].dropna().index)
                    .intersection(yld[col].dropna().index)
                )
                tempDf = pd.DataFrame(
                    {
                        "lag": lagTable[col].loc[common_idx],
                        "gRate": gRate[col].loc[common_idx],
                        "yld": yld[col].loc[common_idx],
                        "antibio_d": antibioD,
                        "inoculum": inoculum,
                        "dropId": common_idx,  # index of the data in this col
                    }
                )
                dfPlots = pd.concat([dfPlots, tempDf], ignore_index=True)

            # remove lines with extreme values
            absYldMin = 5.9
            relYldMax = dfPlots["yld"].max()
            minYield = min(absYldMin, relYldMax * 0.80)
            dfPlots = dfPlots[ (dfPlots["gRate"] < 1)
                & (dfPlots["lag"] > 0)
            ]
            if missingAntibio != []:
                for antibio_d in missingAntibio:
                    # add missing antibio_d
                    newRow = pd.DataFrame(
                        {
                            "lag": None,
                            "gRate": None,
                            "yld": 0,
                            "antibio_d": antibio_d,
                            "inoculum": inoculumList[0],
                            "dropId": None,
                        },
                        index=[0],
                    )
                    dfPlots = pd.concat(
                        [dfPlots, newRow], ignore_index=True, sort=False
                    )
        dfPlots["antibio_c"] = cInitAntib / dfPlots["antibio_d"]
        dfPlots["antibio_c"].replace([np.inf, -np.inf], 0, inplace=True)
        dfPlots["rootPath"] = rootPath
        allPlots.append(dfPlots)

allPlots = pd.concat(allPlots, ignore_index=True)

# In[]
n_root = len(rootPathList)
n_inoc = len(inoculum_symbols)

for parameter, palette in parameterPalettes.items():
    fig, ax = plt.subplots(figsize=(12, 6))
    grouped = (
        allPlots.groupby(["antibio_c", "rootPath", "inoculum"])
        .agg(mean=(parameter, "mean"), std=(parameter, "std"))
        .reset_index()
    )

    # Pour placer les barres côte à côte, on calcule un décalage pour chaque combinaison
    antibio_c_list = sorted(grouped["antibio_c"].unique())
    n_antibio = len(antibio_c_list)
    ax.set_xlim(-0.5, n_antibio - 0.5)
    print(ax.get_xlim())
    ax.set_xticks(range(n_antibio))
    totalWidth = ax.get_xlim()[1] - ax.get_xlim()[0]
    width = totalWidth / (n_inoc*n_antibio)/2

    for i, antibio_c in enumerate(antibio_c_list):
        for j, rootPath in enumerate(rootPathList):
            for k, (inoc, marker) in enumerate(inoculum_symbols.items()):
                # Sélectionne la ligne correspondante
                row = grouped[
                    (grouped["antibio_c"] == antibio_c)
                    & (grouped["rootPath"] == rootPath)
                    & (grouped["inoculum"] == inoc)
                ]
                if not row.empty:
                    mean = row["mean"].values[0]
                    std = row["std"].values[0]
                    color = palette[rootPath]
                    # Position de la barre
                    x = i + (n_inoc + k) * width -0.5 - width/2
                    # Barre d'erreur
                    ax.errorbar(
                        x,
                        mean,
                        yerr=std,
                        fmt=marker,
                        color=color,
                        markerfacecolor="white",
                        markeredgewidth=2,
                        markersize=12,
                        ecolor=color,
                        capsize=5,
                        linewidth=2,
                    )

    # Axe X
    ax.set_xlim(-0.5, len(antibio_c_list) - 0.5)
    ax.set_xticks(range(len(antibio_c_list)))
    ax.set_xticklabels([str(a) for a in antibio_c_list])
    ax.set_xlabel("antibio concentration (µg/mL)")
    ax.set_ylabel("lag time (h)")
    ax.set_ylim(parameterYLims[parameter])
    ax.set_title("Lag time by antibio, rootPath (color), inoculum (symbol)")
    for tick, label in zip(ax.get_xticks(), ax.get_xticklabels()):
        try:
            value = float(label.get_text())
        except ValueError:
            continue
        if dfPlots[dfPlots["antibio_c"] == value].dropna().empty:
            # No data: add a cross symbol
            yMax = ax.get_ylim()[1]
            yMin = ax.get_ylim()[0]
            ax.text(
                tick,
                (yMax - yMin) * 0.5 + yMin,
                "✖",
                color="red",
                fontsize=20,
                ha="center",
                va="bottom",
            )
    # Légendes personnalisées
    from matplotlib.lines import Line2D
    color_legend = [
        Line2D([0], [0], color=color, lw=8, label=path.split("/")[-1])
        for path, color in palette.items()
    ]
    symbol_legend = [
        Line2D([0], [0], marker=marker, color="black", lw=0, markersize=12, markerfacecolor="white", markeredgewidth=2, label=f"{inoc} Inoculum")
        for inoc, marker in inoculum_symbols.items()
    ]
    ax.legend(handles=color_legend + symbol_legend, loc="best", fontsize=10)





# %%
