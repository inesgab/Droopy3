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
import ast
import seaborn as sns
import os
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
    "/Users/inesgabert/Documents/LBE/experiences/GFPmut3b_antib_nal_1machine_39/"
]
source = "/Users/inesgabert/Documents/LBE/experiences/"
channelList = ["GFP"]
channel = "GFP"
missingAntibio = [
    1e2
]  # à remplir si des concentrations d'antibiotiques sont manquantes dans les plots
# =[] sinon
cInitAntib = 1000  # µg/mL
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
        for path in folder:
            lagTable = pd.DataFrame()
            inoculumList = []
            antibioDList = []
            [dropMap, df, labelList] = loadData(path)
            [_, gRate, timeToDetection, stdLag, yld, _, halfTime] = getDataHistogram(
                labelList, path, channel, True
            )
            labelList = sorted(labelList)
            for label in labelList:
                h = halfTime[[label + "_ht", label]].dropna()

                inoculum, antibioD = getValueInLabel(label, path)
                if antibioD not in antibioDList:
                    antibioDList.append(antibioD)

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

            # Récupérer toutes les valeurs d'inoculum et d'antibio
            inoculumList = sorted(
                list({getValueInLabel(col, path)[0] for col in lagTable.columns})
            )
            antibioDList = sorted(
                list({getValueInLabel(col, path)[1] for col in lagTable.columns})
            )

            # Préparer les données
            dfPlots = pd.DataFrame()
            for col in lagTable.columns:
                inoculum, antibioD = getValueInLabel(col, path)
                temp = pd.DataFrame(
                    {
                        "lag": lagTable[col].dropna(),
                        "gRate": gRate[col].dropna(),
                        "yld": yld[col].dropna(),
                        "antibio_d": antibioD,
                        "inoculum": inoculum,
                        "dropId": lagTable[col]
                        .dropna()
                        .index,  # index of the data in this col
                    }
                )
                dfPlots = pd.concat([dfPlots, temp], ignore_index=True)

            # remove lines with extreme values
            absYldMin = 5.9
            relYldMax = dfPlots["yld"].max()
            minYield = min(absYldMin, relYldMax * 0.90)
            dfPlots = dfPlots[
                (dfPlots["yld"] > minYield)
                & (dfPlots["gRate"] < 1)
                & (dfPlots["lag"] > 0)
            ]
            if missingAntibio != []:
                for antibio_d in missingAntibio:
                    # on ajoute une fausse ligne pour plot l'absence de données pour une concentration donnée
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
        # 3. Palette de couleurs pour les inoculums
        palette = dict(
            zip(
                sorted(dfPlots["inoculum"].unique()),
                sns.color_palette("husl", n_colors=len(dfPlots["inoculum"].unique())),
            )
        )
        palette2 = dict(
            zip(
                sorted(dfPlots["inoculum"].unique()),
                sns.color_palette("hls", n_colors=len(dfPlots["inoculum"].unique())),
            )
        )

        # Boxplot seaborn pour lagtime
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

        # Pour chaque valeur, vérifie s'il y a des données, et met une croix si pas de valeurs
        for tick, label in zip(ax.get_xticks(), ax.get_xticklabels()):
            try:
                value = float(label.get_text())
            except ValueError:
                continue
            # Vérifie s'il y a des données pour cette concentration
            if dfPlots[dfPlots["antibio_c"] == value].dropna().empty:
                # Pas de données : ajoute un symbole
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
        plt.tight_layout()
        fig.savefig(
            rootPath + "/boxplotLagAntibio_" + channel + ".pdf",
            format="pdf",
            bbox_inches="tight",
        )

        # Boxplot seaborn pour gRate
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.boxplot(
            x="antibio_c",
            y="gRate",
            hue="inoculum",
            data=dfPlots,
            palette=palette2,
            showfliers=True,
            ax=ax,
        )
        ax.set_title(
            f"gRate in function of antibiotic concentration: {rootPath.split('/')[-2]}"
        )
        ax.set_xlabel("antibio concentration (µg/mL)")
        ax.set_ylabel("gRate")
        ax.legend(title="inoculum (cell/drop)")
        # met une croix si pas de valeurs
        for tick, label in zip(ax.get_xticks(), ax.get_xticklabels()):
            try:
                value = float(label.get_text())
            except ValueError:
                continue
            # Vérifie s'il y a des données pour cette concentration
            if dfPlots[dfPlots["antibio_c"] == value].dropna().empty:
                # Pas de données : ajoute un symbole
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
        plt.tight_layout()
        fig.savefig(
            rootPath + "/boxplotgRateAntibio_" + channel + ".pdf",
            format="pdf",
            bbox_inches="tight",
        )
# %%
# Extraction des deux Series pour inoculum = 1024 et inoculum = 2
series_1024 = dfPlots[dfPlots["inoculum"] == 1024]
series_2 = dfPlots[dfPlots["inoculum"] == 2]

# %%
