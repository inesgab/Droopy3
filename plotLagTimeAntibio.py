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
    "/Users/inesgabert/Documents/LBE/experiences/GFPmut3b_antib_NaHCO3_1machine_39/"
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
        # inoculum palette
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

        # seaborn boxplot for lag time
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
        ax.set_ylim([5, 17])

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

        grouped = (
            dfPlots.groupby(["antibio_c", "inoculum"])["lag"].count().reset_index()
        )
        xticks = ax.get_xticks()
        xticklabels = [float(lbl.get_text()) for lbl in ax.get_xticklabels()]
        hue_levels = sorted(dfPlots["inoculum"].unique())
        n_hue = len(hue_levels)
        width = 0.8 / n_hue  # largeur de chaque box

        for i, x_tick in enumerate(xticks):
            for j, hue in enumerate(hue_levels):
                x = x_tick - 0.4 + width / 2 + j * width
                count = grouped[
                    (grouped["antibio_c"] == xticklabels[i])
                    & (grouped["inoculum"] == hue)
                ]["lag"]
                if not count.empty:
                    # max of lag for this group: position of the text
                    y_max = dfPlots[
                        (dfPlots["antibio_c"] == xticklabels[i])
                        & (dfPlots["inoculum"] == hue)
                    ]["lag"].max()
                    ax.text(
                        x,
                        y_max + 0.1,
                        f"n={int(count.values[0])}",
                        ha="center",
                        va="bottom",
                        fontsize=8,
                        color="black",
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
        # cross for missing data
        for tick, label in zip(ax.get_xticks(), ax.get_xticklabels()):
            try:
                value = float(label.get_text())
            except ValueError:
                continue
            # checks if there are data for this concentration
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
        # number of samples for each group
        grouped = (
            dfPlots.groupby(["antibio_c", "inoculum"])["gRate"].count().reset_index()
        )

        # Gets the positions of the ticks (one per box)
        xticks = ax.get_xticks()
        xticklabels = [float(lbl.get_text()) for lbl in ax.get_xticklabels()]
        hue_levels = sorted(dfPlots["inoculum"].unique())
        n_hue = len(hue_levels)
        width = 0.8 / n_hue  # width of each box

        for i, x_tick in enumerate(xticks):
            for j, hue in enumerate(hue_levels):
                x = x_tick - 0.4 + width / 2 + j * width
                count = grouped[
                    (grouped["antibio_c"] == xticklabels[i])
                    & (grouped["inoculum"] == hue)
                ]["gRate"]
                if not count.empty:
                    y_max = dfPlots[
                        (dfPlots["antibio_c"] == xticklabels[i])
                        & (dfPlots["inoculum"] == hue)
                    ]["gRate"].max()
                    ax.text(
                        x,
                        y_max + 0.01,
                        f"n={int(count.values[0])}",
                        ha="center",
                        va="bottom",
                        fontsize=8,
                        color="black",
                    )
        plt.tight_layout()
        fig.savefig(
            rootPath + "/boxplotgRateAntibio_" + channel + ".pdf",
            format="pdf",
            bbox_inches="tight",
        )
# %%
# tests statistiques

# Liste pour stocker les résultats
mannWhitney = []
print("-------------- Mann Whitney test on inoculum results ----------------")

for antibio_c in sorted(dfPlots["antibio_c"].unique()):
    inoculums = sorted(dfPlots["inoculum"].unique())
    for inoc1, inoc2 in combinations(inoculums, 2):
        # --- Lag time ---
        group1 = dfPlots[
            (dfPlots["antibio_c"] == antibio_c) & (dfPlots["inoculum"] == inoc1)
        ]["lag"].dropna()
        group2 = dfPlots[
            (dfPlots["antibio_c"] == antibio_c) & (dfPlots["inoculum"] == inoc2)
        ]["lag"].dropna()
        if len(group1) > 0 and len(group2) > 0:
            stat, p = mannwhitneyu(group1, group2, alternative="two-sided")
            mannWhitney.append(
                {
                    "measure": "lag",
                    "inoculum_1": inoc1,
                    "inoculum_2": inoc2,
                    "antibio_c": antibio_c,
                    "p_value": p,
                }
            )
            if p < 0.05:
                print(
                    f"[Lag] Distributions different for inoculum {inoc1} vs {inoc2} at antibio_c={antibio_c:.2f} µg/mL (p={p:.3g})"
                )
        # --- gRate ---
        group1 = dfPlots[
            (dfPlots["antibio_c"] == antibio_c) & (dfPlots["inoculum"] == inoc1)
        ]["gRate"].dropna()
        group2 = dfPlots[
            (dfPlots["antibio_c"] == antibio_c) & (dfPlots["inoculum"] == inoc2)
        ]["gRate"].dropna()
        if len(group1) > 0 and len(group2) > 0:
            stat, p = mannwhitneyu(group1, group2, alternative="two-sided")
            mannWhitney.append(
                {
                    "measure": "gRate",
                    "inoculum_1": inoc1,
                    "inoculum_2": inoc2,
                    "antibio_c": antibio_c,
                    "p_value": p,
                }
            )
            if p < 0.05:
                print(
                    f"[gRate] Distributions different for inoculum {inoc1} vs {inoc2} at antibio_c={antibio_c:.2f} µg/mL (p={p:.3g})"
                )


# Convertir en DataFrame
mannWhitneyDf = pd.DataFrame(mannWhitney)

# Afficher ou sauvegarder
print(mannWhitneyDf)

kruskalResults = []
print(
    "\n \n-------------- Kruskal-Willis test on all antibioC for each inoculum --------------"
)
for inoc in [2, 1024]:
    # --- lag ---
    groups = [
        dfPlots[(dfPlots["inoculum"] == inoc) & (dfPlots["antibio_c"] == ab)][
            "lag"
        ].dropna()
        for ab in sorted(dfPlots["antibio_c"].unique())
    ]
    # Garde seulement les groupes non vides
    groups = [g for g in groups if len(g) > 0]
    if len(groups) > 1:
        stat, p = kruskal(*groups)
        kruskalResults.append({"measure": "lag", "inoculum": inoc, "p_value": p})
        if p < 0.05:
            print(
                f"[lag] Kruskal-Wallis: distributions differ between antibio groups for inoculum {inoc} (p={p:.3g})"
            )

    # --- gRate ---
    groups = [
        dfPlots[(dfPlots["inoculum"] == inoc) & (dfPlots["antibio_c"] == ab)][
            "gRate"
        ].dropna()
        for ab in sorted(dfPlots["antibio_c"].unique())
    ]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) > 1:
        stat, p = kruskal(*groups)
        kruskalResults.append({"measure": "gRate", "inoculum": inoc, "p_value": p})
    if p < 0.05:
        print(
            f"[gRate] Kruskal-Wallis: distributions differ between antibio groups for inoculum {inoc} (p={p:.3g})"
        )

kruskalDf = pd.DataFrame(kruskalResults)
print(kruskalDf)
# kruskal_df.to_csv("kruskal_results.csv", index=False)

# test de Kruskal-Wallis pour comparer les groupes
# mannWhitney.to_csv("mannwhitney_mannWhitney.csv", index=False)
# %%
