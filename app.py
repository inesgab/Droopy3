#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#In[]
"""
Created on Fri Apr 16 13:42:40 2021

@author: Maxime Ardré

"""

import shutil
import numpy as np
import csv
import pandas as pd
from os.path import join
from os import listdir
from natsort import natsorted
from PyQt5 import QtCore, QtWidgets, uic
from PyQt5.QtWidgets import (
    QApplication,
    QTabWidget,
    QFileDialog,
    QMessageBox,
)
from PyQt5.QtGui import QPixmap
import ast
import os, glob, sys
import datetime

os.chdir(
    os.path.dirname(os.path.abspath(__file__))
)  # Change to the directory of this script
module_path = os.path.abspath(os.path.join("./"))
if module_path not in sys.path:
    sys.path.append(module_path)
from functionsCleanPipeline import loadData, checkfitDataCalib, wrapFitLoop


qtcreator_file = "./ui/Droopy.ui"
Ui_TabWidget, QtBaseClass = uic.loadUiType(qtcreator_file)


class TabWidget(QTabWidget, Ui_TabWidget):
    def __init__(self):
        QtWidgets.QTabWidget.__init__(self)
        Ui_TabWidget.__init__(self)
        self.setupUi(self)
        self.buttonSetParam.clicked.connect(self.setParam)
        self.loadButton.clicked.connect(self.readCsv)
        self.buttonCheckData.clicked.connect(self.checkData)
        self.spinBoxSelectPlotCheckData.valueChanged.connect(self.loadImage)
        self.pushButtonStartFitting.clicked.connect(self.firstFit)
        self.pushButtonRootPath.clicked.connect(self.selectFolder)
        self.msgBox = QMessageBox()
        self.buttonCheckOutliers.clicked.connect(self.checkOutliers)
        self.pushButtonGood.clicked.connect(lambda: self.writeOutliers(False))
        self.pushButtonOutliers.clicked.connect(lambda: self.writeOutliers(True))
        self.buttonCheckOutliers.clicked.connect(self.setFocus)
        self.pushButtonStartFitting_2.clicked.connect(self.secondFit)
        self.spinBoxSelectCheckFit.valueChanged.connect(self.displayFit)
        self.spinBoxSelectCheckOutlier.valueChanged.connect(self.getBackOutlier)
        self.ButtonGoBackOutliers.clicked.connect(
            lambda: self.spinBoxSelectCheckOutlier.setValue(
                self.spinBoxSelectCheckOutlier.value() - 1
            )
        )
        self.loadedFit = False

    def selectFolder(self):
        name = QFileDialog.getExistingDirectory(self, "Get Dir PAth")
        self.rootPathEdit.setText(name)

        # check, if selected Folder is not directly "Exp-Folder", If yes, give a message that another folder should be selected
        analysis_path = os.path.join(name, "analysis")

        if os.path.exists(analysis_path) and os.path.isdir(analysis_path):
            self.msgBox.setIcon(QMessageBox.Warning)
            self.msgBox.setText(
                "The selected folder contains directly the subfolder 'analysis'. Please select the folder before (containing the experiment folder)"
            )
            self.msgBox.setWindowTitle("wrong folder selection")
            self.msgBox.setStandardButtons(QMessageBox.Ok)

            self.msgBox.buttonClicked.connect(self.handleMessageBoxResponse)
            self.msgBox.exec_()
        else:
            print("Folder to analyze experiment(s) was selected")

    def handleMessageBoxResponse(self, button):
        # Aktionen basierend auf der Benutzerantwort
        if button.text() == "Ok":
            print("new folder can be selected now")
            # Hier können Sie weitere Aktionen für den Fall durchführen, dass der Benutzer "Ok" klickt.

    def readCsv(self):
        if self.loadParam.toPlainText() == "":
            name = QFileDialog.getOpenFileName()
            self.loadParam.setText(name[0])

        self.parameters = pd.read_csv(self.loadParam.toPlainText())
        self.parameters = self.parameters.set_index("channelList")

        self.rootPath = self.parameters.rootPath[0]

        model = DataFrameModel(self.parameters)
        self.tableView.setModel(model)

        for choice in self.parameters.index.values:
            if choice == "GFP":
                self.GFP.setChecked(True)
                self.maxTimeGFP.setText(str(self.parameters["timeThresh"][choice]))
                self.yieldNoGrowthThreshGFP.setText(
                    str(self.parameters["noGrowthThresh"][choice])
                )
                self.maxTimeNoGrowthThreshGFP.setText(
                    str(self.parameters["threshOutlayer"][choice])
                )
                self.startTimeGFP.setText(str(self.parameters["startTime"][choice]))
                calibVal = ast.literal_eval(self.parameters["calib"][choice])
                self.calibAGFP.setText(str(calibVal[0]))
                self.calibBGFP.setText(str(calibVal[1]))
                self.parameters["calib"][choice] = [str(calibVal[0]), str(calibVal[1])]
                self.yMinPlot.setText(str(self.parameters["ymin"][choice]))
                self.yMaxPlot.setText(str(self.parameters["ymax"][choice]))
                self.incCarte.setText(str(self.parameters["incCarte"][choice]))
                self.rootPathEdit.setText(self.parameters["rootPath"][choice])
            if choice == "RFP":
                self.RFP.setChecked(True)
                self.maxTimeRFP.setText(str(self.parameters["timeThresh"][choice]))
                self.yieldNoGrowthThreshRFP.setText(
                    str(self.parameters["noGrowthThresh"][choice])
                )
                self.maxTimeNoGrowthThreshRFP.setText(
                    str(self.parameters["threshOutlayer"][choice])
                )
                self.startTimeRFP.setText(str(self.parameters["startTime"][choice]))
                calibVal = ast.literal_eval(self.parameters["calib"][choice])
                self.calibARFP.setText(str(calibVal[0]))
                self.calibBRFP.setText(str(calibVal[1]))
                self.parameters["calib"][choice] = [str(calibVal[0]), str(calibVal[1])]
                self.yMinPlot.setText(str(self.parameters["ymin"][choice]))
                self.yMaxPlot.setText(str(self.parameters["ymax"][choice]))
                self.incCarte.setText(str(self.parameters["incCarte"][choice]))
                self.rootPathEdit.setText(self.parameters["rootPath"][choice])
            if choice == "PVD":
                self.PVD.setChecked(True)
                self.maxTimePVD.setText(str(self.parameters["timeThresh"][choice]))
                self.yieldNoGrowthThreshPVD.setText(
                    str(self.parameters["noGrowthThresh"][choice])
                )
                self.maxTimeNoGrowthThreshPVD.setText(
                    str(self.parameters["threshOutlayer"][choice])
                )
                self.startTimePVD.setText(str(self.parameters["startTime"][choice]))
                calibVal = ast.literal_eval(self.parameters["calib"][choice])
                self.calibAPVD.setText(str(calibVal[0]))
                self.calibBPVD.setText(str(calibVal[1]))
                self.parameters["calib"][choice] = [str(calibVal[0]), str(calibVal[1])]
                self.yMinPlot.setText(str(self.parameters["ymin"][choice]))
                self.yMaxPlot.setText(str(self.parameters["ymax"][choice]))
                self.incCarte.setText(str(self.parameters["incCarte"][choice]))
                self.rootPathEdit.setText(self.parameters["rootPath"][choice])

            print(self.parameters["calib"])

    def setParam(self):
        channelList = ["GFP", "RFP", "PVD"]
        choice = [self.GFP.isChecked(), self.RFP.isChecked(), self.PVD.isChecked()]
        timeThreshList = [
            self.maxTimeGFP.toPlainText(),
            self.maxTimeRFP.toPlainText(),
            self.maxTimePVD.toPlainText(),
        ]
        threshOutlayerList = [
            self.yieldNoGrowthThreshGFP.toPlainText(),
            self.yieldNoGrowthThreshRFP.toPlainText(),
            self.yieldNoGrowthThreshPVD.toPlainText(),
        ]
        noGrowthThreshList = [
            self.maxTimeNoGrowthThreshGFP.toPlainText(),
            self.maxTimeNoGrowthThreshRFP.toPlainText(),
            self.maxTimeNoGrowthThreshPVD.toPlainText(),
        ]
        yieldNoGrowthThresh = [
            self.yieldNoGrowthThreshGFP.toPlainText(),
            self.yieldNoGrowthThreshRFP.toPlainText(),
            self.yieldNoGrowthThreshPVD.toPlainText(),
        ]
        startTimeList = [
            self.startTimeGFP.toPlainText(),
            self.startTimeRFP.toPlainText(),
            self.startTimePVD.toPlainText(),
        ]
        calibList = [
            [float(self.calibAGFP.toPlainText()), float(self.calibBGFP.toPlainText())],
            [float(self.calibARFP.toPlainText()), float(self.calibBRFP.toPlainText())],
            [float(self.calibAPVD.toPlainText()), float(self.calibBPVD.toPlainText())],
        ]

        d = dict(
            {
                "rootPath": self.rootPathEdit.toPlainText(),
                "channelList": [channelList[i] for i, l in enumerate(choice) if l],
                "incCarte": float(self.incCarte.toPlainText()),
                "timeThresh": [
                    float(timeThreshList[i]) for i, l in enumerate(choice) if l
                ],
                "threshOutlayer": [
                    float(noGrowthThreshList[i]) for i, l in enumerate(choice) if l
                ],
                "display": self.displayPlot.isChecked(),
                "noGrowthThresh": [
                    float(yieldNoGrowthThresh[i]) for i, l in enumerate(choice) if l
                ],
                "startTime": [
                    float(startTimeList[i]) for i, l in enumerate(choice) if l
                ],
                "deletionData": self.deleteResults.checkState(),
                "ymin": [
                    float(self.yMinPlot.toPlainText())
                    for i, l in enumerate(choice)
                    if l
                ],
                "ymax": [
                    float(self.yMaxPlot.toPlainText())
                    for i, l in enumerate(choice)
                    if l
                ],
                "nbReps": 0,
                "calib": [calibList[i] for i, l in enumerate(choice) if l],
            }
        )

        # self.parameters = pd.DataFrame.from_dict(d,orient='index').T
        self.parameters = pd.DataFrame.from_dict(d)
        self.parameters = self.parameters.set_index("channelList")
        print(self.parameters)
        model = DataFrameModel(self.parameters)
        self.tableView.setModel(model)

    def getBackOutlier(self):
        nbFit = self.spinBoxSelectCheckOutlier.value()
        print(self.listPlotFitName)
        if nbFit >= len(self.listPlotFitName) - 1:
            nbFit = len(self.listPlotFitName) - 1
        image = self.listPlotFitName[nbFit]
        path = ""
        for x in image.split("/")[1:-3]:
            path = path + "/" + x
        path = path + "/"
        self.checkOutliersLabelPath.setText(path)
        self.plot(
            self.checkOutliersLabel, self.listPlotFitName, nbFit, self.nameCheckOutliers
        )

    def checkOutliers(self):
        for self.channel in self.parameters.index:
            self.rootPath = self.parameters["rootPath"][self.channel]
            folder = sorted(
                [
                    join(join(self.rootPath, o), "analysis/")
                    for o in listdir(self.rootPath)
                    if (o[:3] == "201" or o[:3] == "EXP")
                    and os.path.isdir(os.path.join(self.rootPath, o))
                ]
            )
            self.listPlotFitName = []
            for self.pathFolder in folder:
                pathRlt = self.pathFolder + "resultIndiv/"
                dirlist = [
                    filename
                    for filename in sorted(os.listdir(pathRlt))
                    if os.path.isdir(os.path.join(pathRlt, filename))
                ]

                if os.path.exists(self.pathFolder + "outliersbyEye.csv"):
                    now = str(datetime.datetime.now())[:19]
                    now = now.replace(":", "_")
                    shutil.move(
                        self.pathFolder + "outliersbyEye.csv",
                        self.pathFolder + "outliersbyEye_save" + now + ".csv",
                    )
                    print(
                        "copied outliersbyEye.csv as outliersbyEye_save" + now + ".csv"
                    )

                names = []
                for d in dirlist:
                    listImg = natsorted(listdir(pathRlt + "/" + d))
                    print(listImg)
                    if ".DS_Store" in listImg:
                        listImg.pop(listImg.index(".DS_Store"))
                    listImg = [Im for Im in listImg if "FromDetectionTime" not in Im]

                    names.extend([pathRlt + d + "/" + l for l in listImg])

                fileOutliersAlgo = pd.DataFrame([])
                if os.path.exists(
                    self.pathFolder + "outliersbyAlgo_" + self.channel + ".csv"
                ):
                    fileOutliersAlgo = pd.read_csv(
                        self.pathFolder + "outliersbyAlgo_" + self.channel + ".csv",
                        header=None,
                    )

                for imName in names:
                    drop = self.getDropFromName(imName)
                    if drop not in fileOutliersAlgo:
                        self.listPlotFitName = self.listPlotFitName + [imName]

            self.plot(
                self.checkOutliersLabel, self.listPlotFitName, 0, self.nameCheckImage
            )

    def getDropFromName(self, imName):
        ii = imName.split("Drop")
        iii = ii[1].split(self.channel)
        drop = int(iii[0])
        return drop

    def writeOutliers(self, outliers):
        self.checkOutliersLabelPath.setText(self.pathFolder)
        self.checkDataLabelChannel_2.setText(self.channel)
        nbFit = self.spinBoxSelectCheckOutlier.value()
        if outliers:
            drop = self.getDropFromName(self.listPlotFitName[nbFit])
            image = self.listPlotFitName[nbFit]
            path = ""
            for x in image.split("/")[1:-3]:
                path = path + "/" + x
            path = path + "/"
            self.dropFitLabel.setText("drop=" + str(drop) + " is an outlayer")
            with open(path + "outliersbyEye.csv", "a") as f:
                writer = csv.writer(f)
                writer.writerow([drop])

        self.progressBarOutliers.setValue(
            np.int_(100 * nbFit / len(self.listPlotFitName))
        )
        app.processEvents()
        nbFit += 1
        self.spinBoxSelectCheckOutlier.setValue(nbFit)
        self.plot(
            self.checkOutliersLabel, self.listPlotFitName, nbFit, self.nameCheckOutliers
        )

    def checkData(self):
        for channel in self.parameters.index:
            self.rootPath = self.parameters["rootPath"][channel]
            folder = sorted(
                [
                    join(join(self.rootPath, o), "analysis/")
                    for o in listdir(self.rootPath)
                    if (o[:3] == "201" or o[:3] == "EXP")
                    and os.path.isdir(os.path.join(self.rootPath, o))
                ]
            )

            for path in folder:
                print(path)
                [dropMap, df, label] = loadData(path)
                # checkfitData(df, dropMap, label, channel, self.rootPath+'/'+path.split('/')[-3],self.parameters)
                checkfitDataCalib(
                    df,
                    dropMap,
                    label,
                    channel,
                    self.rootPath + "/" + path.split("/")[-3],
                    self.parameters,
                )

        self.loadImage()

    def plot(self, view, listImage, nbImage, labelName):
        if nbImage >= len(listImage) - 1:
            nbImage = len(listImage) - 1
        if nbImage < 0:
            nbImage = 0

        print("len list", len(listImage))
        print(listImage[nbImage])
        name = listImage[nbImage].split("/")
        labelName.setText(name[-1])
        pixmap = QPixmap(listImage[nbImage])
        view.setPixmap(pixmap)
        view.setAlignment(QtCore.Qt.AlignCenter)
        view.show()

    def loadImage(self):
        self.checkDataLabelPath.setText(str(self.rootPath))
        self.channelList = self.parameters.index
        self.checkDataLabelChannel.setText(str(self.channelList))
        self.listImage = []
        for channel in self.channelList:
            self.rootPath = self.parameters["rootPath"][channel]
            image = str(self.rootPath + "/*" + "checkDataFit" + channel + ".png")
            self.listImage += glob.glob(image)
        print("list of images loaded")
        print(self.listImage)
        self.spinBoxSelectPlotCheckData.setMinimum(0)
        self.spinBoxSelectPlotCheckData.setMaximum(max(0, len(self.listImage) - 1))
        nbImage = self.spinBoxSelectPlotCheckData.value()
        if nbImage >= len(self.listImage) - 1:
            nbImage = len(self.listImage) - 1
        self.plot(self.checkDataLabel, self.listImage, nbImage, self.nameCheckImage)

    def displayFit(self):
        print(self.loadedFit)
        if self.loadedFit == False:
            for self.channel in self.parameters.index:
                self.rootPath = self.parameters["rootPath"][self.channel]
                folder = sorted(
                    [
                        join(join(self.rootPath, o), "analysis/")
                        for o in listdir(self.rootPath)
                        if (o[:3] == "201" or o[:3] == "EXP")
                        and os.path.isdir(os.path.join(self.rootPath, o))
                    ]
                )
                self.listPlotFitName = []
                for self.pathFolder in folder:
                    self.checkOutliersLabelPath_2.setText(self.pathFolder)

                    pathRlt = self.pathFolder + "resultIndiv/"
                    dirlist = [
                        filename
                        for filename in sorted(os.listdir(pathRlt))
                        if os.path.isdir(os.path.join(pathRlt, filename))
                    ]

                    if os.path.exists(self.pathFolder + "outliersbyEye.csv"):
                        now = str(datetime.datetime.now())[:19]
                        now = now.replace(":", "_")
                        shutil.move(
                            self.pathFolder + "outliersbyEye.csv",
                            self.pathFolder + "outliersbyEye_save" + now + ".csv",
                        )
                        print(
                            "copied outliersbyEye.csv as outliersbyEye_save"
                            + now
                            + ".csv"
                        )

                    names = []
                    for d in dirlist:
                        listImg = natsorted(listdir(pathRlt + "/" + d))
                        if ".DS_Store" in listImg:
                            listImg.pop(listImg.index(".DS_Store"))
                        listImg = [Im for Im in listImg if "FromDetectionTime" in Im]

                        names.extend([pathRlt + d + "/" + l for l in listImg])

                    fileOutliersAlgo = pd.DataFrame([])
                    if os.path.exists(
                        self.pathFolder + "outliersbyAlgo_" + self.channel + ".csv"
                    ):
                        fileOutliersAlgo = pd.read_csv(
                            self.pathFolder + "outliersbyAlgo_" + self.channel + ".csv",
                            header=None,
                        )

                    for imName in names:
                        drop = self.getDropFromName(imName)
                        if drop not in fileOutliersAlgo:
                            self.listPlotFitName = self.listPlotFitName + [imName]

                print("listPlotFitName", self.listPlotFitName)
                self.plot(
                    self.checkFitLabel, self.listPlotFitName, 0, self.nameCheckFit
                )
                self.loadedFit = True

        else:
            nbImage = self.spinBoxSelectCheckFit.value()
            self.plot(
                self.checkFitLabel, self.listPlotFitName, nbImage, self.nameCheckFit
            )

    def firstFit(self):
        print("*******************Start first fit*************************")
        self.channelList = self.parameters.index

        self.parameters["fromLag"] = False
        self.parameters["checkData"] = False
        self.parameters["deletionData"] = True

        i = 0
        for channel in self.channelList:
            self.rootPath = self.parameters["rootPath"][channel]
            folder = sorted(
                [
                    join(join(self.rootPath, o), "analysis/")
                    for o in listdir(self.rootPath)
                    if (o[:3] == "201" or o[:3] == "EXP")
                    and os.path.isdir(os.path.join(self.rootPath, o))
                ]
            )
            for path in folder:
                print(channel)
                i += np.ceil(100 / (len(folder) + 1))
                print(i)
                self.progressBarFit.setValue(i.astype(int))
                app.processEvents()
                wrapFitLoop(self.parameters, self.rootPath, path, channel)

        self.progressBarFit.setValue(100)

    def secondFit(self):
        print("*******************Start second fit*************************")
        i = 0
        self.channelList = self.parameters.index
        self.parameters["fromLag"] = True
        self.parameters["checkData"] = False
        self.parameters["deletionData"] = False

        for channel in self.channelList:
            self.rootPath = self.parameters["rootPath"][channel]
            folder = sorted(
                [
                    join(join(self.rootPath, o), "analysis/")
                    for o in listdir(self.rootPath)
                    if (o[:3] == "201" or o[:3] == "EXP")
                    and os.path.isdir(os.path.join(self.rootPath, o))
                ]
            )

            for path in folder:
                print(channel)
                i += np.ceil(100 / (len(folder) + 1))
                print(i)
                self.progressBarFit_2.setValue(i.astype(int))
                app.processEvents()
                wrapFitLoop(self.parameters, self.rootPath, path, channel)

        self.progressBarFit_2.setValue(100)


class DataFrameModel(QtCore.QAbstractTableModel):
    DtypeRole = QtCore.Qt.UserRole + 1000
    ValueRole = QtCore.Qt.UserRole + 1001

    def __init__(self, df=pd.DataFrame(), parent=None):
        super(DataFrameModel, self).__init__(parent)
        self._dataframe = df

    def setDataFrame(self, dataframe):
        self.beginResetModel()
        self._dataframe = dataframe.copy()
        self.endResetModel()

    def dataFrame(self):
        return self._dataframe

    dataFrame = QtCore.pyqtProperty(pd.DataFrame, fget=dataFrame, fset=setDataFrame)

    @QtCore.pyqtSlot(int, QtCore.Qt.Orientation, result=str)
    def headerData(
        self,
        section: int,
        orientation: QtCore.Qt.Orientation,
        role: int = QtCore.Qt.DisplayRole,
    ):
        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Horizontal:
                return self._dataframe.columns[section]
            else:
                return str(self._dataframe.index[section])
        return QtCore.QVariant()

    def rowCount(self, parent=QtCore.QModelIndex()):
        if parent.isValid():
            return 0
        return len(self._dataframe.index)

    def columnCount(self, parent=QtCore.QModelIndex()):
        if parent.isValid():
            return 0
        return self._dataframe.columns.size

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if not index.isValid() or not (
            0 <= index.row() < self.rowCount()
            and 0 <= index.column() < self.columnCount()
        ):
            return QtCore.QVariant()
        row = self._dataframe.index[index.row()]
        col = self._dataframe.columns[index.column()]
        dt = self._dataframe[col].dtype
        val = self._dataframe.loc[row, col]
        if role == QtCore.Qt.DisplayRole:
            return str(val)
        elif role == DataFrameModel.ValueRole:
            return val
        if role == DataFrameModel.DtypeRole:
            return dt
        return QtCore.QVariant()

    def roleNames(self):
        roles = {
            QtCore.Qt.DisplayRole: b"display",
            DataFrameModel.DtypeRole: b"dtype",
            DataFrameModel.ValueRole: b"value",
        }
        return roles


def main():
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication(sys.argv)
    else:
        app = QtWidgets.QApplication.instance()
    main = TabWidget()
    main.show()

    return app


if __name__ == "__main__":
    # app = main()
    app = QApplication(sys.argv)
    window = TabWidget()
    window.show()
    sys.exit(app.exec_())
