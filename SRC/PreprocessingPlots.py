#!/usr/bin/env python

########################################################################
# PETRUS/SRC/PreprocessingPlots.py:
# This is the PreprocessingPlots Module of PEPPUS tool
#
#  Project:        PEPPUS
#  File:           PreprocessingPlots.py
#  Date(YY/MM/DD): 05/02/21
#
#   Author: GNSS Academy
#   Copyright 2021 GNSS Academy
#
# -----------------------------------------------------------------
# Date       | Author             | Action
# -----------------------------------------------------------------
#
########################################################################

import sys, os
from pandas import unique
from pandas import read_csv
from InputOutput import PreproIdx
from InputOutput import REJECTION_CAUSE_DESC
sys.path.append(os.getcwd() + '/' + \
    os.path.dirname(sys.argv[0]) + '/' + 'COMMON')
from COMMON.Plots import generatePlot
import numpy as np
from collections import OrderedDict
from COMMON import GnssConstants as Const

def initPlot(PreproObsFile, PlotConf, Title, Label):
    PreproObsFileName = os.path.basename(PreproObsFile)
    PreproObsFileNameSplit = PreproObsFileName.split('_')
    Rcvr = PreproObsFileNameSplit[2]
    DatepDat = PreproObsFileNameSplit[3]
    Date = DatepDat.split('.')[0]
    Year = Date[1:3]
    Doy = Date[4:]

    PlotConf["xLabel"] = "Hour of Day %s" % Doy

    PlotConf["Title"] = "%s from %s on Year %s"\
        " DoY %s" % (Title, Rcvr, Year, Doy)

    PlotConf["Path"] = sys.argv[1] + '/OUT/PPVE/figures/%s/' % Label + \
        '%s_%s_Y%sD%s.png' % (Label, Rcvr, Year, Doy)

# Plot Confg
cfg = {
    "SatVisibility"         : 1,
    "NumSatellites"         : 1,
    "PolarView"             : 1,
    "RejFlags"              : 1,
    "CodeRate"              : 1,
    "CodeRateStep"          : 1,
    "PhaseRate"             : 1,
    "PhaseRateStep"         : 1,
    "VTEC"                  : 1,
    "AATR"                  : 1
}

# Plot Satellite Visibility
def plotSatVisibility(PreproObsFile, PreproObsData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Satellite Visibility"

    PlotConf["yLabel"] = "GPS-PRN"
    PlotConf["yTicks"] = sorted(unique(PreproObsData[PreproIdx["PRN"]]))
    PlotConf["yTicksLabels"] = sorted(unique(PreproObsData[PreproIdx["PRN"]]))
    PlotConf["yLim"] = [0, max(unique(PreproObsData[PreproIdx["PRN"]])) + 1]

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '|'
    PlotConf["LineWidth"] = 5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    PlotConf["Color"] = {}

    Label = "status_1"

    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]] / Const.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["PRN"]]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]]

    Label = "status_0"
    
    FilterStatus0 = PreproObsData[PreproIdx["STATUS"]] == 0
    PlotConf["Color"][Label] = "gray"

    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][FilterStatus0] / Const.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["PRN"]][FilterStatus0]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]][FilterStatus0]

    # init plot
    Folder = "SAT_VISIBILITY"
    initPlot(PreproObsFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Number of Satellites
def plotNumSats(PreproObsFile, PreproObsData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Number of Satellites"

    PlotConf["yLabel"] = "Number of Satellites"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '-'
    PlotConf["LineWidth"] = 0.25

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}

    PlotConf["Color"] = {}
    PlotConf["Legend"] = True
    
    Label = ["Raw", "OK"]
    colors = ["orange", "green"]

    sod_list = np.array(sorted(unique(PreproObsData[PreproIdx["SOD"]])), dtype=float)
    raw_list = np.zeros(len(sod_list))
    ok_list = np.zeros(len(sod_list))

    for i, sod in enumerate(sod_list):
        filterSod = PreproObsData[PreproIdx["SOD"]] == sod
        raw = len(PreproObsData[PreproIdx["STATUS"]][filterSod])
        ok = sum(PreproObsData[PreproIdx["STATUS"]][filterSod])
        raw_list[i] = raw
        ok_list[i] = ok

    data = [raw_list, ok_list]
    for idx, label in enumerate(Label):
        PlotConf["Color"][label] = colors[idx]
        PlotConf["xData"][label] = sod_list / Const.S_IN_H
        PlotConf["yData"][label] = data[idx]

    # init plot
    Folder = "NSATS"
    initPlot(PreproObsFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Satellite Polar View - CHALLENGE
def plotSatPolarView(PreproObsFile, PreproObsData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Satellite Polar View"

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarTicks"] = range(0,33)

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    PlotConf["Polar"] = True

    azim = PreproObsData[PreproIdx["AZIM"]]
    elev = PreproObsData[PreproIdx["ELEV"]]

    theta_azim_rad = np.radians(azim)
    r_elev = 90 - elev

    filter_cond = PreproObsData[PreproIdx["STATUS"]] != 0
    Label = 0

    PlotConf["xData"][Label] = theta_azim_rad[filter_cond]
    PlotConf["yData"][Label] = r_elev[filter_cond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["PRN"]][filter_cond]

    # init plot
    Folder = "POLAR_VIEW"
    initPlot(PreproObsFile, PlotConf, Title, Folder)
    # Reset xLabel from initPlot
    PlotConf["xLabel"] = ""
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Rejection Flags
def plotRejectionFlags(PreproObsFile, PreproObsData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Rejection Flags"

    PlotConf["yLabel"] = "Rejection Flags"
    
    PlotConf["yTicks"] = list(REJECTION_CAUSE_DESC.values())
    PlotConf["yTicksLabels"] = list(REJECTION_CAUSE_DESC.keys())
    PlotConf["yLim"] = [0, max(REJECTION_CAUSE_DESC.values()) + 1]
        
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 5

    PlotConf["ColorBar"] = "viridis"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarTicks"] = range(0,33)

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filter_cond = PreproObsData[PreproIdx["REJECT"]] != 0

    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][filter_cond] / Const.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["REJECT"]][filter_cond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["PRN"]][filter_cond]

    # init plot
    Folder = "REJ_FLAGS"
    initPlot(PreproObsFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Code Rate
def plotCodeRate(PreproObsFile, PreproObsData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Code Rate"

    PlotConf["yLabel"] = "Code Rate [m/s]"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filterCond = PreproObsData[PreproIdx["STATUS"]] == 1

    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][filterCond] / Const.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["CODE RATE"]][filterCond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]][filterCond]

    # init plot
    Folder = "CODE_RATE"
    initPlot(PreproObsFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Code Rate Step
def plotCodeRateStep(PreproObsFile, PreproObsData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Code Rate Step"

    PlotConf["yLabel"] = "Code Rate Step [m/s^2]"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filterCond = PreproObsData[PreproIdx["STATUS"]] == 1

    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][filterCond] / Const.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["CODE ACC"]][filterCond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]][filterCond]

    # init plot
    Folder = "CODE_RATE_STEP"
    initPlot(PreproObsFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Phase Rate
def plotPhaseRate(PreproObsFile, PreproObsData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Phase Rate"

    PlotConf["yLabel"] = "Phase Rate [m/s]"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filterCond = PreproObsData[PreproIdx["STATUS"]] == 1

    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][filterCond] / Const.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["PHASE RATE"]][filterCond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]][filterCond]

    # init plot
    Folder = "PHASE_RATE"
    initPlot(PreproObsFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Phase Rate Step
def plotPhaseRateStep(PreproObsFile, PreproObsData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Phase Rate Step"

    PlotConf["yLabel"] = "Phase Rate Step [m/s^2]"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filterCond = PreproObsData[PreproIdx["STATUS"]] == 1

    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][filterCond] / Const.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["PHASE ACC"]][filterCond]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]][filterCond]

    # init plot
    Folder = "PHASE_RATE_STEP"
    initPlot(PreproObsFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# VTEC Gradient
def plotVtecGradient(PreproObsFile, PreproObsData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "VTEC Gradient"

    PlotConf["yLabel"] = "VTEC Gradient [mm/s]"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filterCond = PreproObsData[PreproIdx["STATUS"]] == 1
    filterCond2 = PreproObsData[PreproIdx["VTEC RATE"]] != Const.NAN

    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][filterCond][filterCond2] / Const.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["VTEC RATE"]][filterCond][filterCond2]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]][filterCond][filterCond2]

    # init plot
    Folder = "VTEC"
    initPlot(PreproObsFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# AATR index
def plotAatr(PreproObsFile, PreproObsData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "AATR"

    PlotConf["yLabel"] = "AATR [mm/s]"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filterCond = PreproObsData[PreproIdx["STATUS"]] == 1
    filterCond2 = PreproObsData[PreproIdx["iAATR"]] != Const.NAN

    PlotConf["xData"][Label] = PreproObsData[PreproIdx["SOD"]][filterCond][filterCond2] / Const.S_IN_H
    PlotConf["yData"][Label] = PreproObsData[PreproIdx["iAATR"]][filterCond][filterCond2]
    PlotConf["zData"][Label] = PreproObsData[PreproIdx["ELEV"]][filterCond][filterCond2]

    # init plot
    Folder = "AATR"
    initPlot(PreproObsFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Generate PreProPlots
def generatePreproPlots(PreproObsFile):
    # Purpose: generate output plots regarding Preprocessing results

    # Parameters
    # ==========
    # PreproObsFile: str
    #         Path to PREPRO OBS output file

    # Returns
    # =======
    # Nothing
    # ----------------------------------------------------------
    # Satellite Visibility
    # ----------------------------------------------------------
    if (cfg["SatVisibility"] == 1):
        # Read the cols we need from PREPRO OBS file
        PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[PreproIdx["SOD"], PreproIdx["CONST"], PreproIdx["PRN"], PreproIdx["ELEV"], PreproIdx["STATUS"]])

        print( '\nPlot Satellite Visibility Periods ...')

        # Call Plot Function
        plotSatVisibility(PreproObsFile, PreproObsData)
      
    # Number of Satellite
    # ----------------------------------------------------------
    if (cfg["NumSatellites"] == 1):
        # Read the cols we need from PREPRO OBS file
        PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[PreproIdx["SOD"], PreproIdx["STATUS"]])

        print( '\nPlot Number of Satellites ...')

        # Call Plot Function
        plotNumSats(PreproObsFile, PreproObsData)

    # Polar View
    # ----------------------------------------------------------
    if (cfg["PolarView"] == 1):
        # Read the cols we need from PREPRO OBS file
        PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[PreproIdx["ELEV"], PreproIdx["AZIM"], PreproIdx["PRN"], PreproIdx["STATUS"]])

        print( '\nPlot Satellite Polar View ...')

        # Call Plot Function
        plotSatPolarView(PreproObsFile, PreproObsData)

    # Rejection Flags
    # ----------------------------------------------------------
    if (cfg["RejFlags"] == 1):
        # Read the cols we need from PREPRO OBS file
        PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[PreproIdx["SOD"], PreproIdx["PRN"], PreproIdx["REJECT"]])

        print( '\nPlot Rejection Flags ...')

        # Call Plot Function
        plotRejectionFlags(PreproObsFile, PreproObsData)

    # Code Rate
    # ----------------------------------------------------------
    if (cfg["CodeRate"] == 1):
        # Read the cols we need from PREPRO OBS file
        PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[PreproIdx["SOD"], PreproIdx["ELEV"], PreproIdx["STATUS"], PreproIdx["CODE RATE"]])

        print( '\nPlot Code Rate for OK measurements...')

        # Call Plot Function
        plotCodeRate(PreproObsFile, PreproObsData)
    
    # Code Rate Step
    # ----------------------------------------------------------
    if (cfg["CodeRateStep"] == 1):
        # Read the cols we need from PREPRO OBS file
        PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[PreproIdx["SOD"], PreproIdx["ELEV"], PreproIdx["STATUS"], PreproIdx["CODE ACC"]])

        print( '\nPlot Code Rate Step for OK measurements...')

        # Call Plot Function
        plotCodeRateStep(PreproObsFile, PreproObsData)
    
    # Phase Rate
    # ----------------------------------------------------------
    if (cfg["PhaseRate"] == 1):
        # Read the cols we need from PREPRO OBS file
        PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[PreproIdx["SOD"], PreproIdx["ELEV"], PreproIdx["STATUS"], PreproIdx["PHASE RATE"]])

        print( '\nPlot Phase Rate for OK measurements...')

        # Call Plot Function
        plotPhaseRate(PreproObsFile, PreproObsData)
    
    # Phase Rate Step
    # ----------------------------------------------------------
    if (cfg["PhaseRateStep"] == 1):
        # Read the cols we need from PREPRO OBS file
        PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[PreproIdx["SOD"], PreproIdx["ELEV"], PreproIdx["STATUS"], PreproIdx["PHASE ACC"]])

        print( '\nPlot Phase Rate Step for OK measurements...')

        # Call Plot Function
        plotPhaseRateStep(PreproObsFile, PreproObsData)
    
    # VTEC Gradient
    # ----------------------------------------------------------
    if (cfg["VTEC"] == 1):
        # Read the cols we need from PREPRO OBS file
        PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[PreproIdx["SOD"], PreproIdx["ELEV"], PreproIdx["STATUS"], PreproIdx["VTEC RATE"]])

        print( '\nPlot VTEC gradients for OK measurements...')

        # Call Plot Function
        plotVtecGradient(PreproObsFile, PreproObsData)
    
    # Instantaneous AATR
    # ----------------------------------------------------------
    if (cfg["AATR"] == 1):
        # Read the cols we need from PREPRO OBS file
        PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[PreproIdx["SOD"], PreproIdx["ELEV"], PreproIdx["STATUS"], PreproIdx["iAATR"]])

        print( '\nPlot Instantaneous AATR for OK measurements...')

        # Call Plot Function
        plotAatr(PreproObsFile, PreproObsData)


