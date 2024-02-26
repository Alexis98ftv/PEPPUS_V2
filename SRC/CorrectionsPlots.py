#!/usr/bin/env python

########################################################################
# PEPPUS/SRC/CorrectionsPlots.py:
# This is the CorrectionsPlots Module of PEPPUS tool
#
#  Project:        PEPPUS
#  File:           CorrectionsPlots.py
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
from InputOutput import CorrIdx
sys.path.append(os.getcwd() + '/' + \
    os.path.dirname(sys.argv[0]) + '/' + 'COMMON')
from COMMON.Plots import generatePlot
import numpy as np
from collections import OrderedDict
from COMMON import GnssConstants as Const
from COMMON.Coordinates import xyz2llh

def initPlot(CorrFile, PlotConf, Title, Label):
    CorrFileName = os.path.basename(CorrFile)
    CorrFileNameSplit = CorrFileName.split('_')
    Rcvr = CorrFileNameSplit[1]
    DatepDat = CorrFileNameSplit[2]
    Date = DatepDat.split('.')[0]
    Year = Date[1:3]
    Doy = Date[4:]

    if "xLabel" not in PlotConf:
        PlotConf["xLabel"] = "Hour of Day %s" % Doy

    PlotConf["Title"] = "%s from %s on Year %s"\
        " DoY %s" % (Title, Rcvr, Year, Doy)

    PlotConf["Path"] = sys.argv[1] + '/OUT/PCOR/figures/%s/' % Label + \
        '%s_%s_Y%sD%s.png' % (Label, Rcvr, Year, Doy)

# Plot Confg
cfg = {
    "SatTrack"          : 1,
    "FlightTime"        : 1,
    "Dtr"               : 1,
    "STD"               : 1,
    "STD_elev"          : 1,
    "STropo_elev"       : 1,
    "SMp_elev"          : 1,
    "SNoiDiv"           : 1,
    "SAirb"             : 1,
    "SUERE"             : 1,
    "RCVR"              : 1,
    "C-RES"             : 1,
    "PH-RES"            : 1,
    "PH-RES_zoom"       : 1
}

def plotSatMonitoringTrack(CorrFile, CorrData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Satellite Tracks"
    PlotConf["xLabel"] = ""

    PlotConf["LonMin"] = -180
    PlotConf["LonMax"] = 180
    PlotConf["LatMin"] = -90
    PlotConf["LatMax"] = 90
    PlotConf["LonStep"] = 15
    PlotConf["LatStep"] = 10

    # PlotConf["yLabel"] = "Latitude [deg]"
    PlotConf["yTicks"] = range(PlotConf["LatMin"],PlotConf["LatMax"]+1,10)
    PlotConf["yLim"] = [PlotConf["LatMin"], PlotConf["LatMax"]]

    # PlotConf["xLabel"] = "Longitude [deg]"
    PlotConf["xTicks"] = range(PlotConf["LonMin"],PlotConf["LonMax"]+1,15)
    PlotConf["xLim"] = [PlotConf["LonMin"], PlotConf["LonMax"]]

    PlotConf["Grid"] = True

    PlotConf["Map"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    filter_flag = CorrData[CorrIdx["FLAG"]] == 1
    # Transform ECEF to Geodetic
    satX = CorrData[CorrIdx["SAT-X"]][filter_flag].to_numpy()
    satY = CorrData[CorrIdx["SAT-Y"]][filter_flag].to_numpy()
    satZ = CorrData[CorrIdx["SAT-Z"]][filter_flag].to_numpy()
    DataLen = len(satX)
    Longitude = np.zeros(DataLen)
    Latitude = np.zeros(DataLen)
    # transformer = Transformer.from_crs('epsg:4978', 'epsg:4326')
    for index in range(DataLen):
        x = satX[index]
        y = satY[index]
        z = satZ[index]
        Longitude[index], Latitude[index], h = xyz2llh(x, y, z)
        # Latitude[index], Longitude[index], h = transformer.transform(x, y, z)

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    Label = 0

    PlotConf["xData"][Label] = Longitude
    PlotConf["yData"][Label] = Latitude
    PlotConf["zData"][Label] = CorrData[CorrIdx["ELEV"]][filter_flag]

    # init plot
    Folder = "MONSAT_TRACKS"
    initPlot(CorrFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

def plotFlightTime(CorrFile, CorrData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Flight Time"

    PlotConf["yLabel"] = "Flight Time [miliseconds]"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '|'
    PlotConf["LineWidth"] = 1.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filter_flag = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"][Label] = CorrData[CorrIdx["SOD"]][filter_flag] / Const.S_IN_H
    PlotConf["yData"][Label] = CorrData[CorrIdx["FLIGHT-TIME"]][filter_flag]
    PlotConf["zData"][Label] = CorrData[CorrIdx["ELEV"]][filter_flag]

    # init plot
    Folder = "FLIGHT-TIME"
    initPlot(CorrFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

def plotDtr(CorrFile, CorrData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Relativistic Corrections (Dtr)"

    PlotConf["yLabel"] = "Relativistic Corrections (Dtr) [m]"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '|'
    PlotConf["LineWidth"] = 1.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filter_flag = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"][Label] = CorrData[CorrIdx["SOD"]][filter_flag] / Const.S_IN_H
    PlotConf["yData"][Label] = CorrData[CorrIdx["DTR"]][filter_flag]
    PlotConf["zData"][Label] = CorrData[CorrIdx["ELEV"]][filter_flag]

    # init plot
    Folder = "DTR"
    initPlot(CorrFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

def plotSlantTropoDelay(CorrFile, CorrData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Slant Tropo Delay"

    PlotConf["yLabel"] = "Slant Tropo Delay [m]"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '|'
    PlotConf["LineWidth"] = 1.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filter_flag = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"][Label] = CorrData[CorrIdx["SOD"]][filter_flag] / Const.S_IN_H
    PlotConf["yData"][Label] = CorrData[CorrIdx["STD"]][filter_flag]
    PlotConf["zData"][Label] = CorrData[CorrIdx["ELEV"]][filter_flag]

    # init plot
    Folder = "STD"
    initPlot(CorrFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

def plotSTD_Elev(CorrFile, CorrData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Slant Tropo Delay vs Elevation"

    PlotConf["yLabel"] = "Slant Tropo Delay vs Elevation [m]"

    PlotConf["xLabel"] = "Elevation [deg]"
    PlotConf["xLim"] = [0, 90]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '|'
    PlotConf["LineWidth"] = 1.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "PRN"
    PlotConf["ColorBarTicks"] = range(0,33)

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filter_flag = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][filter_flag]
    PlotConf["yData"][Label] = CorrData[CorrIdx["STD"]][filter_flag]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][filter_flag]

    # init plot
    Folder = "STDvsELEV"
    initPlot(CorrFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

def plotSTropo_Elev(CorrFile, CorrData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Sigma Tropo"

    PlotConf["yLabel"] = "Sigma Tropo [m]"

    PlotConf["xLabel"] = "Elevation [deg]"
    PlotConf["xLim"] = [0, 90]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '|'
    PlotConf["LineWidth"] = 1.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "PRN"
    PlotConf["ColorBarTicks"] = range(0,33)

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filter_flag = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][filter_flag]
    PlotConf["yData"][Label] = CorrData[CorrIdx["STROPO"]][filter_flag]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][filter_flag]

    # init plot
    Folder = "SigmaTROPOvsELEV"
    initPlot(CorrFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

def plotSMP_elev(CorrFile, CorrData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Sigma Multipath"

    PlotConf["yLabel"] = "Sigma Multipath [m]"

    PlotConf["xLabel"] = "Elevation [deg]"
    PlotConf["xLim"] = [0, 90]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '|'
    PlotConf["LineWidth"] = 1.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "PRN"
    PlotConf["ColorBarTicks"] = range(0,33)

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filter_flag = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][filter_flag]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SMP"]][filter_flag]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][filter_flag]

    # init plot
    Folder = "SigmaMPvsELEV"
    initPlot(CorrFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

def plotSNoiDiv_elev(CorrFile, CorrData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Sigma Noise + Divergence"

    PlotConf["yLabel"] = "Sigma Noise + Divergence [m]"

    PlotConf["xLabel"] = "Elevation [deg]"
    PlotConf["xLim"] = [0, 90]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '|'
    PlotConf["LineWidth"] = 1.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "PRN"
    PlotConf["ColorBarTicks"] = range(0,33)

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filter_flag = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][filter_flag]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SNOISEDIV"]][filter_flag]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][filter_flag]

    # init plot
    Folder = "SigmaNOISEDIVvsELEV"
    initPlot(CorrFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

def plotSAirb_elev(CorrFile, CorrData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Sigma Airborne"

    PlotConf["yLabel"] = "Sigma Airborne [m]"

    PlotConf["xLabel"] = "Elevation [deg]"
    PlotConf["xLim"] = [0, 90]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '|'
    PlotConf["LineWidth"] = 1.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "PRN"
    PlotConf["ColorBarTicks"] = range(0,33)

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filter_flag = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][filter_flag]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SAIR"]][filter_flag]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][filter_flag]

    # init plot
    Folder = "SigmaAIRvsELEV"
    initPlot(CorrFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

def plotSUERE_elev(CorrFile, CorrData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Sigma UERE"

    PlotConf["yLabel"] = "Sigma UERE [m]"

    PlotConf["xLabel"] = "Elevation [deg]"
    PlotConf["xLim"] = [0, 90]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '|'
    PlotConf["LineWidth"] = 1.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "PRN"
    PlotConf["ColorBarTicks"] = range(0,33)

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filter_flag = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][filter_flag]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SUERE"]][filter_flag]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][filter_flag]

    # init plot
    Folder = "SigmaUEREvsELEV"
    initPlot(CorrFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

def plotRCVR(CorrFile, CorrData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Receiver Clock Estimation"

    PlotConf["yLabel"] = "Receiver Clock Estimation [m]"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.5

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}

    Label = 0

    filter_flag = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"][Label] = CorrData[CorrIdx["SOD"]][filter_flag] / Const.S_IN_H
    PlotConf["yData"][Label] = CorrData[CorrIdx["RCVR-CLK"]][filter_flag]

    # init plot
    Folder = "RCVR-CLK"
    initPlot(CorrFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

def plotCodeResiduals(CorrFile, CorrData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Code Residuals"

    PlotConf["yLabel"] = "Code Residuals [m]"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "PRN"
    PlotConf["ColorBarTicks"] = range(0,33)

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filter_flag = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"][Label] = CorrData[CorrIdx["SOD"]][filter_flag] / Const.S_IN_H
    PlotConf["yData"][Label] = CorrData[CorrIdx["CODE-RES"]][filter_flag]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][filter_flag]

    # init plot
    Folder = "CODE-RES"
    initPlot(CorrFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

def plotPhaseResiduals(CorrFile, CorrData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Phase Residuals"

    PlotConf["yLabel"] = "Phase Residuals [m]"

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "PRN"
    PlotConf["ColorBarTicks"] = range(0,33)

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filter_flag = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"][Label] = CorrData[CorrIdx["SOD"]][filter_flag] / Const.S_IN_H
    PlotConf["yData"][Label] = CorrData[CorrIdx["PHASE-RES"]][filter_flag]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][filter_flag]

    # init plot
    Folder = "PHASE-RES"
    initPlot(CorrFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

def plotPhaseResidualsZoom(CorrFile, CorrData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    Title = "Phase Residuals with Zoom"

    PlotConf["yLabel"] = "Phase Residuals with Zoom[m]"
    PlotConf["yLim"] = [-8, 8]

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "PRN"
    PlotConf["ColorBarTicks"] = range(0,33)

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0

    filter_flag = CorrData[CorrIdx["FLAG"]] == 1

    PlotConf["xData"][Label] = CorrData[CorrIdx["SOD"]][filter_flag] / Const.S_IN_H
    PlotConf["yData"][Label] = CorrData[CorrIdx["PHASE-RES"]][filter_flag]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][filter_flag]

    # init plot
    Folder = "PHASE-RES-ZOOM"
    initPlot(CorrFile, PlotConf, Title, Folder)
    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Generate CorrPlots
def generateCorrPlots(CorrFile, Rcvr):
    # Purpose: generate output plots regarding Corrections results

    # Parameters
    # ==========
    # CorrFile: str
    #         Path to CORR PREPRO output file
    # Rcvr: str
    #           Receiver information

    # Returns
    # =======
    # Nothing
    # ----------------------------------------------------------
    # Satellite Tracks
    # ----------------------------------------------------------
    if (cfg["SatTrack"] == 1):
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["FLAG"], CorrIdx["ELEV"], \
                 CorrIdx["SAT-X"], CorrIdx["SAT-Y"], CorrIdx["SAT-Z"]])
    
        print( '\nPlot Satellite Monitoring Tracks ...')

        # Call Plot Function
        plotSatMonitoringTrack(CorrFile, CorrData)
    
    # Time of Flight
    # ----------------------------------------------------------
    if (cfg["FlightTime"] == 1):
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["SOD"], CorrIdx["ELEV"], CorrIdx["FLIGHT-TIME"], CorrIdx["FLAG"]])

        print( '\nPlot Time of Flight ...')

        # Call Plot Function
        plotFlightTime(CorrFile, CorrData)
    
    # Relativistic Corrections DTR
    # ----------------------------------------------------------
    if (cfg["Dtr"] == 1):
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["SOD"], CorrIdx["ELEV"], CorrIdx["DTR"], CorrIdx["FLAG"]])

        print( '\nPlot Relativistic Corrections (Dtr) ...')

        # Call Plot Function
        plotDtr(CorrFile, CorrData)
    
    # Slant Tropo Delay
    # ----------------------------------------------------------
    if (cfg["STD"] == 1):
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["SOD"], CorrIdx["ELEV"], CorrIdx["STD"], CorrIdx["FLAG"]])

        print( '\nPlot Slant Tropo Delay ...')

        # Call Plot Function
        plotSlantTropoDelay(CorrFile, CorrData)

    # STD VS Elevation
    # ----------------------------------------------------------
    if (cfg["STD_elev"] == 1):
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["ELEV"], CorrIdx["STD"], CorrIdx["PRN"], CorrIdx["FLAG"]])

        print( '\nPlot STD vs Elevation ...')

        # Call Plot Function
        plotSTD_Elev(CorrFile, CorrData)
    
    # SigmaTropo vs Elevation
    # ----------------------------------------------------------
    if (cfg["STropo_elev"] == 1):
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["ELEV"], CorrIdx["STROPO"], CorrIdx["PRN"], CorrIdx["FLAG"]])

        print( '\nPlot Sigma Tropo vs Elevation ...')

        # Call Plot Function
        plotSTropo_Elev(CorrFile, CorrData)

    # SigmaMultipath vs Elevation
    # ----------------------------------------------------------
    if (cfg["SMp_elev"] == 1):
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["ELEV"], CorrIdx["SMP"], CorrIdx["PRN"], CorrIdx["FLAG"]])

        print( '\nPlot Sigma Multipath vs Elevation ...')

        # Call Plot Function
        plotSMP_elev(CorrFile, CorrData)

    # Sigma Noise + Divergence
    # ----------------------------------------------------------
    if (cfg["SNoiDiv"] == 1):
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["ELEV"], CorrIdx["SNOISEDIV"], CorrIdx["PRN"], CorrIdx["FLAG"]])

        print( '\nPlot Sigma Noise + Divergence vs Elevation ...')

        # Call Plot Function
        plotSNoiDiv_elev(CorrFile, CorrData)
    
    # Sigma Airborne
    # ----------------------------------------------------------
    if (cfg["SAirb"] == 1):
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["ELEV"], CorrIdx["SAIR"], CorrIdx["PRN"], CorrIdx["FLAG"]])

        print( '\nPlot Sigma Airborne vs Elevation ...')

        # Call Plot Function
        plotSAirb_elev(CorrFile, CorrData)

    # Sigma UERE
    # ----------------------------------------------------------
    if (cfg["SUERE"] == 1):
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["ELEV"], CorrIdx["SUERE"], CorrIdx["PRN"], CorrIdx["FLAG"]])

        print( '\nPlot Sigma UERE vs Elevation ...')

        # Call Plot Function
        plotSUERE_elev(CorrFile, CorrData)
    
    # Receiver Clock Estimation
    # ----------------------------------------------------------
    if (cfg["RCVR"] == 1):
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["SOD"], CorrIdx["RCVR-CLK"], CorrIdx["FLAG"]])

        print( '\nPlot Receiver Clock Estimation ...')

        # Call Plot Function
        plotRCVR(CorrFile, CorrData)
    
    # Code Residuals
    # ----------------------------------------------------------
    if (cfg["C-RES"] == 1):
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["SOD"], CorrIdx["CODE-RES"], CorrIdx["PRN"], CorrIdx["FLAG"]])

        print( '\nPlot Code Residuals ...')

        # Call Plot Function
        plotCodeResiduals(CorrFile, CorrData)
    
    # Phase Residuals
    # ----------------------------------------------------------
    if (cfg["PH-RES"] == 1):
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["SOD"], CorrIdx["PHASE-RES"], CorrIdx["PRN"], CorrIdx["FLAG"]])

        print( '\nPlot Phase Residuals ...')

        # Call Plot Function
        plotPhaseResiduals(CorrFile, CorrData)
    
    # Phase Residuals Zoomed
    # ----------------------------------------------------------
    if (cfg["PH-RES_zoom"] == 1):
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["SOD"], CorrIdx["PHASE-RES"], CorrIdx["PRN"], CorrIdx["FLAG"]])

        print( '\nPlot Phase Residuals with Zoom ...')

        # Call Plot Function
        plotPhaseResidualsZoom(CorrFile, CorrData)

