#!/usr/bin/env python

########################################################################
# PETRUS/SRC/Corrections.py:
# This is the Corrections Module of PEPPUS tool
#
#  Project:        PEPPUS
#  File:           Corrections.py
#  Date(YY/MM/DD): 16/02/21
#
#   Author: GNSS Academy
#   Copyright 2021 GNSS Academy
#
# -----------------------------------------------------------------
# Date       | Author             | Action
# -----------------------------------------------------------------
#
########################################################################

# Import External and Internal functions and Libraries
#----------------------------------------------------------------------
import sys, os
# Add path to find all modules
Common = os.path.dirname(os.path.dirname(
    os.path.abspath(sys.argv[0]))) + '/COMMON'
sys.path.insert(0, Common)
from collections import OrderedDict
from COMMON import GnssConstants as Const
from COMMON.Misc import findSun, crossProd
# from COMMON.Tropo import computeTropoMpp, computeZtd, computeSigmaTropo
from InputOutput import RcvrIdx, ObsIdx, SatPosIdx, SatClkIdx, SatApoIdx
import numpy as np
from bisect import bisect_left, bisect_right


def runCorrectMeas(Conf, Rcvr, ObsInfo, PreproObsInfo, 
SatPosInfo, SatClkInfo, SatApoInfo, SatComPos_1, Sod_1):

    # Purpose: correct GNSS preprocessed measurements and compute the first
    #          pseudo range residuals

    #          More in detail, this function handles the following:
    #          tasks:

    #             *  Compute the Satellite Antenna Phase Center position at the transmission time and corrected from the Sagnac
    #                effect interpolating the SP3 file positions
    #             *  Compute the Satellite Clock Bias interpolating the biases coming from the RINEX CLK file and
    #                applying the Relativistic Correction (DTR)
    #             *  Estimate the Slant Troposphere delay (STD) using MOPS model (ZTD) and its mapping function. 
    #             *  Correct the Pre-processed measurements from Geometrical Range, Satellite clock and Troposphere. 
    #             *  Build the Corrected Measurements and Measurement Residuals
    #             *  Build the Sigma UERE


    # Parameters
    # ==========
    # Conf: dict
    #         Configuration dictionary
    # Rcvr: list
    #         Receiver information: position, masking angle...
    # ObsInfo: list
    #         OBS info for current epoch
    #         ObsInfo[1][1] is the second field of the 
    #         second satellite
    # PreproObsInfo: dict
    #         Preprocessed observations for current epoch per sat
    #         PreproObsInfo["G01"]["C1"]
    # SatPosInfo: dict
    #         containing the SP3 file info
    # SatClkInfo: dict
    #         containing the RINEX CLK file info
    # SatApoInfo: dict
    #         containing the ANTEX file info
    # SatComPos_1: dict
    #         containing the previous satellite positions
    # Sod_1: dict
    #         containing the time stamp of previous satellite positions

    # Returns
    # =======
    # CorrInfo: dict
    #         Corrected measurements for current epoch per sat
    #         CorrInfo["G01"]["CorrectedPsr"]

    # Initialize output
    CorrInfo = OrderedDict({})

    # Initialize some values
    ResSum = 0.0
    ResN = 0

    # Get SoD
    Sod = int(float(ObsInfo[0][ObsIdx["SOD"]]))

    # Get DoY
    Doy = int(float(ObsInfo[0][ObsIdx["DOY"]]))

    # Get Year
    Year = int(float(ObsInfo[0][ObsIdx["YEAR"]]))

    # Find Sun position
    SunPos = findSun(Year, Doy, Sod)

    # Get receiver reference position
    RcvrRefPosXyz = np.array(\
                            (\
                                Rcvr[RcvrIdx["XYZ"]][0],
                                Rcvr[RcvrIdx["XYZ"]][1],
                                Rcvr[RcvrIdx["XYZ"]][2],
                            )
                        )

    # Loop over satellites
    for SatLabel, SatPrepro in PreproObsInfo.items():
        # Initialize output info
        SatCorrInfo = {
            "Sod": 0.0,             # Second of day
            "Doy": 0,               # Day of year
            "Elevation": 0.0,       # Elevation
            "Azimuth": 0.0,         # Azimuth
            "Flag": 1,              # 0: Not Used 1: Used for PA 2: Used for NPA
            "SatX": 0.0,            # X-Component of the Satellite CoP Position 
                                    # at transmission time and corrected from Sagnac
            "SatY": 0.0,            # Y-Component of the Satellite CoP Position  
                                    # at transmission time and corrected from Sagnac
            "SatZ": 0.0,            # Z-Component of the Satellite CoP Position  
                                    # at transmission time and corrected from Sagnac
            "ApoX": 0.0,            # X-Component of the Satellite APO in ECEF
            "ApoY": 0.0,            # Y-Component of the Satellite APO in ECEF
            "ApoZ": 0.0,            # Z-Component of the Satellite APO in ECEF
            "SatClk": 0.0,          # Satellite Clock Bias
            "FlightTime": 0.0,      # Signal Flight Time
            "Dtr": 0.0,             # Relativistic correction
            "Std": 0.0,             # Slant Tropospheric Delay
            "CorrCode": 0.0,        # Code corrected from delays
            "CorrPhase": 0.0,       # Phase corrected from delays
            "GeomRange": 0.0,       # Geometrical Range (distance between Satellite 
                                    # Position and Receiver Reference Position)
            "CodeResidual": 0.0,    # Code Residual
            "PhaseResidual": 0.0,   # Phase Residual
            "RcvrClk": 0.0,         # Receiver Clock estimation
            "SigmaTropo": 0.0,      # Sigma of the Tropo Delay Error
            "SigmaAirborne": 0.0,   # Sigma Airborne Error
            "SigmaNoiseDiv": 0.0,   # Sigma of the receiver noise + divergence
            "SigmaMultipath": 0.0,  # Sigma of the receiver multipath
            "SigmaUere": 0.0,       # Sigma User Equivalent Range Error (Sigma of 
                                    # the total residual error associated to the 
                                    # satellite)
            "TropoMpp": 0.0,        # Tropospheric mapping function

        } # End of SatCorrInfo

        # Prepare outputs
        # Get SoD
        SatCorrInfo["Sod"] = Sod
        # Get DoY
        SatCorrInfo["Doy"] = Doy
        # Get Elevation
        SatCorrInfo["Elevation"] = SatPrepro["Elevation"]
        # Get Azimuth
        SatCorrInfo["Azimuth"] = SatPrepro["Azimuth"]


        # Only for those Satellites with Status OK
        if SatPrepro["Status"] == 1:
            # Get Components of SatLabel to access InfoFiles
            Constel = SatLabel[0]
            Prn = int(SatLabel[1:])

            # Compute Satellite Clock Bias (linear interpolation between closer inputs) RELOJES GAP !!!
            #-----------------------------------------------------------------------
            clkBias = SatClkInfo[Constel][Prn][Sod]

            #SatClkBias = computeSatClkBias(Sod, SatClkInfo, SatLabel)
            #SatClkBias2 = lagrangeInterpolation(Sod, SatClkInfo, SatLabel, 2)
            SatClkPrn = SatClkInfo[Constel][Prn]
            SodList = np.array(list(SatClkPrn.keys()))
            SodData = np.array(list(SatClkPrn.values()))

            #if Sod == 300:
                #print("debug")

            point = bisect_right(SodList, Sod)
            
            # Compute Delta t
            #-----------------------------------------------------------------------
            DeltaT = SatPrepro["C1"] / Const.SPEED_OF_LIGHT

            # Compute Transmission Time
            TransmissionTime = Sod - DeltaT - clkBias

            # Compute Satellite CoM Position at Transmission Time
            # 10-point Lagrange interpolation between closer inputs (SP3 positions)
            #-----------------------------------------------------------------------
           
             

        else:
            SatCorrInfo["Flag"] == 0

        # prepare output
        CorrInfo[SatLabel] = SatCorrInfo
   
    return CorrInfo

# Linear interpolation to compute Satellite Clock Bias
def computeSatClkBias(Sod, SatClkInfo, SatLabel):
    t1 = 0
    t2 = 30

    if not (t1 <= Sod and t2 >= Sod):
        t1 = t1 + 30
        t2 = t2 + 30

    SatClkBias = SatClkInfo[SatLabel[0]][int(SatLabel[1:])][t1] + \
        (Sod-t1) * (SatClkInfo[SatLabel[0]][int(SatLabel[1:])][t2] - \
            SatClkInfo[SatLabel[0]][int(SatLabel[1:])][t1]) / (t2-t1)

    return SatClkBias

# Lagrange Interpolation
def lagrangeInterpolation(x, info, SatLabel, n):

    # Ordena la lista de acuerdo a la distancia entre cada punto y x
    lista_ordenada = sorted(info[SatLabel[0]][int(SatLabel[1:])].keys(), key=lambda punto: pointDistance(punto, (x, 0)))

    # Extrae las posiciones y valores ordenados
    posiciones_ordenadas = [punto[0] for punto in lista_ordenada]
    valores_ordenados = [punto[1] for punto in lista_ordenada]

     # Selecciona los primeros n elementos de las listas ordenadas
    x_interpolate_points = posiciones_ordenadas[:n]
    y_interpolate_points = valores_ordenados[:n]

    # Perform Lagrange Interpolation
    result = 0
    for i in range(len(x_interpolate_points)):
        term = y_interpolate_points[i]
        for j in range(len(x_interpolate_points)):
            if j != i:
                term = term * (x - x_interpolate_points[j]) / (x_interpolate_points[i] - x_interpolate_points[j])
        result += term

    return result

def pointDistance (point1, point2):
    return abs(point1 - point2)



