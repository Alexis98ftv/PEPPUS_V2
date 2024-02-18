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



    return CorrInfo
