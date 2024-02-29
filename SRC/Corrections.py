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
from COMMON.Tropo import computeGeoidHeight
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
    # Variables associated with
    # First estimation of RcvrClock
    CodeResiduals = []
    WeightedUEREs = []


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
            # Caring of gaps
            try:
                if not Sod_1[SatLabel]:
                    SatCorrInfo["Flag"] = 0
                else:
                    if Sod - Sod_1[SatLabel] > 30:
                        SatCorrInfo["Flag"] = 0
            except KeyError:
                Sod_1[SatLabel] = {}
                SatComPos_1[SatLabel] = {}
                SatCorrInfo["Flag"] = 0
            
            # Compute Satellite Clock Bias (linear interpolation between closer inputs)
            #-----------------------------------------------------------------------
            clkBias, gap = computeSatClkBias(Sod, SatClkInfo, SatLabel)
            # SatClkBias in meters
            SatClkBias = clkBias*Const.SPEED_OF_LIGHT

            if gap == True:
                # Can't use these measurements
                SatCorrInfo["Flag"] = 0
                SatCorrInfo["SatClk"] = SatClkBias
            else:

                # Compute Delta t
                #-----------------------------------------------------------------------
                DeltaT = SatPrepro["C1"] / Const.SPEED_OF_LIGHT

                # Compute Transmission Time
                TransmissionTime = Sod - DeltaT - clkBias

                # Compute Satellite CoM Position at Transmission Time
                # 10-point Lagrange interpolation between closer inputs (SP3 positions)
                #-----------------------------------------------------------------------
                SatComPos = computeSatComPos(TransmissionTime, SatPosInfo, SatLabel)
                
                # Compute Flight Time
                #-----------------------------------------------------------------------
                FlightTime = computeFlightTime(SatComPos, RcvrRefPosXyz)

                # Apply Sagnac correction
                #-----------------------------------------------------------------------
                SatComPos = applySagnac(SatComPos, FlightTime)
                
                # Compute APO in ECEF from ANTEX APOs in
                # satellite-body reference frame
                #-----------------------------------------------------------------------
                Apofree = computeSatApo(SatComPos, SunPos, SatApoInfo, SatLabel)

                # Apply APOs to the Satellite Position
                SatCopPos = SatComPos + Apofree

                # Compute Dtr (Relativistic correction)
                #-----------------------------------------------------------------------
                Dtr = computeDtr(SatComPos_1, SatComPos, Sod, Sod_1, SatLabel)

                # If no prev position and time reject measurement
                if Dtr == Const.NAN or SatCorrInfo["Flag"] == 0:
                    SatCorrInfo["Flag"] = 0
                    #-----------------------------------
                    # Prepare output
                    #-----------------------------------
                    SatCorrInfo["SatX"] = SatCopPos[0]
                    SatCorrInfo["SatY"] = SatCopPos[1]
                    SatCorrInfo["SatZ"] = SatCopPos[2]
                    SatCorrInfo["ApoX"] = Apofree[0]
                    SatCorrInfo["ApoY"] = Apofree[1]
                    SatCorrInfo["ApoZ"] = Apofree[2]
                    SatCorrInfo["SatClk"] = SatClkBias
                    SatCorrInfo["FlightTime"] = FlightTime
                    SatCorrInfo["Dtr"] = Const.NAN
                    #-----------------------------------
                    # Store values for next epoch
                    #-----------------------------------
                    Sod_1[SatLabel] = Sod
                    SatComPos_1[SatLabel] = SatComPos
                # else: Continue
                else: 

                    # Apply Dtr to Sat Clock Bias
                    SatClkBias = SatClkBias + Dtr

                    # Compute the STD: Slant Tropo Delay and associated SigmaTROPO
                    # Refer to MOPS guidelines in Appendix A section A.4.2.4
                    #-----------------------------------------------------------------------
                    # Compute Tropospheric Mapping Function
                    TropoMpp = computeTropoMpp(SatCorrInfo["Elevation"])

                    # Compute the Slant Tropospheric Delay Error Sigma (Tropospheric Vertical Error = 0.12 meters)
                    SigmaTROPO = 0.12 * TropoMpp

                    # Compute the Slant Tropospheric Delay
                    STD = computeSlantTropoDelay(TropoMpp, Rcvr, Doy)

                    # Compute User Airborne Sigma. Ref: MOPS-DO-229D Section J.2.4
                    #-----------------------------------------------------------------------
                    # Consider Maximum Signal Level when satellite elevation is greater
                    # than Conf.ELEV_NOISE_TH=20, and Minimum Signal Level otherwise
                    # Apply Conf.SIGMA_AIR_DF factor to both the MP and Noise components

                    SigmaAir, Sigma_MP, Sigma_noise_divg = computeSigmaAir(SatCorrInfo["Elevation"], Conf)

                    # Compute Sigma UERE by combining all Sigma contributions
                    #-----------------------------------------------------------------------
                    SigmaUERE = computeSigmaUERE(Conf, SigmaTROPO, SigmaAir)

                    # Corrected Measurements from previous information
                    #-----------------------------------------------------------------------
                    CorrCode = SatPrepro["IF_C"] + SatClkBias - STD
                    CorrPhase = SatPrepro["IF_L"] + SatClkBias - STD
                    
                    # Compute the Geometrical Range
                    #-----------------------------------------------------------------------
                    GeomRange = computeGeomRange(SatCopPos, RcvrRefPosXyz)

                    # Compute the first Residual removing the geometrical range
                    # They include Receiver Clock estimation
                    #-----------------------------------------------------------------------
                    CodeResidual = CorrCode - GeomRange
                    PhaseResidual = CorrPhase - GeomRange

                    # Prepare First estimation of RcvrClock
                    #-----------------------------------------------------------------------
                    CodeResiduals.append(CodeResidual)
                    WeightedUEREs.append(1/SigmaUERE**2)

                    #-----------------------------------
                    # Prepare output
                    #-----------------------------------
                    SatCorrInfo["SatX"] = SatCopPos[0]
                    SatCorrInfo["SatY"] = SatCopPos[1]
                    SatCorrInfo["SatZ"] = SatCopPos[2]
                    SatCorrInfo["ApoX"] = Apofree[0]
                    SatCorrInfo["ApoY"] = Apofree[1]
                    SatCorrInfo["ApoZ"] = Apofree[2]
                    SatCorrInfo["SatClk"] = SatClkBias
                    SatCorrInfo["FlightTime"] = FlightTime
                    SatCorrInfo["Dtr"] = Dtr
                    SatCorrInfo["Std"] = STD
                    SatCorrInfo["CorrCode"] = CorrCode
                    SatCorrInfo["CorrPhase"] = CorrPhase
                    SatCorrInfo["GeomRange"] = GeomRange
                    SatCorrInfo["CodeResidual"] = CodeResidual
                    SatCorrInfo["PhaseResidual"] = PhaseResidual

                    SatCorrInfo["SigmaTropo"] = SigmaTROPO
                    SatCorrInfo["SigmaAirborne"] = SigmaAir
                    SatCorrInfo["SigmaNoiseDiv"] = Sigma_noise_divg
                    SatCorrInfo["SigmaMultipath"] = Sigma_MP
                    SatCorrInfo["SigmaUere"] = SigmaUERE
                    SatCorrInfo["TropoMpp"] = TropoMpp
                    #-----------------------------------
                    # Store values for next epoch
                    #-----------------------------------
                    Sod_1[SatLabel] = Sod
                    SatComPos_1[SatLabel] = SatComPos
                    
        else:
            SatCorrInfo["Flag"] = 0
            SatComPos_1[SatLabel] = {}
            Sod_1[SatLabel] = {}

        # prepare output
        CorrInfo[SatLabel] = SatCorrInfo
    
    # Estimate the Receiver Clock first guess as a weighted average of the Residuals
    # (with the weights W=1/UERE2)(CodeResidual, SigmaUERE)          
    RcvrClock = estimateRcvrClk(CodeResiduals, WeightedUEREs)
    # Loop over satellites
    for SatLabel, SatCorr in CorrInfo.items():
        if SatCorr["Flag"] != 0:
            SatCorr["RcvrClk"] = RcvrClock
            SatCorr["CodeResidual"] = SatCorr["CodeResidual"]-RcvrClock
            SatCorr["PhaseResidual"] = SatCorr["PhaseResidual"]-RcvrClock
    
        
    return CorrInfo


#######################################################
# EXTERNAL FUNCTIONS
#######################################################

# Lagrange Interpolation
def lagrangeInterpolation(x, y, t, n):
    # Perform Lagrange Interpolation
    result = 0.0
    for i in range(n):
        term = y[i]
        for j in range(n):
            if j != i:
                term = term * (t - x[j]) / (x[i] - x[j])
        result += term
    return result

# Compute Satellite Clock Bias 
def computeSatClkBias(Sod, SatClkInfo, SatLabel):
    # Get Components of SatLabel to access InfoFiles
    Constel = SatLabel[0]
    Prn = int(SatLabel[1:])
    SatClkPrn = SatClkInfo[Constel][Prn]
    SodList = np.array(list(SatClkPrn.keys()))

    # Care of gaps in SatClkInfo
    gap = False
    
    # No need of interpolation if it's in the infoFile
    if Sod in SodList:
        clkBias = SatClkPrn[Sod]
    else:
        n = 2
        # Obtain the Sod of near points
        position = bisect_right(SodList, Sod)
        x1 = SodList[position-1]
        x2 = SodList[position]
        x = [x1, x2]
        # Take into account the possible gap in SatClkInfo
        if ((abs(x1-Sod) > 300) or (abs(Sod-x2) > 300)):
            gap = True
            clkBias = Const.NAN
            return clkBias, gap

        # Obtain the ClkBiases of near points to interpolate
        y1 = SatClkPrn[x1]
        y2 = SatClkPrn[x2]
        y = [y1, y2]

        # Lagrange Interpolation (n=2: Linear Interpolation)
        clkBias = lagrangeInterpolation(x, y, Sod, n)

    return clkBias, gap

# Compute Satellite CoM Position at Transmission Time with Lagrange Interpolation
def computeSatComPos(TransmissionTime, SatPosInfo, SatLabel):
    # Get Components of SatLabel to access InfoFiles
    Constel = SatLabel[0]
    Prn = int(SatLabel[1:])
    SatPosPrn = SatPosInfo[Constel][Prn]

    SodList = np.array(list(SatPosPrn.keys()))
    tt = TransmissionTime
    # No need of interpolation if it's in the infoFile
    if tt in SodList:
        SatComPos = SatPosPrn[tt]  
    else:
        #SodData = np.array(list(SatPosPrn.values())) 
        n = 10
        sod_inter = [0.0] * n

        # Fill the list of points to interpolate
        if tt < Const.S_IN_H:
            for i in range(n):
                sod_inter[i] = SodList[i]
            
        elif tt > 81900.0:
            for i in range(n):
                sod_inter[-i-1] = SodList[-i-1]
            
        else:
            # Obtain positions of near points 
            pos5 = bisect_right(SodList, tt)
            pos4 = pos5 - 1
            # Add to the array the nearest 10 points
            for i in range(5): # (n/2)
                sod_inter[5+i] = SodList[pos5 + i]
                sod_inter[4-i] = SodList[pos4 - i]

        # Values to interpolate
        xCM_inter = [0.0] * n
        yCM_inter = [0.0] * n
        zCM_inter = [0.0] * n

        # Obtain values for each position of points
        for i in range(n):
            xCM_inter[i] = SatPosPrn[sod_inter[i]][0]
            yCM_inter[i] = SatPosPrn[sod_inter[i]][1]
            zCM_inter[i] = SatPosPrn[sod_inter[i]][2]
        
        # Interpolate for each coordinate (X,Y,Z)
        xCM_interpolated = lagrangeInterpolation(sod_inter, xCM_inter, tt, n)
        yCM_interpolated = lagrangeInterpolation(sod_inter, yCM_inter, tt, n)
        zCM_interpolated = lagrangeInterpolation(sod_inter, zCM_inter, tt, n)

        SatComPos = [xCM_interpolated, yCM_interpolated, zCM_interpolated]

    return SatComPos

# Compute Flight Time
def computeFlightTime(SatComPos, RcvrRefPosXyz):

    # Calculate distance
    # Coordinates in meters
    distance = np.sqrt((SatComPos[0]*Const.M_IN_KM - RcvrRefPosXyz[0])**2 + \
                       (SatComPos[1]*Const.M_IN_KM - RcvrRefPosXyz[1])**2 + \
                       (SatComPos[2]*Const.M_IN_KM - RcvrRefPosXyz[2])**2)
   
    # FlightTime in ms
    FlightTime = (distance / Const.SPEED_OF_LIGHT)*Const.MS_IN_S

    return FlightTime

# Apply Sagnac Corrections to the satellite
def applySagnac(SatComPos, FlightTime):
    # Angle that satellite should rotate
    theta = Const.OMEGA_EARTH*FlightTime/Const.MS_IN_S

    RotationMatrix = [
    [np.cos(theta), np.sin(theta), 0],
    [-np.sin(theta), np.cos(theta), 0],
    [0, 0, 1]
    ]

    RotatedSatPos = np.dot(np.array(RotationMatrix), np.array(SatComPos)*1000)

    return RotatedSatPos

# Compute APO in ECEF from ANTEX APOs in satellite-body reference frame
def computeSatApo(SatComPos, SunPos, SatApoInfo, SatLabel):
    # Get Components of SatLabel to access InfoFiles
    Constel = SatLabel[0]
    Prn = int(SatLabel[1:])
    FreqL1 = "L1"
    FreqL2 = "L2"
    ApoInfoL1 = SatApoInfo[Constel][Prn][FreqL1]/Const.MM_IN_M
    ApoInfoL2 = SatApoInfo[Constel][Prn][FreqL2]/Const.MM_IN_M

    e = (SunPos - SatComPos) / np.linalg.norm(SunPos - SatComPos) 
    k = -(SatComPos / np.linalg.norm(SatComPos))

    j = crossProd(k, e)
    i = crossProd(j, k)

    # Unitary vectors
    j_mag = np.sqrt(j[0]**2+j[1]**2+j[2]**2)
    i_mag = np.sqrt(i[0]**2+i[1]**2+i[2]**2)
    j_u = j*(1/j_mag)
    i_u = i*(1/i_mag)
    
    R1 = [i_u,j_u,k]
    R = np.transpose(R1)

    ApoL1 = np.dot(R, ApoInfoL1)
    ApoL2 = np.dot(R, ApoInfoL2)
    # Ionofree comb of Antenna Phase Offsets of L1 and L2
    Apofree = (ApoL2 - Const.GPS_GAMMA_L1L2*ApoL1) / (1-Const.GPS_GAMMA_L1L2)

    return Apofree

# Compute Dtr (Relativistic correction)
def computeDtr(SatComPos_1, SatComPos, Sod, Sod_1, SatLabel):
    # Check last values
    try:
        if not Sod_1[SatLabel]:
            Dtr = Const.NAN
        else:
            velocity = (SatComPos-SatComPos_1[SatLabel])/(Sod-Sod_1[SatLabel])
            Dtr = -2*(np.dot(SatComPos, velocity)/Const.SPEED_OF_LIGHT)
    except KeyError:
        Sod_1[SatLabel] = {}
        SatComPos_1[SatLabel] = {}
        Dtr = Const.NAN

    return Dtr

# Compute Tropospheric Mapping Function
def computeTropoMpp(elev_deg):
    elev_rad = elev_deg * (np.pi/180)
    mpp = 1.001 / np.sqrt(0.002001 + np.sin(elev_rad)**2)

    return mpp

# Compute the Slant Tropospheric Delay: STD = (d_hyd + d_wet) * TropoMpp
def computeSlantTropoDelay(TropoMpp, Rcvr, Doy):      
    # METEOROLOGICAL PARAMETERS FOR TROPOSPHERIC DELAYS
    # Dependency on Latitude(º) [15 or less, 30, 45, 60, 75 or greater]
    #-----------------------------------------------------------------------
    # Average
    #---------------------
    # Pressure [mbar]
    P_0 = [1013.25, 1017.25, 1015.75, 1011.75, 1013.00]
    # Temperature [K]
    T_0 = [299.65, 294.15, 283.15, 272.15, 263.65]
    # Water vapor pressure [mbar]
    e_0 = [26.31, 21.79, 11.66, 6.78, 4.11]
    # Temperrature lapse rate [K/m]
    Beta_0 = [6.30e-3, 6.05e-3, 5.58e-3, 5.39e-3, 4.53e-3] 
    # Water vapor "lapse rate" [dimensionless]
    Lambda_0 = [2.77, 3.15, 2.57, 1.81, 1.55]
    
    # Seasonal Variation
    #---------------------
    delta_P_0 = [0.00, -3.75, -2.25, -1.75, -0.50]
    delta_T_0 = [0.00, 7.00, 11.00, 15.00, 14.50]
    delta_e_0 = [0.00, 8.85, 7.24, 5.36, 3.39]
    delta_B_0 = [0.00e-3, 0.25e-3, 0.32e-3, 0.81e-3, 0.62e-3]
    delta_L_0 = [0.00, 0.33, 0.46, 0.74, 0.30]

    # Latitudes
    #---------------------
    latitudes = [15.00, 30.00, 45.00, 60.00, 75.00]

    # Latitude, Longitud and Altitude of Receiver
    LatRx = Rcvr[RcvrIdx["LAT"]]
    LonRx = Rcvr[RcvrIdx["LON"]]
    AltRx = Rcvr[RcvrIdx["ALT"]]

    # No need to interpolation if lat < 15º
    if LatRx <= latitudes[0]:
        P_interp = P_0[0]
        T_interp = T_0[0]
        e_interp = e_0[0]
        B_interp = Beta_0[0]
        L_interp = Lambda_0[0]

        delta_P_interp = delta_P_0[0]
        delta_T_interp = delta_T_0[0]
        delta_e_interp = delta_e_0[0]
        delta_B_interp = delta_B_0[0]
        delta_L_interp = delta_L_0[0]

    # No need to interpolation if lat > 75º
    elif LatRx >= latitudes[4]:
        P_interp = P_0[4]
        T_interp = T_0[4]
        e_interp = e_0[4]
        B_interp = Beta_0[4]
        L_interp = Lambda_0[4]

        delta_P_interp = delta_P_0[4]
        delta_T_interp = delta_T_0[4]
        delta_e_interp = delta_e_0[4]
        delta_B_interp = delta_B_0[4]
        delta_L_interp = delta_L_0[4]

    # Linear interpolation between values for the two closest latitudes
    else:
        
        n = 2
        # Obtain the position of near points
        x2 = bisect_right(latitudes, LatRx)
        x1 = x2 - 1
        x = [latitudes[x1], latitudes[x2]]    
        # Obtain the parameters of near points to interpolate
        yP = [P_0[x1], P_0[x2]]
        yT = [T_0[x1], T_0[x2]]
        ye = [e_0[x1], e_0[x2]]
        yB = [Beta_0[x1], Beta_0[x2]]
        yL = [Lambda_0[x1], Lambda_0[x2]]

        ydP = [delta_P_0[x1], delta_P_0[x2]]
        ydT = [delta_T_0[x1], delta_T_0[x2]]
        yde = [delta_e_0[x1], delta_e_0[x2]]
        ydB = [delta_B_0[x1], delta_B_0[x2]]
        ydL = [delta_L_0[x1], delta_L_0[x2]]

        P_interp = lagrangeInterpolation(x, yP, LatRx, n)
        T_interp = lagrangeInterpolation(x, yT, LatRx, n)
        e_interp = lagrangeInterpolation(x, ye, LatRx, n)
        B_interp = lagrangeInterpolation(x, yB, LatRx, n)
        L_interp = lagrangeInterpolation(x, yL, LatRx, n)

        delta_P_interp = lagrangeInterpolation(x, ydP, LatRx, n)
        delta_T_interp = lagrangeInterpolation(x, ydT, LatRx, n)
        delta_e_interp = lagrangeInterpolation(x, yde, LatRx, n)
        delta_B_interp = lagrangeInterpolation(x, ydB, LatRx, n)
        delta_L_interp = lagrangeInterpolation(x, ydL, LatRx, n)

    # Value of the 5 parameters  depended on rx Latitude and day of year (Doy)
    # Const 
    #---------------------
    Dmin = 28 # Northern latitudes (Dmin = 211 for Southern latitudes)
    # Compute 5 Meteorological parameters
    #---------------------
    P = P_interp - delta_P_interp*np.cos((2*np.pi*(Doy-Dmin))/365.25)
    T = T_interp - delta_T_interp*np.cos((2*np.pi*(Doy-Dmin))/365.25)
    e = e_interp - delta_e_interp*np.cos((2*np.pi*(Doy-Dmin))/365.25)
    B = B_interp - delta_B_interp*np.cos((2*np.pi*(Doy-Dmin))/365.25)
    L = L_interp - delta_L_interp*np.cos((2*np.pi*(Doy-Dmin))/365.25)

    # Zero-altitude zenith delays terms [z_hyd, z_wet (m)]
    #--------------------------------------------------------
    # Consts:
    # [K/mbar]
    k1 = 77.604
    # [K^2/mbar]
    k2 = 382000
    # [J/(kg*K)]
    Rd = 287.054
    # [m/s^2]
    gm = 9.784

    z_hyd = ((10**-6)*k1*Rd*P) / gm
    z_wet = (((10**-6)*k2*Rd) / (gm*(L+1)-B*Rd))*(e/T)

    # [d_hyd and d_wet] estimated range delays for a satellite 
    # at 90º elevation angle (gases and water vapor)
    #--------------------------------------------------------
    # Consts:
    # [m/s^2]
    g = 9.80665

    geoidH = computeGeoidHeight(LonRx, LatRx)
    # Receiver's height (H) above mean-sea-level
    H = AltRx - geoidH 

    d_hyd = z_hyd*((1-((B*H)/T))**(g/(Rd*B)))
    d_wet = z_wet*((1-((B*H)/T))**(((g*(L+1))/(Rd*B))-1))

    # Tropospheric delay correction TC for satellite
    #--------------------------------------------------------
    TC = (d_hyd + d_wet) * TropoMpp

    # return STD
    return TC

# Compute User Airborne Sigma Ref: MOPS-DO-229D Section J.2.4
def computeSigmaAir(elev_deg, Conf):
    # k factor SIGMA_AIR_DF
    k = int(Conf["SIGMA_AIR_DF"])

    # Sigma MultiPath   [meters]    
    sigma_MP = (0.13+0.53*np.exp(-elev_deg/10))*k

    # Consider minimum or maximum signal level
    # Minimum signal level
    if float(Conf["ELEV_NOISE_TH"]) >= elev_deg:  
        # Airborne Accuracy Designator [A,B]
        if Conf["AIR_ACC_DESIG"] == 'A':
            sigma_noise_divg = 0.36*k
        else:
            sigma_noise_divg = 0.15*k
    else:
        # Maximum signal level
        if Conf["AIR_ACC_DESIG"] == 'A':
            sigma_noise_divg = 0.15*k
        else:
            sigma_noise_divg = 0.11*k

    # Depended on Equipment Class [1,2,3,4]
    if int(Conf["EQUIPMENT_CLASS"]) == 1:
        sigmaAir = 5 # [m^2]
    else:
        sigmaAir = np.sqrt(sigma_MP**2 + sigma_noise_divg**2)
    
    # return sigma Air
    return sigmaAir, sigma_MP, sigma_noise_divg

# Compute Sigma UERE by combining all Sigma contributions
def computeSigmaUERE(Conf, SigmaTROPO, SigmaAir):
    ##----------- sp3 (cm) | clk (ns) | tropo (m) | air (m^2) -------------##
    SigmaSP3 = Conf["SP3_ACC"]/100
    SigmaCLK = (Conf["CLK_ACC"]/1e9)*Const.SPEED_OF_LIGHT
    
    Sigma_UERE = np.sqrt(SigmaSP3**2 + SigmaCLK**2 + SigmaTROPO**2 + SigmaAir**2)

    return Sigma_UERE

# Compute the Geometrical Range
def computeGeomRange(SatCopPos, RcvrRefPosXyz):
    distance = np.sqrt((SatCopPos[0] - RcvrRefPosXyz[0])**2 + \
                       (SatCopPos[1] - RcvrRefPosXyz[1])**2 + \
                       (SatCopPos[2] - RcvrRefPosXyz[2])**2)
    return distance

# Estimate the Receiver Clock first guess as a weighted average of the Residuals
# (with the weights W=1/UERE2)
def estimateRcvrClk(CodeResiduals, WeightedUEREs):
    if len(CodeResiduals) != 0:
        numerator = sum(value*weight for value,weight in zip(CodeResiduals,WeightedUEREs))
        denominator = sum(WeightedUEREs)
        RcvrClock = numerator/denominator
    else: 
        RcvrClock = 0.0
    return RcvrClock

#-----------------------------END------------------------------------------#