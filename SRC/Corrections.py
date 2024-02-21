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
            '''
            if (Sod==3000):
                print("debug")
            '''
            # Compute Satellite Clock Bias (linear interpolation between closer inputs) ######RELOJES GAP !!!
            #-----------------------------------------------------------------------
            clkBias = computeSatClkBias(Sod, SatClkInfo, SatLabel)

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
            #SatComPos = applySagnac(SatComPos, FlightTime)
            
            theta = Const.OMEGA_EARTH*FlightTime

            # Compute APO in ECEF from ANTEX APOs in
            # satellite-body reference frame
            #-----------------------------------------------------------------------
            #Apo = computeSatApo(SatComPos, RcvrPos, SunPos, SatComPos, SatApoInfo)

            #

            # Apply APOs to the Satellite Position
            #SatCopPos = SatComPos + Apo

            #

            # Compute Dtr (Relativistic correction)
            #-----------------------------------------------------------------------
            #Dtr = computeDtr(SatComPos_1, SatComPos, Sod, Sod_1)

            #

            # Apply Dtr to Clock Bias
            #SatClkBias = SatClkBias + Dtr

            #

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



        else:
            SatCorrInfo["Flag"] == 0

        # prepare output
        CorrInfo[SatLabel] = SatCorrInfo
   
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
    # No need of interpolation if it's in the infoFile
    if Sod in SodList:
        clkBias = SatClkPrn[Sod]
    else:
        #SodData = np.array(list(SatClkPrn.values()))
        n = 2
        # Obtain the Sod of near points
        x1 = SodList[bisect_left(SodList, Sod)]
        x2 = SodList[bisect_right(SodList, Sod)]
        x = [x1, x2]    
        # Obtain the ClkBiases of near points to interpolate
        y1 = SatClkPrn[x1]
        y2 = SatClkPrn[x2]
        y = [y1, y2]

        # Lagrange Interpolation (n=2: Linear Interpolation)
        clkBias = lagrangeInterpolation(x, y, Sod, n)

    return clkBias

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
        if tt < 3600:
            for i in range(n):
                sod_inter[i] = SodList[i]
        elif tt > 82800:
            for i in range(n):
                sod_inter[-i-1] = SodList[-i-1]
        else:
            # Obtain positions of near points 
            pos4 = bisect_left(SodList, tt)
            pos5 = bisect_right(SodList, tt)
            # Add to the array the nearest 10 points
            for i in range((n/2)-1):
                sod_inter[(n/2)+i] = SodList[pos5 + i]
                sod_inter[((n/2)-1)-i] = SodList[pos4 - i]

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
    distance = np.sqrt((SatComPos[0] - RcvrRefPosXyz[0])**2 + \
                       (SatComPos[1] - RcvrRefPosXyz[1])**2 + \
                       (SatComPos[2] - RcvrRefPosXyz[2])**2)
   
    FlightTime = distance / Const.SPEED_OF_LIGHT

    return FlightTime

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
    P_0 = [1013.25, 1017.25, 1011.75, 1011.75, 1013.00]
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

    # Latitude and longitud of Receiver
    LatRx = Rcvr[RcvrIdx["LAT"]]
    LonRx = Rcvr[RcvrIdx["LON"]]

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
        P_interp = lagrangeInterpolation(latitudes, P_0, LatRx, 2)
        T_interp = lagrangeInterpolation(latitudes, T_0, LatRx, 2)
        e_interp = lagrangeInterpolation(latitudes, e_0, LatRx, 2)
        B_interp = lagrangeInterpolation(latitudes, Beta_0, LatRx, 2)
        L_interp = lagrangeInterpolation(latitudes, Lambda_0, LatRx, 2)

        delta_P_interp = lagrangeInterpolation(latitudes, delta_P_0, LatRx, 2)
        delta_T_interp = lagrangeInterpolation(latitudes, delta_T_0, LatRx, 2)
        delta_e_interp = lagrangeInterpolation(latitudes, delta_e_0, LatRx, 2)
        delta_B_interp = lagrangeInterpolation(latitudes, delta_B_0, LatRx, 2)
        delta_L_interp = lagrangeInterpolation(latitudes, delta_L_0, LatRx, 2)

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

    z_hyd = (10e-6*k1*Rd*P) / gm
    z_wet = ((10e-6*k2*Rd) / (gm*(L+1)-B*Rd))*(e/T)

    # [d_hyd and d_wet] estimated range delays for a satellite 
    # at 90º elevation angle (gases and water vapor)
    #--------------------------------------------------------
    # Consts:
    # [m/s^2]
    g = 9.80665

    geoidH = computeGeoidHeight(LonRx, LatRx)
    # Receiver's height (H) above mean-sea-level
    H = geoidH - Rcvr[RcvrIdx["ALT"]]

    d_hyd = z_hyd*((1-((B*H)/T))**(g/(Rd*B)))
    d_wet = z_wet*((1-((B*H)/T))**(((g*(L+1))/(Rd*B))-1))

    # Tropospheric delay correction TC for satellite
    #--------------------------------------------------------
    TC = -(d_hyd + d_wet) * TropoMpp

    # return STD
    return TC



