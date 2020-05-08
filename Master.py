# Project Name: Master
# Revision: V1.0
# Author: Liam Mallamo
# Description: Central script for inverter Volt-VAr optimisation project. 
# Dependencies: numpy, scipy, pyswarms, pfAutomation

import sys
sys.path.append(r"E:\Programs\DIgSILENT\PowerFactory 2019 SP3\Python\3.7") # Setting up path to run PowerFactory remotely
import powerfactory  
import os
import numpy as np
import math
from pyswarm import pso

import pfAutomation as pfAuto

import time

def fitnessFunction(newCurveSetpoints):
    global exportPath, activeProj, opScen, PV, Vmin, Vmax, loadProfile, pvProfile, NMIs, pvSize, vNom

    pfAuto.setChars(loadProfile,pvProfile) # Reset PV & Load profiles to defaults
    newCurve = pfAuto.generateCurve2(pvSize,vNom,newCurveSetpoints)

    for line in PV:
        pfAuto.setVVcurve(newCurve, line)

    pfAuto.runSim()
    pfAuto.deratePV2(PV, user)

    pfAuto.runSim()
    Ploss = pfAuto.calcCurtailment(PV, user)

    Vmin, Vmax = pfAuto.calcVrange(NMIs)

    return Ploss

def constraintFunction(newCurveSetpoints):
    global Vmax
    return [253-Vmax]
    
start = time.time()

Vmax = 0
Vmin = 0

####################### SET INPUT PARAMETERS #################################
pvSize = 5
vNom = 230

ts129Setpoints = [207, 220, 248, 253, 0.9, 0.95] # Curve setpoints as per TS129
ts129Curve = pfAuto.generateCurve(pvSize, vNom, ts129Setpoints) # Generate baseline curve

#################### START POWERFACTORY & OPEN MODELS ###########################
app=powerfactory.GetApplication() # Start PowerFactory in Unattended Mode
user=app.GetCurrentUser() 

# Get load and PV profiles from master project library
charLib = user.GetContents(r"Inverter Voltage Control\Library\Characteristics")[0]
chars = charLib.GetContents('*.ChaTime')
loadProfile = chars[0]
pvProfile = chars[1]

# Choosing folders
projFolder = user.GetContents(r"Inverter Voltage Control\Complete Models")[0]
models = projFolder.GetChildren(0)

resultsFiles = []

charLib = user.GetContents(r"Inverter Voltage Control\Library\Characteristics")[0]
chars = charLib.GetContents('*.ChaTime')
loadProfile = chars[0]
pvProfile = chars[1]
    
# Run for all models
for model in models:
    ############# Choose project / operational scenario ###################
    model.Activate()
    activeProj = app.GetActiveProject()
    print("Current Model: " + activeProj.loc_name)
    exportPath = os.path.join(r"E:\Liam\OneDrive for Business\OneDrive - University of South Australia\2019\Honours\Modelling\Results", activeProj.loc_name) # Directory for exporting results

    opScenFolder = app.GetProjectFolder('scen')
    opScens = opScenFolder.GetChildren(0)

    opScenNames = []
    for opScen in opScens:
        opScenNames.append(opScen.loc_name)

    PV = app.GetCalcRelevantObjects('*.ElmPVsys')
    terms = app.GetCalcRelevantObjects('*.ElmTerm')
    loads = app.GetCalcRelevantObjects('*.ElmLod')

    NMIs = []
    for term in terms:
        if('Terminal' not in term.loc_name):
            NMIs.append(term)

    initVmax = 0
    initVmin = 0

    for opScen in opScens:
        opScen = opScens[opScenNames.index("60% PVPen 0% Volt-VAr")] # 60% PV 50% Volt-VAr
        opScen.Activate()
        print("Current OpScen: " + opScen.loc_name)

        if(opScen.loc_name == "0% PV"): # Test without PV
            pfAuto.runSim() # Run sim
            resultsFile = pfAuto.exportResults(exportPath, activeProj, opScen, "")

        else:
            if(" 0% Volt-VAr" in opScen.loc_name): # Test without Volt-VAr
                pfAuto.runSim() # Run sim
                initVmin, initVmax = pfAuto.calcVrange(NMIs)
                pfAuto.exportResults(exportPath, activeProj, opScen, "") # Export results to CSV

            else:
                ######################### Test with baseline curve ############
                for line in PV: 
                    pfAuto.setVVcurve(ts129Curve, line) #Apply base curve to all systems
                pfAuto.runSim() # Run the 24hr sim
                pfAuto.deratePV2(PV, user) # apply derated PV profile
                pfAuto.runSim() # re-run the sim
                minV, maxV = pfAuto.calcVrange(NMIs)

                Ploss = pfAuto.calcCurtailment(PV, user)
                print("Total circuit curtailment with TS129 is " + str(round(Ploss,3)) + "kWh")

                pfAuto.exportResults(exportPath, activeProj, opScen, "TS129") # Export curtailed results to CSV
                pfAuto.setChars(loadProfile,pvProfile) # Reset PV & Load profiles to defaults

                #######################################################################################
                ################################## OPTIMISATION #######################################
                #######################################################################################
                if(initVmax >= 252 or initVmin <= 216):
                    print('Peak voltage of ' + str(round(initVmax, 3)) + 'V detected, proceeding to optimise with Volt-VAr.')

                    lb = [249, 251, 0.8]
                    ub = [250, 253,  1]

                    xopt, fopt = pso(fitnessFunction, lb, ub, f_ieqcons=constraintFunction, swarmsize=24, maxiter = 25, omega = 0.9, phip = 0.5, phig = 0.3, minfunc = 0.01,  debug = True)
                    
                    optCurve = pfAuto.generateCurve2(pvSize,vNom,xopt)

                    for line in PV:
                        pfAuto.setVVcurve(optCurve, line)

                    pfAuto.runSim()
                    pfAuto.deratePV2(PV, user)

                    pfAuto.runSim()
                    pfAuto.exportResults(exportPath, activeProj, opScen, 'optCurve')

                    pfAuto.setCompliance(PV, 0) # Reset so no PV systems have Volt-VAr
                    pfAuto.setChars(loadProfile,pvProfile) # Reset PV & Load profiles to defaults

                else:
                    print("Max voltage for this PV penetration without Volt-VAr is " + str(round(initVmax, 3)) + "V, no optimisation needed.")

        opScen.Save()
        opScen.Deactivate()

    activeProj.Deactivate()

    # Print optimisation runtime
    end = time.time()
    print(end - start)


