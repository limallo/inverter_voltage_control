# Project Name: circuitSetup
# Revision: V1.0
# Author: Liam Mallamo
# Description: Script to intialise power-flow models for testing and optimisation in PowerFactory. 
# Dependencies: pfAutomation


import sys
sys.path.append(r"E:\Programs\DIgSILENT\PowerFactory 2019 SP3\Python\3.7") # Setting up path to run PowerFactory remotely
import powerfactory  

import pfAutomation as pfAuto

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

PVPenTests = [0, 30, 60, 100]
complianceTests = [0, 10, 25, 50, 75, 100] # Set of compliance % to test at

baseLoad = 1.2
threePhasePercent = 30

for model in models:
    model.Activate()
    activeProj = app.GetActiveProject()
    print("Current Model: " + activeProj.loc_name)

    opScenFolder=app.GetProjectFolder('scen')
    opScens = opScenFolder.GetChildren(0)

    for opScen in opScens: # clear existing scenarios
            opScen.Delete()

    PV = app.GetCalcRelevantObjects('*.ElmPVsys')
    numPV = len(PV)

    baseLoad = 1.2
    threePhasePercent = 30
    pfAuto.setLoads(baseLoad, threePhasePercent, loadProfile)
    pfAuto.setChars(loadProfile, pvProfile)

    pfAuto.setLoads(baseLoad, threePhasePercent, loadProfile)
    pfAuto.setChars(loadProfile, pvProfile)

    for PVPen in PVPenTests:
        sites = pfAuto.setPVPen(PV, PVPen)

        if(PVPen!= 0): # only make compliance scenarios for non-zero levels of PVPen
            for comp in complianceTests:
                opScenName = str(PVPen) + "% PVPen " + str(comp) + "% Volt-VAr"
                opScen = opScenFolder.CreateObject('IntScenario', opScenName)

                opScen.Activate()
                print("Current OpScen: " + opScen.loc_name)

                for i, line in enumerate(PV):
                    line.outserv = sites[i]

                pfAuto.setCompliance(PV, comp)

                opScen.Save()
                opScen.Deactivate()

        else: 
            opScenName = "0% PV"
            opScen = opScenFolder.CreateObject('IntScenario', opScenName)
            opScen.Activate()
            print("Current OpScen: " + opScen.loc_name)
            for i, line in enumerate(PV):  
                line.outserv = sites[i]
            opScen.Save()
            opScen.Deactivate()

    activeProj.Deactivate()





