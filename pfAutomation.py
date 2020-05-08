# Project Name: pfAutomation
# Revision: V1.0
# Author: Liam Mallamo
# Description: Function script for interaction with PowerFactory
# Dependencies: numpy, matplotlib, scipy, pandas

import random
import os
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt

import sys
sys.path.append(r"E:\Programs\DIgSILENT\PowerFactory 2019 SP3\Python\3.7")

import powerfactory  
app=powerfactory.GetApplication() # Start PowerFactory in Unattended Mode

#################################### CIRCUIT SETUP ####################################################

# Function to set PV penetration
def setPVPen(PV, PVPen):  
    numPV = len(PV)
    sites = []
    PVCount = 0

    if(PVPen != 100 and PVPen != 0):
            while ((PVCount < int(0.01 * PVPen * numPV)) & (len(sites) != numPV)):
                    PVCount = 0
                    sites.clear()

                    while (PVCount < int(0.01 * PVPen * numPV)):
                            isPV = random.getrandbits(1)

                            if isPV == 1:
                                    sites.append(0)
                                    PVCount += 1
                            else:
                                    sites.append(1)

                    while (len(sites) < numPV + 1) :
                            sites.append(1)

                    # for j, line in enumerate(PV):
                    #         line.outserv = sites[j] 
                                              
    elif(PVPen == 100):
        for line in PV:
            PVCount += 1
            sites.append(0)

    else:
        for line in PV:
            sites.append(1)

    print("PV Penetration for OpScenario " + str(PVPen) +" = " + str((PVCount/numPV) * 100) + "%")

    return sites

# Function to set required Volt-VAr compliance
def setCompliance(PV, comp):
    
    numActivePV = 0
    vvCount = 0

    pvSites = np.zeros(len(PV), dtype = 'object')
    vvSites = []

    for i, line in enumerate(PV):
        if line.outserv == 0:
            numActivePV += 1
            vvCount += 1
            pvSites[i] = "qvchar"
            vvSites.append(i)
        else:
            pvSites[i] = "constq"

    if(comp != 100):
        while (vvCount > round(0.01 * comp * numActivePV)): 
            removeVV = random.choice(vvSites)
            vvSites.remove(removeVV)
            vvCount-=1
            app.PrintPlain("Compliance =  " + str((vvCount/numActivePV) * 100) + "%")                

    for j, line in enumerate(pvSites):
        if j in vvSites:
            pvSites[j] = "qvchar"
        else:
            pvSites[j] = "constq"

    for k, line in enumerate(PV):
        line.av_mode = pvSites[k]

    if(comp!=0):
        print("Volt-VAr compliance for " + str(comp) + "% test = " + str((vvCount/numActivePV) * 100) + "%")     

    return vvCount  

# Function to set & scale loads and phasing
def setLoads(baseLoad, threePhasePercent , loadProfile):
        loads = app.GetCalcRelevantObjects('*.ElmLod') #get list of all loads
        numLoads = len(loads)

        loadPhases = []
        loadScale = []

        TPCount = 0

        for line in loads:
                while ((TPCount < (threePhasePercent * 0.01 * numLoads)) & (len(loadPhases) != numLoads)):
                        TPCount = 0
                        loadPhases.clear()

                        while (TPCount < (threePhasePercent * 0.01 * numLoads)):
                                is3PH = random.getrandbits(1)
                                loadPhases.append(is3PH)

                                if is3PH == 0:
                                        TPCount += 1
                        
                        while (len(loadPhases) < numLoads + 1) :
                                loadPhases.append(1)

        #Setting balanced/unbalanced for loads
        j = 0
        for line in loads:
                line.i_sym = loadPhases[j]
                j += 1

        #Setting phases for unbalanced loads
        phase1 = []
        phase2 = []
        phase3 = []

        j = 0
        while (j < numLoads):
                phase1.insert(j, baseLoad)
                phase1.insert(j + 1, 0)
                phase1.insert(j + 2, 0)

                phase2.insert(j, 0)
                phase2.insert(j + 1, baseLoad)
                phase2.insert(j + 2, 0)

                phase3.insert(j, 0)
                phase3.insert(j + 1, 0)
                phase3.insert(j + 2, baseLoad)

                j = j + 3
        j = 0
        for line in loads:
                        line.plini = baseLoad
                        line.plinir = phase1[j]
                        line.plinis = phase2[j]
                        line.plinit = phase3[j]

                        j+=1

        ####################### Scaling Loads #####################################
        while(len(loadScale) < numLoads + 1):
                loadScale.append(0.75 + random.uniform(0, 1))

        j = 0

        for line in loads:
                line.scale0 = loadScale[j]
                j += 1

# Function to set generation and load profiles
def setChars(loadProfile, pvProfile):
        loads = app.GetCalcRelevantObjects('*.ElmLod') #get list of all loads
        PV = app.GetCalcRelevantObjects('*.ElmPVsys') #get list of all PV systems

        for line in loads:
                # Deleting old characteristic
                oldChar = line.GetContents('*.Cha*')
                for i in oldChar:
                        i.Delete()
                                        
                # Setting new characteristic
                newCharBal = line.CreateObject('ChaRef','plini')
                newCharUnBal1 = line.CreateObject('ChaRef', 'plinir')
                newCharUnBal2 = line.CreateObject('ChaRef', 'plinis')
                newCharUnBal3 = line.CreateObject('ChaRef', 'plinit')

                newCharBal.typ_id = loadProfile
                newCharUnBal1.typ_id = loadProfile
                newCharUnBal2.typ_id = loadProfile
                newCharUnBal3.typ_id = loadProfile

        for line in PV:
                # Deleting old characteristic
                oldChar = line.GetContents('*.Cha*')
                for i in oldChar:
                        i.Delete()
                                        
                # Setting new characteristic
                newChar = line.CreateObject('ChaRef','pgini')
                newChar.typ_id = pvProfile


#################################### VOLT-VAR CURVES AND SIMULATION ########################################

# Function to generate a baseline curve from input constraints
def generateSeedCurve(pvSize, vNom, curveConstraints):
    
    V1 = np.median(curveConstraints[0])
    V2 = np.median(curveConstraints[1]) 
    V3 = np.median(curveConstraints[2])
    V4 = np.median(curveConstraints[3])
    PFlag = -1*np.median(curveConstraints[4])
    PFlead = np.median(curveConstraints[4])
    
    Qmin_per = -1*math.sqrt(1**2-PFlag**2)*100
    Qmax_per = math.sqrt(1**2-PFlead**2)*100

    Qmin_kvar = round(-1*math.sqrt(pvSize**2 - (PFlag * pvSize)**2),3)
    Qmax_kvar = round(math.sqrt(pvSize**2 - (PFlead * pvSize)**2),3)

    dbLow = round(V2/vNom,3)
    dbHigh = round(V3/vNom,3)

    indDroop = round((Qmin_kvar/pvSize)/(V4/vNom-dbHigh),3)
    capDroop = round((Qmax_kvar/pvSize)/(dbLow-V1/vNom),3)

    seedCurve = [Qmin_per, Qmax_per, Qmin_kvar, Qmax_kvar, dbLow, dbHigh, indDroop, capDroop]

    return seedCurve  

# Function to generate a curve using setpoints 2
def generateCurve(pvSize, vNom, curveSetpoints):
    V1 = curveSetpoints[0]
    V2 = curveSetpoints[1]
    V3 = curveSetpoints[2]
    V4 = curveSetpoints[3]

    PFlag = -1*curveSetpoints[4]
    PFlead = curveSetpoints[5]
   
    Qmin_per = -1*math.sqrt(1**2-PFlag**2)*100
    Qmax_per = math.sqrt(1**2-PFlead**2)*100

    Qmin_kvar = round(-1*math.sqrt(pvSize**2 - (PFlag * pvSize)**2),3)
    Qmax_kvar = round(math.sqrt(pvSize**2 - (PFlead * pvSize)**2),3)

    dbLow = round(V2/vNom,3)
    dbHigh = round(V3/vNom,3)

    indDroop = round((Qmin_per)/(V4-V3),3)
    capDroop = round((Qmax_per)/(V2-V1),3)

    curveParams = [Qmin_per, Qmax_per, Qmin_kvar, Qmax_kvar, dbLow, dbHigh, indDroop, capDroop]

    return curveParams

# Function to generate a curve using setpoints 2
def generateCurve2(pvSize, vNom, curveSetpoints):
    #V1 = curveSetpoints[0]
    #V2 = curveSetpoints[1]
    # V3 = curveSetpoints[2]
    # V4 = curveSetpoints[3]
    V1 = 207   
    V2 = 220
    V3 = curveSetpoints[0]
    V4 = curveSetpoints[1]
    
    #PFlag = -1*curveSetpoints[4]
    #PFlead = curveSetpoints[5]
    PFlag = -1*curveSetpoints[2]
    PFlead = 0.9

    Qmin_per = -1*math.sqrt(1**2-PFlag**2)*100
    Qmax_per = math.sqrt(1**2-PFlead**2)*100

    Qmin_kvar = round(-1*math.sqrt(pvSize**2 - (PFlag * pvSize)**2),3)
    Qmax_kvar = round(math.sqrt(pvSize**2 - (PFlead * pvSize)**2),3)

    dbLow = round(V2/vNom,3)
    dbHigh = round(V3/vNom,3)

    indDroop = round((Qmin_per)/(V4-V3),3)
    capDroop = round((Qmax_per)/(V2-V1),3)

    curveParams = [Qmin_per, Qmax_per, Qmin_kvar, Qmax_kvar, dbLow, dbHigh, indDroop, capDroop]

    return curveParams

# Function to apply given Volt-VAr settings
def setVVcurve(curveParams, PVsys):
    PVsys.iopt_tdr = 1
    PVsys.ddroopue = curveParams[7]
    PVsys.ddroopoe = curveParams[6]
    PVsys.Qfu_max = curveParams[3]
    PVsys.Qfu_min = curveParams[2]
    PVsys.udeadblow = curveParams[4]
    PVsys.udeadbup = curveParams[5]

# Function to run a 24hr load-flow simulation
def runSim(): 
    dynSim=app.GetFromStudyCase('ComStatsim')
    dynSim.calcPeriod=0 #Time period = Complete Day
    dynSim.stepSize=30#Step Size Step = 30
    dynSim.stepUnit=1 #Step Size Unit = minutes
    dynSim.iEnableParal = 1
    dynSim.Execute()

############################################## RESULTS CALCULATIONS ######################################

# Function to export simulation results to CSV - note: export variables must be set up manually in PowerFactory, maybe look into automating this, time permitting
def exportResults(exportPath, activeProj, opScen, currentCurve):
    if not os.path.exists(exportPath): 
        os.makedirs(exportPath) 

    comRes = app.GetFromStudyCase("ComRes")
    elmRes = app.GetFromStudyCase("Quasi-Dynamic Simulation AC unbalanced.ElmRes")

    comRes.pResult = elmRes
    comRes.iopt_exp = 6 # to export as csv
    comRes.f_name = os.path.join(exportPath , activeProj.loc_name + " " + opScen.loc_name + " " + currentCurve + ".csv") #File Name
    comRes.iopt_time = 3 # set time format to hh:mm:ss
    comRes.iopt_sep = 1 # to use the system seperator
    comRes.iopt_honly = 0 # to export data and not only the header
    comRes.iopt_csel = 1 # export only selected variables


    # Only get NMIs, not midline terminals - this only works for models with NMIs, so the Tesla sites - REVISE IT!!!
    terms = app.GetCalcRelevantObjects("*.ElmTerm")
    loads = app.GetCalcRelevantObjects('*.ElmLod')
    PV = app.GetCalcRelevantObjects("*.ElmPVSys")
    numNMIs=len(loads)
    numTerms=len(terms)

    del terms[numNMIs:numTerms]

    elements = []
    variable = []
    resultObject = []

    elements.append(elmRes)
    variable.append("b:ucttime")
    resultObject.append(elmRes)

    for line in terms:
        elements.append(line)
        resultObject.append(elmRes)
        variable.append("m:U")

    for line in loads:
        elements.append(line)
        resultObject.append(elmRes)
        variable.append("m:Psum:bus1")

        elements.append(line)
        resultObject.append(elmRes)
        variable.append("m:Qsum:bus1")

        elements.append(line)
        resultObject.append(elmRes)
        variable.append("m:Ssum:bus1")

    for line in PV:
        elements.append(line)
        resultObject.append(elmRes)
        variable.append("m:Psum:bus1")

        elements.append(line)
        resultObject.append(elmRes)
        variable.append("m:Qsum:bus1")

        elements.append(line)
        resultObject.append(elmRes)
        variable.append("m:Ssum:bus1")

    comRes.resultobj = resultObject
    comRes.element = elements
    comRes.variable = variable

    comRes.Execute()

    return comRes.f_name

# Function to import CSV results and produce voltage and curtailment matrices
def importResults(resultsFile, nodes):
    pv_profile = [0,0,0,0,0,0,0,0,0,0,0,0,0,0.003,0.026,0.084,0.191,0.314,0.439,0.551,0.66,0.744,0.783,0.892,0.943,
    0.958,0.99,1,0.949,0.899,0.835,0.724,0.633,0.524,0.403,0.273,0.144,0.052,0.007,0.001,0.001,0,0,0,0,0,0,0]

    data = pd.read_csv(resultsFile, header = [0,1])
    data.rename(columns={"Quasi-Dynamic Simulation AC unbalanced": "Time"}, inplace = True)
    data.set_index("Time", inplace = True)

    # Get the node name of the worst voltage on the circuit at 13:30 (peak PV production)
    eolNode = (data.idxmax(axis = 1)[28])[0]

    eolVoltage = data[eolNode, 'U, Magnitude in V']

    totalCurtailment = 0

    for i, node in enumerate(nodes):
        if(node.outserv == 0):
            Prated = node.scale0 * node.sgn

            siteP = data[node.loc_name, 'Total Active Power in kW']
            
            for j, P in enumerate(siteP):
                if(node.av_mode == 'qvchar'):
                    totalCurtailment += abs(((pv_profile[j] * Prated) - P))/2 # curtailment of that point in kWh
            
    return (eolVoltage.min(), eolVoltage.max(), totalCurtailment)

# Function to apply modified PV profiles to constrained Volt-VAr systems
def deratePV(resultsFile, PV, user): 
    data = pd.read_csv(resultsFile, header = [0,1])
    data.rename(columns={"Quasi-Dynamic Simulation AC unbalanced": "Time"}, inplace = True)
    data.set_index("Time", inplace = True)

    ################### Getting baseline PV profile ###################
    charLib = user.GetContents(r"Inverter Voltage Control\Library\Characteristics")[0]
    chars = charLib.GetContents('*.ChaTime')
    pvProfile = chars[1].vector

    charFold=app.GetProjectFolder('chars') # folder to store derated characteristics

    oldChars = charFold.GetContents('*.Cha*')
    for i in oldChars:
            i.Delete()

    for line in PV:
        if(line.av_mode == 'qvchar'):
            S_rated = line.scale0 * line.sgn

            siteP = data[line.loc_name, 'Total Active Power in kW']
            siteQ = data[line.loc_name, 'Total Reactive Power in kvar']
            siteS = data[line.loc_name, 'Total Apparent Power in kVA']
        
            pvProfileCurtailed = pvProfile

            isCurtailed = 0

            for i, S in enumerate(siteS):
                if (S > S_rated):
                    P_derated = round(math.sqrt(S_rated**2 - abs(siteQ[i])**2),2)
                    pvProfileCurtailed[i] = P_derated/S_rated
                    isCurtailed = 1

                    #print(line.loc_name + ' is curtailed at index ' + str(i))

            if (isCurtailed == 1):
                curtailedChar = charFold.CreateObject('ChaTime', line.loc_name + ' Profile')
                curtailedChar.stepSize = 30
                curtailedChar.usage = 1
                curtailedChar.approx = 0
                curtailedChar.vector = pvProfileCurtailed

                newChar = line.GetContents('*.Cha*')[0]
                newChar.typ_id = curtailedChar

def deratePV2(PV, user):
    ################### Getting baseline PV profile ###################
    charLib = user.GetContents(r"Inverter Voltage Control\Library\Characteristics")[0]
    chars = charLib.GetContents('*.ChaTime')
    pvProfile = chars[1].vector

    charFold=app.GetProjectFolder('chars') # folder to store derated characteristics

    oldChars = charFold.GetContents('*.Cha*')
    for i in oldChars:
            i.Delete()

    for line in PV:
        if(line.av_mode == 'qvchar'):
            S_rated = line.scale0 * line.sgn
            powers = getPowers(line)
            siteS = powers[:,2]
            siteQ = powers[:,1]

            pvProfileCurtailed = pvProfile

            isCurtailed = 0
            
            for i, S in enumerate(siteS):
                if (S > S_rated):
                    P_derated = round(math.sqrt(S_rated**2 - abs(siteQ[i])**2),2)
                    pvProfileCurtailed[i] = P_derated/S_rated
                    isCurtailed = 1

                    #print(line.loc_name + ' is curtailed at index ' + str(i))

            if (isCurtailed == 1):
                curtailedChar = charFold.CreateObject('ChaTime', line.loc_name + ' Profile')
                curtailedChar.stepSize = 30
                curtailedChar.usage = 1
                curtailedChar.approx = 0
                curtailedChar.vector = pvProfileCurtailed

                newChar = line.GetContents('*.Cha*')[0]
                newChar.typ_id = curtailedChar

def getPowers(element):
    resultsFile = app.GetFromStudyCase("*.ElmRes")
    resultsFile.Load()
    nRows = app.ResGetValueCount(resultsFile, 0)
    values = np.empty([nRows, 3])

    variables = ["m:Psum:bus1", "m:Qsum:bus1", "m:Ssum:bus1"]

    col = 0
    for var in variables:
        colIndex = app.ResGetIndex(resultsFile, element, var)
        for i in range(nRows):
            values[i,col] = resultsFile.GetValue(i, colIndex)[1]
        col += 1

    return values

def calcVrange(terminals):
    resultsFile = app.GetFromStudyCase("*.ElmRes")
    resultsFile.Load()
    nRows = app.ResGetValueCount(resultsFile, 0)
    nCols = len(terminals)
    values = np.empty([nRows, nCols])

    col = 0
    for term in terminals:
        colIndex = app.ResGetIndex(resultsFile, term, 'm:U')
        for i in range(nRows):
            values[i,col] = resultsFile.GetValue(i, colIndex)[1]
        col += 1

    minV = np.amin(values)
    maxV = np.amax(values)

    return minV, maxV
    #return [252 - Vmax]

def calcCurtailment(PVsystems, user):
    charLib = user.GetContents(r"Inverter Voltage Control\Library\Characteristics")[0]
    chars = charLib.GetContents('*.ChaTime')
    pvProfile = chars[1].vector

    totalCurtailment = 0

    for i, sys in enumerate(PVsystems):
        powers = getPowers(sys)
        if(sys.outserv == 0):
            Prated = sys.scale0 * sys.sgn

            siteP = powers[:,0]
            
            for j, P in enumerate(siteP):
                if(sys.av_mode == 'qvchar'):
                    totalCurtailment += abs(((pvProfile[j] * Prated) - P))/2 # curtailment of that point in kWh
    
    return totalCurtailment


################################################ PLOTS ########################################################
 
# Function to plot curve
def plotCurve(curveParams):
    Qmin = curveParams[0]
    Qmax = curveParams[1]
    dbLow_pu = curveParams[4]
    dbHigh_pu = curveParams[5]
    dbHigh_volts = round(dbHigh_pu*230,1)
    dbLow_volts = round(dbLow_pu*230,1)
    indDroop = curveParams[6]
    capDroop = curveParams[7]

    V = np.arange(170,270,0.1)
    Q  = np.zeros([len(V),1])

    for i, voltage in enumerate(V):

        if voltage < dbLow_volts:     # capacitive side
            if (dbLow_volts-voltage)*capDroop >= Qmax: 
                Q[i] = Qmax 
            else:
                Q[i] = (dbLow_volts-voltage)*capDroop

        elif dbLow_volts <= voltage <= dbHigh_volts: # deadband
            Q[i] = 0
 
        else: # inductive side
            if(dbHigh_volts-voltage)*-indDroop <= Qmin:
                Q[i] = Qmin

            else:
                Q[i] = (dbHigh_volts-voltage)*-indDroop

    plt.plot(V,Q)
    plt.axis([200, 270, -80, 80])
    plt.title("Volt-VAr Curves")
    plt.grid(True)
    plt.xlabel("Voltage (V)")
    plt.ylabel("VAr/Rated VA (%)")
    plt.show()

# Function to make a matplotlib plot of results
# Make one to plot no VV vs each compliance level for each circuit - so 4 plots per circuit
def createPlot(resultsFile, terms):
    data = pd.read_csv(resultsFile, header = [0,1])
    data.rename(columns={"Quasi-Dynamic Simulation AC unbalanced": "Time"}, inplace=True)
    data.set_index("Time", inplace= True)

    eol = data.idxmax(axis = 1)[28]

    nodes = []
    for node in terms:
        nodes.append(node.loc_name)

    data.plot(y=nodes, legend = None)
    plt.show()

