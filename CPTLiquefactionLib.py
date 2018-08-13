
import numpy as np
import matplotlib.pyplot as plt
from math import exp, log, sin, radians, isnan, pi
import matplotlib.gridspec as gridspec
from SBTn import returnSoilType

finesDict = {0:1,9:0,8:10,7:40,6:70,5:85,4:95,3:100,2:100,1:99} # from Yi, 2014 Table 1
dryDensity = {0:14.5,9:17,8:17,7:14.5,6:14,5:14.5,4:16,3:17,2:14.5,1:16} # estimates from local knowledge, based on typical soil types
wetDensity = {0:16,9:19,8:19,7:16.5,6:16,5:16.5,4:18,3:19,2:15.5,1:17.5} # estimates from local knowledge, based on typical soil types


def CalculateLiquefaction(a, WaterLevel, SoilsList, GraphName):
    ar = a*9.81 # peak ground accelleration
    M = 6.0 # Earthqake magnitude
    MSF = 6.9*exp((-1*M)/4)-0.058 # pg 12 eqn 2.17 
    Cfc = 0.0
    # primes some variables
    SigmaV = 0
    SigmaPrimeV = 0

    # Primes some lists to allow graphing later
    CRR_graph = []
    CSR_graph = []
    depthGraph = []
    FS_graph = []
    qc_graph = []
    strainDepth = []

    TempSoilsList = []

    for thisLineNo in range(1,len(SoilsList)):
        H, qc, fs, u2, soil = SoilsList[thisLineNo] # takes the items from the cpt list
        qc = qc * 1000
        fs = fs * 1000
        u2 = u2 * 1000

        qt = qc+u2*(1-0.85)
        Rf = (float(fs)/float(qt))*100
        if soil == 0:
            point = [Rf, qt]
            soil = returnSoilType(point)
        
        qc_graph.append(qc/1000)
        if H < WaterLevel:
            deltaDepth = SoilsList[thisLineNo][0]-SoilsList[thisLineNo-1][0]
            deltaWeight = 14.5*deltaDepth
            u0 = 0
            #eltaWeight = dryDensity[soil]*(H-SoilsList[thisLineNo-1][0]) #uses the dry density from above
        else:
            u0 = (H - WaterLevel)*9.81

            deltaDepth = SoilsList[thisLineNo][0]-SoilsList[thisLineNo-1][0]
            deltaWeight = 14.5*deltaDepth

        SigmaV = SigmaV + deltaWeight # adds the increment to the stress. Still in kN

        SigmaVu = SigmaV # convert to MN to keep in step with the rest of the calcs

        SigmaPrimeV = SigmaV - u0 # Vertical effective stress Rajapaske pg 47 ff
        SigmaPrimeVu = SigmaPrimeV 
        
            
        # some values to allow the while loop below
        # from Gregg, guide to cone penetration testing, 2015
        checkVal = 0
        count_bail = 0
        qt = qc+u2*(1-0.85)
        
        n = 1
        Q = (qt - SigmaVu)/SigmaPrimeVu
        F = (fs/(qt-SigmaVu))*100

        Ic = ((3.47-np.log10(Q))**2+(np.log10(F)+1.22)**2)**0.5
        checkN = 0
        oldN = n
        nCount = 0
        while checkN == 0:
            n = 0.381*Ic+0.05*(SigmaPrimeVu/100)-0.15
            nCount +=1
            if 0.99 < n/oldN and 1.01 > n/oldN:
                checkN = 1
            if nCount > 50:
                checkN = 1
            if n > 1:
                n = 1
            if n < 0.5:
                n = 0.5
            oldN = n
            
            Q = ((qt-SigmaVu)/100)*(100/SigmaPrimeVu)**n
            
            Ic = ((3.47-np.log10(Q))**2+(np.log10(F)+1.22)**2)**0.5
        #if Ic > 2.7:
        #    Ic = 2.7

        
        FC = 80*(Ic+Cfc)-137 

        if FC < 0:
            FC = 0
        if FC > 100:
            FC = 100

        Cn = 1 
        qc1N = qc/100
        deltaqc1N=(11.9+qc1N/14.6)*exp(1.63-9.7/(FC+2)-(15.7/(FC+2))**2) 
        qc1NCS = deltaqc1N + qc1N 
        
        oldM = 100

        while checkVal == 0:
            m = 1.338-0.249*(qc1NCS)**0.264 #2.15 b 
            Cn = (100/SigmaPrimeVu)**m # Pg 10 2.15a  
            if Cn > 1.7:
                Cn = 1.7

            qc1N = Cn*(qc/100) # 2.15b atmospheres in 
            deltaqc1N=(11.9+qc1N/14.6)*exp(1.63-9.7/(FC+2)-(15.7/(FC+2))**2) #2.22 pg 15
            qc1NCS2 = qc1N + deltaqc1N # pg 7 2.10
            #get a close on qc1NCS
            if 0.9999 < qc1NCS/qc1NCS2 and 1.0001 > qc1NCS/qc1NCS2:
                qc1NCS = qc1NCS2
                checkVal = 1
            else:
                qc1NCS = qc1NCS2
                count_bail += 1
            if count_bail > 60:
                #print ("Bailed")
                checkVal = 1

        qc1N = Cn*(qc/100)
        deltaqc1N=(11.9+qc1N/14.6)*exp(1.63-9.7/(FC+2)-(15.7/(FC+2))**2)
        qc1NCS = qc1N + deltaqc1N
        #if qc1NCS < 21:
        #    qc1NCS = 21
        if qc1NCS > 254:
            qc1NCS = 254

        #print (deltaqc1N)

        Csigma = 1/(37.3-8.27*qc1NCS**0.264) # 2.16 b pg 11
        if Csigma > 0.3:
            Csigma = 0.3

        Ksigma = 1-Csigma*log(SigmaPrimeVu/101.3) # 2.16 a pg 11

        alphaZ = -1.012-1.126*sin(H/11.73+5.133) #2.14 b pg 8
        betaZ = 0.106 + 0.118*sin(H/11.28+5.142) #2.14 c pg 8
        rd = exp(alphaZ + betaZ * M) # 2.14a pg 8
        
        
        CSR_1 = 0.65*(SigmaVu/SigmaPrimeVu)*(ar/9.81)*rd #2.2 pg 5  
        
        CSR = CSR_1/(MSF * Ksigma) # 2.6 pg 6
       
         # Pg 10, qc1NCS must be between 21 and 254
        if qc1NCS < 21:
            qc1NCS = 21
        if qc1NCS > 254:
            qc1NCS = 254
        

        if H < WaterLevel:
            #pass
            CRR = 4 # very high value so that the dry layers are non-liquefiable
        else:
        #    if qc1NCS/113+(qc1NCS/1000)**2-(qc1NCS/140)**3+(qc1NCS/137)**4-2.8 > 100:
        #        CRR = qc1NCS/113+(qc1NCS/1000)**2-(qc1NCS/140)**3+(qc1NCS/137)**4-2.8 
        #    else:
             CRR = exp(qc1NCS/113+(qc1NCS/1000)**2-(qc1NCS/140)**3+(qc1NCS/137)**4-2.8) # 2.24 pg 17
        if FC > 90:
            CRR = 4 
       
         
        FS = round(CRR/CSR,1)
       
        #print (H, u0, SigmaVu, SigmaPrimeVu, rd, MSF, Ksigma, CSR, Ic, m, Cn, qc1NCS, CRR, FS) 
        # From Zhang, Robertson, Brachman, 2008 appendix A
        if FS >= 2:
            EpsV = 0

        if FS >= 1.3 and FS < 2:
            EpsV = 7.6*qc1NCS**(-0.71)

        if FS == 1.2:
            EpsV = 9.7*qc1NCS**(-0.69)

        if FS == 1.1:
            EpsV = 11*qc1NCS**(-0.65)

        if FS == 1.0:
            EpsV = 64*qc1NCS**(-0.93)

        if FS == 0.9 and qc1NCS > 60:
            EpsV = 1430*qc1NCS**(-1.48)
        if FS == 0.9 and qc1NCS < 60:
            EpsV = 102*qc1NCS**(-0.82)

        if FS == 0.8 and qc1NCS > 80:
            EpsV = 1690*qc1NCS**(-1.46)
        if FS == 0.8 and qc1NCS < 80:
            EpsV = 102*qc1NCS**(-0.82)

        if FS == 0.7 and qc1NCS > 110:
            EpsV = 1701*qc1NCS**(-1.42)
        if FS == 0.7 and qc1NCS < 110:
            EpsV = 102*qc1NCS**(-0.82)

        if FS == 0.6 and qc1NCS > 147:
            EpsV = 2411*qc1NCS**(-1.45)

        if FS == 0.6 and qc1NCS < 147:
            EpsV = 102*qc1NCS**(-0.82)

        if FS <= 0.5:
            EpsV = 102*qc1NCS**(-0.82)
        if EpsV > 100:
            EpsV = 95
        if H < WaterLevel:
            EpsV = 0
                
        if Ic <= 1.64:
            Kc = 1
        if Ic > 1.64:
            Kc = 5.581*Ic**3-0.430*Ic**4-21.63*Ic**2+33.75*Ic-17.88

        N1_60_cs = (Q*Kc)/(8.5*(1-Ic/4.6))
        Nc = (M-4)**2.17
        EpsilonVol15 = wetDensity[soil]/1000*(N1_60_cs/20)**(-1.20)
        EpsilonVol = EpsilonVol15*(Nc/15)**0.45
        
        CRR_graph.append(CRR)
        CSR_graph.append(CSR)
        depthGraph.append(-1*H)
        FS_graph.append(FS)
        Strain = EpsV/100*(H-SoilsList[thisLineNo-1][0])*1000 # pg 1171, 1
        
        
        if Strain > 20: #just making sure that the strain does not exceed the layer thickness
            Strain = 10

        if isnan(Strain):
            strainDepth.append(0)
        else:
            strainDepth.append(Strain)
        if H < WaterLevel:
            CRRReturn = -1
        else:
            CRRReturn = CRR

#        if H < 10:
#            print (Strain)

        shortList = [H, CRRReturn, CSR, FS, Strain]
        TempSoilsList.append(shortList)
        print (H, u0, SigmaV, SigmaPrimeV, MSF, Ksigma, CSR, qt, Ic, Cn, qc1N, qc1NCS, CRR, FS)
    print (" ")

    CRR_graph = np.array(CRR_graph)
    CSR_graph = np.array(CSR_graph)
    depthGraph = np.array(depthGraph)
    qc_graph = np.array(qc_graph)
    strainDepth = np.array(strainDepth)
    SumSettle = np.sum(strainDepth)
    FS_graph = np.array(FS_graph)
   
    # adds up the settlements
    ReturnList = []
    count = -1
    SettleForGraph = []
    for thisLine in TempSoilsList:
        count += 1
        SettleItem = np.sum(strainDepth[count:])
        SettleForGraph.append(SettleItem)
        thisLine.append(SettleItem)
        ReturnList.append(thisLine)

    SettleForGraph = np.array(SettleForGraph)
    fig = plt.figure()

    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    ax1.plot(CRR_graph, depthGraph)
    ax1.plot(CSR_graph, depthGraph)
    ax1.set_xlim(0,0.6)
    ax1.set_title("CRR and CSR")

    ax2.plot(FS_graph, depthGraph)
    ax2.set_xlim(0,2)
    ax2.set_title("Factor of Safety")

    ax3.plot(SettleForGraph, depthGraph)
    ax3.set_title("Settlements (mm)")

    TitleText = "Liquefaction and settlement plots for "+str(round(ar/9.81,3))+"g"
    plt.suptitle(TitleText, fontsize=20)
    plt.savefig(GraphName)
    

    return(SumSettle, ReturnList)

