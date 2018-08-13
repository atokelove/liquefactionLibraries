
import numpy as np
import math as m
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from SBTn import returnSoilType

def Boussineq(length, breadth, pressure, depth):
    # Newmark's method, translated, in part, from the fortran in Bowles, 1988
    force = 0
    cenPtLen = length/2
    cenPtWidth = breadth/2
    M = cenPtWidth/depth
    N = cenPtLen/depth
    V = M**2 + N**2 + 1
    V1 = (M*N)**2
    C4 = m.atan((2*M*N*V**0.5)/(V-V1))

    if C4 <  0:
        C4 = C4 + m.pi
        
    qv = pressure*1/m.pi*((2*M*N*V**0.5)/(V+V1)*(V+1)/V+C4)

    return(qv)

def CalcSettlement(WaterLevel, load, SoilsList, footingWidth, footingLength, longTermTime, GraphName):
    dryDensity = {0:12,9:17,8:17,7:14,6:14,5:14,4:16,3:17,2:14,1:16} # estimates from local knowledge, based on typical soil types
    wetDensity = {0:15,9:19,8:19,7:16,6:16,5:16,4:18,3:19,2:15.5,1:17} # estimates from local knowledge, based on typical soil types
    load = float(load) / 1000.000
    pa = 0.1013
    settlements = []
    longTermSettlements = []
    MList = []
    qtList = []
    depthList = []
    SigmaVo = 0

    SoilsData = SoilsList 
    
    for thisLine in np.arange(1,len(SoilsData)):
        thisSoil = SoilsData[thisLine]
        depthList.append(thisSoil[0])
        dZ = SoilsData[thisLine][0]-SoilsData[thisLine-1][0]
        soil = thisSoil[4]
        if soil == 0:
            qc = thisSoil[1]
            u2 = thisSoil[3]
            fs = thisSoil[2]
            qt = qc+u2*(1-0.85)
            Rf = (float(fs)/float(qt))*100
            point = [Rf, qt]
            soil = returnSoilType(point)
        
        if thisSoil[0] > WaterLevel:
            SoilGamma = wetDensity[soil]
        else:
            SoilGamma = dryDensity[soil]

        SigmaVo = SigmaVo + dZ * SoilGamma/1000
        
        fs = thisSoil[2]
        qt = thisSoil[1]+thisSoil[3]*(1-0.85)
        
        qtList.append(qt)
        Uo = 9.81*(thisSoil[0]-WaterLevel)/1000
        if Uo < 0:
            Uo = 0
        SigmaPrimeVo = SigmaVo - Uo
        if SigmaPrimeVo == 0:
            SigmaPrimeVo = 0.000001
        Qt1 = (qt-SigmaVo)/SigmaPrimeVo
        if Qt1 < 0:
            Qt1 = 0.000001
        Fr = fs/(qt-SigmaVo)*100
        if Fr < 0.000001:
            Fr = 0.000001
        if  ((3.47-m.log(Qt1,10))**2+ (m.log(Fr,10)+1.22)**2) < 0.000001:
            Ic = 0.000001
            print(Ic)
        else:
            Ic = ((3.47-m.log(Qt1,10))**2+ (m.log(Fr,10)+1.22)**2)**0.5
        
        if Ic < 1.64:
            n = 0.5
        if Ic > 3.3:
            n = 1
        if Ic >= 1.64 and Ic <= 3.3:
            n = (Ic - 1.64)* 0.3 + 0.5
        Qtn = ((qt-SigmaVo)/pa)*(pa/SigmaPrimeVo)**n
        IcOld = Ic
        QtnOld = Qtn
        dn = 1
        count = 0
        checkVal = 0
        while checkVal == 0:
            nOld = n
            count = count + 1
            if Qtn < 0.000001:
                Qtn = 0.000001
                
            if  ((3.47-m.log(Qtn,10))**2+ (m.log(Fr,10)+1.22)**2) < 0.000001:
                Ic = 0.000001            
            else:
                Ic = ((3.47-m.log(Qtn,10))**2+ (m.log(Fr,10)+1.22)**2)**0.5

            if Ic < 1.64:
                n = 0.5
            if Ic > 3.3:
                n = 1
            if Ic >= 1.64 and Ic <= 3.3:
                n = (Ic - 1.64)* 0.3 + 0.5
            Qtn = ((qt-SigmaVo)/pa)*(pa/SigmaPrimeVo)**n
            if count == 1000:
                checkVal = 1
                print("counted out")
            deltaN = abs(n - nOld)
            if deltaN < 0.001:
                checkVal = 1
        if Ic > 2.2:
            if Qtn < 14:
                alphaM = Qtn
            if Qtn> 14:
                alphaM = 14
        if Ic < 2.2:
            alphaM = 0.0188*(10**(0.55*Ic+1.68))
            
        M = alphaM*(qt-SigmaVo)
        MList.append(M)
        deltaSigmaV = Boussineq(footingLength, footingWidth, load, float(thisSoil[0]))
        #rint (SigmaVo, deltaSigmaV)
        #print thisSoil[0], Qtn, M
        S1 = (deltaSigmaV/M) * dZ
       
        cAlpha = 0.1*((SigmaPrimeVo)/M)
        Ss = (cAlpha*dZ*m.log(longTermTime,10))/5
        if Ss+S1 > dZ:
            S1 = 0
            Ss = 0
        settlements.append(S1)
        longTermSettlements.append(Ss+S1)
        cOld = 1
        CheckC = 0
        count = 0
        while CheckC == 0:
            count += 1
            c = ((((qt-SigmaVo)/(SigmaVo/0.1))**2+(m.log10(Fr)+55.42)**2-272.38)/(275.19-272.38))
            if abs(c-cOld) < 0.01:
                cOld = c
                CheckC = 1
            else:
                cOld = c
            if count > 500:
                print ("Counted Out ",c, cOld)
                CheckC = 1
        
        if c > 1:
            c = 1
        if c < 0.25:
            c = 0.25

        if SigmaPrimeVo == 0:
            SigmaPrimeVo = 0.001
        if (SigmaPrimeVo/0.1)**c < 10**-5:
            SecondPart = 10**-5
        else:
            SecondPart = (SigmaPrimeVo/0.1)**c 

        qt1Net = (qt-SigmaVo)/SecondPart
        wL = 10**(1.506+0.310*m.log10(Fr)-m.log10(qt1Net)/2.526)
        Cc = 0.009*(wL - 10)
        if Cc > 1.6:
            Cc = 1.6
        
        dH = dZ*Cc/(1+0.8)*m.log10((SigmaPrimeVo+deltaSigmaV)/SigmaPrimeVo)
        if dH < 0:
            dH = 0
        if dH > dZ:
            dH = dZ

        print (dH)

    settlementArray = np.array(settlements)
    longSettlementArray = np.array(longTermSettlements)
    qtArray = np.array(qtList)
    depthGraph = -1*np.array(depthList)
    MGraph = np.array(MList)
    SummedShortSettles = []
    SummedLongSettles = []

    thisSettle = 0
    for thisItem in np.arange(0,len(settlements)):
        SummedShortSettles.append(np.sum(settlementArray[thisItem:])*1000)
        SummedLongSettles.append(np.sum(longSettlementArray[thisItem:])*1000)
    SummedShortSettles = np.array(SummedShortSettles)    
    SummedLongSettles = np.array(SummedLongSettles)

    fig = plt.figure()

    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    ax1.plot(qtArray, depthGraph)
    #ax1.set_xlim(0,0.6)
    ax1.set_title("qt (MPa)")

    ax2.plot(MGraph, depthGraph)
    ax2.set_title("M(CPT) (MPa)")

    ax3.plot(SummedShortSettles, depthGraph)
    ax3.plot(SummedLongSettles, depthGraph)
    ax3.set_title("Settlements (mm)")

    TitleText = "Settlements Calculation according to theory of elasticity"
    plt.suptitle(TitleText, fontsize=16)
    plt.savefig(GraphName)
        
    return(np.sum(settlementArray)*1000,np.sum(longSettlementArray)*1000)

