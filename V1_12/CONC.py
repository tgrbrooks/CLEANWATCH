import Iso
from math import *
import numpy as np
vol = np.pi*(51/2)*(pow(13,2)-pow(12.5,2))-(np.pi/2)*pow(13,2) #m^3
den = 2300
mass = 2959812.248671454 #vol*den
#print(mass)
#defPPM = [61, 30, 493]
defPPM = [1.0451, 2.711, 7.26e3]
IsoAct = defPPM
revIsoAct = defPPM
IsoList = Iso.CONC
IType = ['Bq/kg' for i in range(len(IsoList))]
IsoDecay = [Iso.U238,
            Iso.Th232,
            Iso.K40]
IsoEff =   [[0, 0, 0, 0, 0], #U238
            [0, 0, 0, 0],    #Th232
            [0]]             #K40
EffErr = IsoEff


def Activity(PPM):
    IAct = []
    for i in range(len(PPM)):
        #IAct.append(PPM[i]*mass)
        IAct.append((Iso.Lam[i]*Iso.Abs[i])/(Iso.Ms[i]*1e6)*mass*PPM[i])
        print('Activity due to ' + Iso.CONC[i] + ' = %.5e' % IAct[i])
    return IAct


def revActivity(BG, Eff, ratio):
    rIsoAct = [0 for i in range(len(IsoList))]
    for i in range(len(BG)):
        maxbg = max(BG[i])
        x = BG[i].index(maxbg)
        if Eff[i][x] != 0:
            #rIsoAct[i] = maxbg/Eff[i][x]/mass/sqrt(ratio)
            rIsoAct[i] = maxbg/Eff[i][x]/mass*(Iso.Ms[i]*1e6)/(Iso.Lam[i]*Iso.Abs[i])/ratio
        else:
            rIsoAct[i] = 0
    return rIsoAct
