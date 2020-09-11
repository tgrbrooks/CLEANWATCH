import Iso, Eff
from math import *
import numpy as np

r = 10026.35e-3
vol = np.pi*pow(r, 3)*2
defPPM = [1e-6] # 1e-6 Bq/kg->Bq/m^3
IsoAct = defPPM
revIsoAct = defPPM
IsoList = Iso.WATER
IType = ['Bq/kg']
IsoDecay = [Iso.Rn222]
IsoEff = [Eff.WATERRn222]
EffErr = [Eff.WATERRn222Err]


def Activity(PPM):
    IAct = []
    for i in range(len(PPM)):
        IAct.append(PPM[i]*vol*1e3)
        print('Activity for ' + Iso.WATER[i] + ' = %5e' % IAct[i])
    return IAct


def revActivity(BG, Eff, ratio):
    rIsoAct = []
    for i in range(len(IsoList)):
        maxbg = max(BG[i])
        x = BG[i].index(maxbg)
        if Eff[i][x] != 0:
            rIsoAct.append(maxbg/Eff[i][x]/(vol*1e3)/sqrt(ratio))
        else:
            revIsoAct.append(0)
    #print(rIsoAct)
    return rIsoAct
#defAct = Activity(defPPM)
