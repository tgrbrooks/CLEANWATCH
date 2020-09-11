import Iso, Eff
import numpy as np
from math import *
den = 2165
vol = np.pi*((pow(18,2)*35.5)-pow(13,2)*25.5)
mass = den*vol
#print(mass)
#defPPM = [10e-3, 220e-3, 750]
defPPM = [0.9, 0.42, 0.156]
IType = ['PPM', 'PPM', 'PPM']
IsoAct = defPPM
revIsoAct = defPPM
IsoList = Iso.ROCK
IsoDecay = [Iso.U238,
            Iso.Th232,
            Iso.K40]
IsoEff = [[0, 0, 0, 0, 0],
          [0, 0, 0, 0],
          [0]]
EffErr = [[0, 0, 0, 0, 0],
          [0, 0, 0, 0],
          [0]]


def Activity(PPM):
    IAct = []
    for i in range(len(PPM)):
        IAct.append((Iso.Lam[i]*Iso.Abs[i])/(Iso.Ms[i]*1e6)*mass*PPM[i])
        print('Activity for ' + Iso.ROCK[i] + ' = %5e' % IAct[i])
    return IAct


def revActivity(BG, Eff, ratio):
    rIsoAct = [0 for i in range(len(IsoList))]
    for i in range(len(BG)):
        maxbg = max(BG[i])
        x = BG[i].index(maxbg)
        if Eff[i][x] != 0:
            rIsoAct[i] = maxbg/Eff[i][x]/mass*(Iso.Ms[i]*1e6)/(Iso.Lam[i]*Iso.Abs[i])/ratio
        else:
            rIsoAct[i] = 0
    return rIsoAct
#defAct = Activity(defPPM)
