from ast import literal_eval
#Properties
Ms = [3.953e-25, 3.853145e-25, 6.636286e-26] #Mass [kg] [U238, Th232, K40]
Lam = [4.916e-18, 1.57e-18, 1.842e-18] #decay constant [U238, Th232, K40]
Abs = [1, 1, 0.00117] #natural abundance [U238, Th232, K40]
#Component Isotopes
PMT =  ['U238', 'Th232', 'K40']
VETO = ['U238', 'Th232', 'K40']
TANK = ['U238', 'Th232', 'K40', 'Co60', 'Cs137']
CONC = ['U238', 'Th232', 'K40']
#TANK = ['U238', 'Th232', 'K40']
ROCK = ['U238', 'Th232', 'K40', 'Fn']
WATER= ['Rn222', 'Rn']
GD =   ['U238', 'Th232', 'U235', 'U238_l', 'Th232_l', 'U235_l']
#Decay Chains
U238 =  ['Pa234', 'Pb214', 'Bi214', 'Bi210', 'Tl210']
Th232 = ['Ac228', 'Pb212', 'Bi211', 'Tl207']
U235 =  ['Th231', 'Fr223', 'Pb211', 'Bi211', 'Tl207']
K40 =   ['K40']
Co60 =  ['Co60']
Cs137 = ['Cs137']
Rn222 = ['Pb214', 'Bi214', 'Bi210', 'Tl210']
FN =    ['Fn']
RN =    ['Rn']
def setPPM(Iso, PPM):
    output = []
    for i in range(len(Iso)):
        try:
            a = literal_eval(input('Input PPM for ' + Iso[i] + ': '))
            a >= 0
            output.append(a)
            print('PPM for ' + Iso[i] + ' = %.5e' % output[i])
        except:
            output.append(PPM[i])
            print('PPM for ' + Iso[i] + ' set to default value of = %.5e' % output[i])
    return output
def disdef(Iso, val, t):
    for i in range(len(Iso)):
        print(t + ' for ' + Iso[i] + ' set to default value = %.5e' % val[i]) 
def setEff(IsoDecay, Iso, IsoEff, IsoErr):
    IEff = IsoEff
    IErr = IsoErr
    for i in range(len(IsoDecay)):
        print('##########################################')
        print(Iso[i] + ' chain')
        for x in range(len(IsoDecay[i])):
            try:
                a = literal_eval(input('Input Efficiency for ' + IsoDecay[i][x] + ': '))
                a >= 0
                IEff[i][x] = a
            except:
                print('Efficiency for ' + IsoDecay[i][x] + ' set to default value')
            try:
                b = literal_eval(input('Input Error for Efficiency for ' + IsoDecay[i][x] @ ': ')
                b >= 0
                IErr[i][x] = b
                print('Efficiency for ' + IsoDecay[i][x] + ' = %.5e +/- %.5e' % (IEff[i][x], IErr[i][x]))
            except:
                IErr[i][x] = 0
                print('Error for ' + IsoDecay[i][x] + ' set to 0')
                print('Efficiency for ' + IsoDecay[i][x] + ' = %.5e +/- %.5e' % (IEff[i][x], IErr[i][x]))
    return IEff, IErr
def BGrate(Act, Eff, Decay):
    t = 0
    Err = Eff
    for i in range(len(IsoDecay)):
        for x in range(len(IsoDecay[i])):
            if IsoDecay[i][x] == 'Tl210':
                IsoBG[i].append(IsoAct[i]*IsoEff[i][x]*0.002)
            else:
                IsoBG[i].append(IsoAct[i]*IsoEff[i][x])
            Err[i].append(Eff.ErrProp(EffErr[i][x], IsoEff[i][x], IsoBG[i][x]))
            print('BG rate for ' + IsoDecay[i][x] + ' = %.5e +/- %.5e' % (IsoBG[i][x], Err[i][x]))
        t += sum(IsoBG[i])
    return t

