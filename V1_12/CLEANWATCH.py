#import component
import PMT
import VETO
import TANK
import CONC
import ROCK
import WATER
import GD
#imports
import Iso
import Eff
#import Prate
import Nrate
import os
from ast import literal_eval
from math import pow, sqrt, isnan
import copy
import numpy as np
import matplotlib.pyplot as plt
import csv

#from uncertainties import ufloat
#Vars
compList = ['PMT', 'VETO', 'TANK', 'CONC', 'ROCK', 'WATER', 'GD']
#PMT
PMTPPM = PMT.defPPM
PMTAct = [] 
PMTEff = PMT.IsoEff
PMTErr = PMT.EffErr
PMT_Nr = [Nrate.PMTU238,
          Nrate.PMTTh232,
          Nrate.PMTK40]
#VETO
VETOPPM = VETO.defPPM
VETOAct = []
VETOEff = VETO.IsoEff
VETOErr = VETO.EffErr
VETO_Nr = [Nrate.VETOU238,
           Nrate.VETOTh232,
           Nrate.VETOK40]
#TANK
TANKPPM = TANK.defPPM
TANKAct = []
TANKEff = TANK.IsoEff
TANKErr = TANK.EffErr
TANK_Nr = [Nrate.TANKU238,
           Nrate.TANKTh232,
           Nrate.TANKK40,
           [Nrate.TANKSTEEL[0]],
           [Nrate.TANKSTEEL[1]]]
#CONC
CONCPPM = CONC.defPPM
CONCAct = []
CONCEff = CONC.IsoEff
CONCErr = CONC.EffErr
CONC_Nr = [[0, 0, 0, 0, 0],
           [0, 0, 0, 0],
           [0]]
#ROCK
ROCKPPM = ROCK.defPPM
ROCKAct = []
ROCKEff = [[0, 0, 0, 0, 0],
           [0, 0, 0, 0],
           [0]]
ROCKErr = [[0, 0, 0, 0, 0],
           [0, 0, 0, 0],
           [0],
           [0]]
ROCK_Nr = [[0, 0, 0, 0, 0],
           [0, 0, 0, 0],
           [0],
           [0]]
#WATER
WATERPPM = WATER.defPPM
WATERAct = []
WATEREff = WATER.IsoEff
WATERErr = WATER.EffErr
WATER_Nr = [Nrate.WATERRn222]
#GD
GDPPM = GD.defPPM
GDAct = []
GDEff = GD.IsoEff
GDErr = GD.EffErr
GD_Nr = [Nrate.GDU238,
         Nrate.GDTh232,
         Nrate.GDU235]
#         Nrate.GDU238,
#         Nrate.GDTh232,
#         Nrate.GDU235]
#input check
ai = False
ei = False
compAct = []
compEff = []
tot = 0
timeD = 0
dtCut = 0.0001 # prompt-delayed time coincidence = 100 us
dsCut = 0.05 # Distance coincidence of < 2 m removes 95%

fast_neutron_rate = 0.023 #per day
radionucleotide_rate = 0.034 #per day

names = ['PMT', 'VETO', 'TANK', 'CONC', 'ROCK', 'WATER', 'GD']
isotopes = [Iso.PMT, Iso.VETO, Iso.TANK, Iso.CONC, Iso.ROCK, Iso.WATER, Iso.GD]
ppm = [PMTPPM, VETOPPM, TANKPPM, CONCPPM, ROCKPPM, WATERPPM, GDPPM]
activity = [PMTAct, VETOAct, TANKAct, CONCAct, ROCKAct, WATERAct, GDAct]
components = [PMT, VETO, TANK, CONC, ROCK, WATER, GD]
prompt_eff = [PMTEff, VETOEff, TANKEff, CONCEff, ROCKEff, WATEREff, GDEff]
error = [PMTErr, VETOErr, TANKErr, CONCErr, ROCKErr, WATERErr, GDErr]
delayed_eff = [PMT_Nr, VETO_Nr, TANK_Nr, CONC_Nr, ROCK_Nr, WATER_Nr, GD_Nr]

#funcs
def menu(): #menu text
    """
    Displays options
    """
    a = ''
    options = ['a', 'e', 'bgr', 'exit', 'td', 'maxbg', 'cb', 'plt']
    while a.lower() not in options:
        print('##################################################')
        print('CLEANWATCH, V1.11')
        print('Alex Healey, UoS, 2020')
        print('Options: ')
        print('- Input Values for Activity      [a]')
        print('- Input Values for Efficiency    [e]')
        print('- Calculate Background Rate      [bgr]')
        print('- Calculate Time Detection       [td]')
        #print('- Minimum Time Detection         [mintd]')
        print('- Calculate Maximum Background   [maxbg]')
        print('- Cleanliness Budget             [cb]')
        print('- Make Plots                     [plt]')
        print('- Exit software                  [exit]')
        print('##################################################')
        a = str(input('Select an option: '))
        if a.lower() in options and a.lower() != 'exit':
            print('Option selected')
            print('Loading...')
            break
    return a


def inputVal(itype):
    iput = []
    a = input('What components would you like to input values of ' + itype  + ' for? [PMT/VETO/TANK/CONC/ROCK/WATER/GD]  ')
    iput = a.split()
    return iput


def clear():
    """
    Clears output
    """
    ui = ""
    while ui.lower() != 'y' or ui.lower() != 'n':
        ui = input('Do you want to clear the output? [y/n] ')
        if ui.lower() == 'y':
            os.system('clc' if os.name == 'nt' else 'clear')
            break
        if ui.lower() == 'n':
            break


def ActDefault():
    for i, name in enumerate(names):
        if name not in compAct:
            print('##########################################')
            print('Default values for PPM for Iso in ',name)
            Iso.disdef(isotopes[i], ppm[i], 'PPM')
            activity[i] = components[i].Activity(ppm[i])
    return activity
    

def EffDefault():
    for comp_i, name in enumerate(names):
        if name not in compEff:
            print('##########################################')
            print('Default values for Efficiency for Iso in ',name)
            for i in range(len(isotopes[comp_i])):
                print('##########################################')
                print(isotopes[comp_i][i] + ' chain')
                for x in range(len(prompt_eff[comp_i][i])):
                    print('Efficiency for ' + components[comp_i].IsoDecay[i][x] 
                          + ' set to default value of = %.5e \u00B1 %.5e' %
                          (prompt_eff[comp_i][i][x], error[comp_i][i][x]))


def totErr(BGErr):
    a = 0
    for i in range(len(BGErr)):
        for x in range(len(BGErr[i])):
            a += pow(BGErr[i][x], 2)
    a = np.sqrt(a)
    return a


def bgrate():
    totBG_P = 0
    totBG_N = 0
    P = []
    N = []
    for comp_i, name in enumerate(names):
        temp_P = [[0 for x in range(len(prompt_eff[comp_i][i]))] for i in range(len(prompt_eff[comp_i]))]
        temp_N = [[0 for x in range(len(prompt_eff[comp_i][i]))] for i in range(len(prompt_eff[comp_i]))]
        BGr_P = 0
        BGr_N = 0
        BGrErr = 0
        print('##########################################')
        print('BG for ', name)
        for i in range(len(prompt_eff[comp_i])):
            print('##########################################')
            print(isotopes[comp_i][i] + ' chain')
            for x in range(len(prompt_eff[comp_i][i])):
                if components[comp_i].IsoDecay[i][x] == 'Tl210':
                    temp_P[i][x] = (prompt_eff[comp_i][i][x]*activity[comp_i][i]*0.002)
                    temp_N[i][x] = (delayed_eff[comp_i][i][x]*activity[comp_i][i]*0.002)
                else:
                    temp_P[i][x] = (activity[comp_i][i]*prompt_eff[comp_i][i][x])
                    temp_N[i][x] = (activity[comp_i][i]*delayed_eff[comp_i][i][x])
                print('BG due to ' + components[comp_i].IsoDecay[i][x] 
                      + ' = %.5e Hz \u00B1 %.5e' %
                      (temp_P[i][x], error[comp_i][i][x]))
            BGrErr = totErr(error[comp_i])
            print('Total BG for ' + isotopes[comp_i][i] + ' = %.5e Hz'% (sum(temp_P[i])))
        for i in range(len(temp_P)):
            BGr_P += sum(temp_P[i])
        for i in range(len(temp_N)):
            BGr_N += sum(temp_N[i])
        print('##########################################')
        print('Total BG due to ',name,' = %.5e Hz \u00B1 %.5e' % (BGr_P, BGrErr))
        totBG_P += BGr_P
        totBG_N += BGr_N
        P.append(temp_P)
        N.append(temp_N)

    totAcc = totBG_P * totBG_N * dtCut * dsCut * (pow(60,2)*24)
    totBG = totAcc + fast_neutron_rate + radionucleotide_rate
    print('##########################################')
    print('Total Accidental rate = %.5e /day' % totAcc)
    print('Total BG rate including FN + RN = %.5e /day' % totBG)
    print('Total prompt BG rate = %.5e Hz' % totBG_P)
    return totBG_P, totAcc, P, N


def tdcalc(totAcc):
    #get signal rate
    try:
        signal = literal_eval(input('Input signal rate: '))
        signal > 0
        print('Signal rate set to = %.5e' % signal)
    except:
        signal = 0.485 #0.387
        print('Signal set to default value = %.5e' % signal)
    #calculate td
    S = signal*0.9
    B = totAcc + fast_neutron_rate + radionucleotide_rate + S*1.15
    sigma = 4.65
    R_onoff = 1.5
    td = pow(sigma, 2)*((B+(B+S)/R_onoff))*(1/pow(S,2))*(1+R_onoff) #/((60**2)*24) #[days] 
    #convert to days
    #td /= (pow(60,2)*24)
    print('Reactor off time to detection @ 3 sigma rate = %.5e' % td + ' days')
    return td


def maxBG():
    try:
        signal = literal_eval(input('Input signal rate: '))
        signal > 0
        print('Signal rate set to = %.5e' % signal)
    except:
        signal = 0.485 #0.387
        print('Signal set to default value = %.5e' % signal)
    try:
        days = literal_eval(input('Input time to dection in days: '))
        days > 0
        print('Time to detection set to = %.5e' % days)
    except:
        days = 390 #1.92437e+02
        print('Time to detection set to default value of = %.5e' % days)
    sigma = 4.65
    S = signal*0.9
    R_onoff = 1.5
    #B = (1.5*days*pow(S, 2))/(2.5*pow(sigma, 2)) - S/2.5
    MBG = np.power(S, 2.)*days*R_onoff/(np.power(sigma, 2.)*np.power(1+R_onoff, 2.)) - S/(R_onoff+1.) 
    MBG = MBG - (S*1.15) - fast_neutron_rate - radionucleotide_rate
    print('Maximum BG rate for time detection of %.5e days = %.5e' % (days, MBG))
    return MBG


def share(maxBG, totAcc, IsoBG_P, CompIso, scales=None):
    IsoShare = IsoBG_P
    for i in range(len(IsoBG_P)):
        scale = 1
        if scales is not None:
            scale = scales[i]
        for x in range(len(IsoBG_P[i])):
            IsoShare[i][x] = IsoBG_P[i][x] + (maxBG-totAcc)/totAcc * scale * IsoBG_P[i][x] # prompt rate for given isotope in Hz
            print('BG due to ' + CompIso[i][x] + '  = %.5e' % IsoShare[i][x])
    return IsoShare # Hz


def CBOUT(IsoAct, BGIsoCB, Iso): #BGIso, Iso): #(, , ,COMP.IsoList)
    for i in range(len(IsoAct)):
        print('Singles Budget for %.7s = %.5e Hz' % (Iso[i], sum(BGIsoCB[i]))) #TODO change this to give the singles rate for prompt events
    # FIXME doesn't make sense to have accidentals budget per isotope
    for i in range(len(IsoAct)):
        print('Accidentals Budget for %.7s = %.5e Hz' % (Iso[i], (sum(BGIsoCB[i]) * dsCut * dtCut)))
    for i in range(len(IsoAct)):
        print('Radioactivity Budget for %.7s = %.5e' % (Iso[i], IsoAct[i]))
        #print('Nominal singles rate for %.7s = %.5e Hz' % (Iso[i], sum(BGIso[i]))) #????


def CBTable(components, prompt_eff, ratio, CB_BG):
    with open('cleanliness_results.csv', 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(['Isotope', 'Radioactivity [Bq]', 'Singles Rate [Hz]', 'Singles Budget [Hz]', 'Radioactivity Budget'])
        # Loop over components
        for i, name in enumerate(names):
            writer.writerow([name, '', '', '', ''])
            # Determine the radioactivity budget assuming it's shared equally amongst components and isotopes
            CB_Act = components[i].revActivity(CB_BG[i], prompt_eff[i], ratio)
            # Loop over isotopes in component
            for j in range(len(components[i].IsoList)):
                activity_exponent = np.floor(np.log10(np.abs(activity[i][j]))).astype(int)
                activity_string = ('%1.2f$\\times10^{%i}$' % (activity[i][j]/(1.*10.**activity_exponent), activity_exponent))
                if sum(CB_BG[i][j]) == 0:
                    writer.writerow([components[i].IsoList[j], activity_string, '0', '0', '0'])
                    continue
                rate = sum(CB_BG[i][j])/ratio
                rate_exponent = np.floor(np.log10(np.abs(rate))).astype(int)
                rate_string = ('%1.2f$\\times10^{%i}$' % (rate/(1.*10.**rate_exponent), rate_exponent))
                sb = (sum(CB_BG[i][j])/sqrt(ratio))
                sb_exponent = np.floor(np.log10(np.abs(sb))).astype(int)
                sb_string = ('%1.2f$\\times10^{%i}$' % (sb/(1.*10.**sb_exponent), sb_exponent))
                rb_exponent = np.floor(np.log10(np.abs(CB_Act[j]))).astype(int)
                rb_string = ('%1.2f$\\times10^{%i}$ %s' % (CB_Act[j]/(1.*10.**rb_exponent), rb_exponent, components[i].IType[j]))
                writer.writerow([components[i].IsoList[j],
                                 activity_string,
                                 rate_string,
                                 sb_string,
                                 rb_string])


def PNRates(Eff, Iso, IsoDecay, Nr, Act):
    BG_P = 0
    BG_N = 0
    N = [[0 for x in range(len(Eff[i]))] for i in range(len(Eff))]
    P = [[0 for x in range(len(Eff[i]))] for i in range(len(Eff))]
    for i in range(len(Iso)):
        for x in range(len(IsoDecay[i])):
            if IsoDecay[i][x] == 'Tl210':
                P[i][x] = (Eff[i][x]*Act[i]*0.002)
                N[i][x] = (Nr[i][x]*Act[i]*0.002)
            else:
                P[i][x] = (Act[i]*Eff[i][x])
                N[i][x] = (Act[i]*Nr[i][x])
    for i in range(len(P)):
        BG_P += sum(P[i])
    for i in range(len(N)):
        BG_N += sum(N[i])
    return BG_P, BG_N


def Background(signal=0.485, act=None):

    totBG_P = 0
    totBG_N = 0
    global activity

    if act is not None:
        activity = act

    ## Calculate the background prompt and delayed delayed_eff
    for i, name in enumerate(names):
        temp_P, temp_N = PNRates(prompt_eff[i], isotopes[i],
                                 components[i].IsoDecay, delayed_eff[i], activity[i])
        totBG_P += temp_P
        totBG_N += temp_N
    
    S = signal*0.9

    # Total accidental rate
    totAcc = totBG_P * totBG_N * dtCut * dsCut * (pow(60,2)*24)
    # + fast neutrons + radionucleotides
    totBG = totAcc + fast_neutron_rate + radionucleotide_rate
    # + other core + 15% other neutrinos
    background = totBG + S*1.15
    #return background
    return totAcc


def BackgroundRatio(signal, days, sigma, act=None):

    background = Background(signal, act)

    # Maximum allowed background
    S = signal*0.9
    R_onoff = 1.5
    MBG = np.power(S, 2.)*days*R_onoff/(np.power(sigma, 2.)*np.power(1+R_onoff, 2.)) - S/(R_onoff+1.)
    MBG = MBG - S*1.15 - fast_neutron_rate - radionucleotide_rate

    return background/MBG


def DwellTime(signal, sigma):

    # Maximum allowed background
    S = signal*0.9
    B = Background(signal) + S*1.15 + fast_neutron_rate + radionucleotide_rate
    R_onoff = 1.5
    D = np.power(sigma, 2.)*((S + B)/R_onoff + B)*(1. + R_onoff)/np.power(S, 2.)
    #D = np.power(sigma, 2.)*(B*(R_onoff+1.) + S*R_onoff)/((1. + R_onoff)*np.power(S, 2.))
    return D


def InvGradients():
    gradients = []
    norm = 0
    activity = [PMTAct, VETOAct, TANKAct, CONCAct, ROCKAct, WATERAct, GDAct]
    if ai == False:
        activity = ActDefault()
    orig_activities = copy.deepcopy(activity)
    for comp_i, name in enumerate(names):
        grad = []
        for i in range(len(isotopes[comp_i])):
            activity = copy.deepcopy(orig_activities)
            activity[comp_i][i] = orig_activities[comp_i][i]*0.5
            print(activity[comp_i][i])
            y1 = BackgroundRatio(0.485, 400, 4.65, activity)
            activity[comp_i][i] = orig_activities[comp_i][i]*1.5
            print(activity[comp_i][i])
            y2 = BackgroundRatio(0.485, 400, 4.65, activity)
            print(y1, y2)
            if (y2-y1) == 0:
                grad.append(0)
                norm += 0
            else:
                grad.append(1./(y2-y1))
                norm += 1./(y2-y1)
        gradients.append(grad)
    activity = orig_activities
    return gradients, norm

######################################################################################################
ans = menu()
while ans.lower() != 'exit':
    if ans.lower() == 'a':
        ai = True
        #get list of compoents
        compAct = inputVal('PPM')
        #change PPM values
        for i in range(len(compAct)):
            for comp_i, name in enumerate(names):
                if compAct[i].upper() == name:
                    print('##########################################')
                    print('Input values for PPM for Iso in ', name)
                    ppm[comp_i] = Iso.setPPM(isotopes[comp_i], components[comp_i].defPPM, components[comp_i].IType)
                    activity[comp_i] = components[comp_i].Activity(ppm[comp_i])

        #set to default
        ActDefault()
        #clear()
        ans = menu()

    if ans.lower() == 'e':
        ei = True
        compEff = inputVal('Efficiency')
        for i in range(len(compEff)):
            print(compEff[i].upper())
            for comp_i, name in enumerate(names):
                if compEff[i].upper() == name:
                    print('##########################################')
                    print('Input values for Efficiency for Iso in ', name)
                    prompt_eff[comp_i], error[comp_i] = Iso.setEff(components[comp_i].IsoDecay, isotopes[comp_i], components[comp_i].IsoEff, components[comp_i].EffErr)

        #set to default
        EffDefault()
        #clear()
        ans = menu()

    if ans.lower() == 'bgr':
        if ai == False:
            activity = ActDefault()
        if ei == False:
            EffDefault()
        totBG_P, totAcc, P, N = bgrate()
        #clear()
        ans = menu()

    if ans.lower() == 'td':
        if ai == False:
            activity = ActDefault()
        if ei == False:
            EffDefault()
        totBG_P, totAcc, P, N = bgrate()
        timeD = tdcalc(totAcc) 
        #clear()
        ans = menu()

    if ans.lower() == 'maxbg':
        i = maxBG()
        #clear()
        ans = menu()

    if ans.lower() == 'cb':
        dist_type = ''
        options = ['e', 'c']
        while dist_type.lower() not in options:
            print('##################################################')
            print('Choose how to distribute the budget: ')
            print('- Equally amongst isotopes       [e]')
            print('- Proportional to 1/contribution [c]')
            print('##################################################')
            dist_type = str(input('Select an option: '))
        #check if activity has been changed
        if ai == False:
            activity = ActDefault()
        #check if prompt_eff has been changed
        if ei == False:
            EffDefault()
        #calculate BG for comps
        totBG_P, totAcc, P, N = bgrate()
        #calculate max BG for signal rate and td
        MBG = maxBG()
        #calculate the shares
        CB_BG = []
        grads = None
        if dist_type.lower() == 'c':
            grads, norm = InvGradients()
            sumg = 0
            for gr in grads:
                for g in gr:
                    sumg += g/norm
            print('norm = ',norm)
            print('sum = ',sumg)
        for i, name in enumerate(names):
            print('##########################################')
            print('Calculated BG for ', name)
            scales = None
            if dist_type.lower() == 'c':
                scales = [x/norm for x in grads[i]]
            temp_CB_BG = share(MBG, totAcc, P[i], components[i].IsoDecay, scales)
            CB_BG.append(temp_CB_BG)
        for i, name in enumerate(names):
            print('##########################################')
            print('CB for ', name)
            CB_Act = components[i].revActivity(CB_BG[i], prompt_eff[i], MBG/totAcc)
            CBOUT(CB_Act, CB_BG[i], components[i].IsoList)
        CBTable(components, prompt_eff, MBG/totAcc, CB_BG)
        #reset
        #clear()
        ans = menu()

    if ans.lower() == 'plt':

        # Get the plot type
        plt_type = ''
        options = ['b', 'd', 'c', 'exit']
        while plt_type.lower() not in options:
            print('##################################################')
            print('Plot type: ')
            print('- Budget vs isotope              [b]')
            print('- Dwell time vs isotope          [d]')
            print('- Contribution breakdown         [c]')
            print('- Exit software                  [exit]')
            print('##################################################')
            plt_type = str(input('Select an option: '))
        if plt_type.lower() == 'exit':
            break

        #check if activity has been changed
        activity = [PMTAct, VETOAct, TANKAct, CONCAct, ROCKAct, WATERAct, GDAct]
        if ai == False:
            activity = ActDefault()

        #check if prompt_eff has been changed
        if ei == False:
            EffDefault()

        # Get the isotope and component to plot in
        component_input = ''
        comp_options = ['pmt', 'veto', 'tank', 'conc', 'rock', 'water', 'gd', 'all']
        while component_input.lower() not in comp_options:
            print('##################################################')
            print('Select a component to study: ')
            print(comp_options)
            component_input = str(input('Select an option: '))

        iso_input = ''
        iso_options = Iso.GetIsotopes(component_input)
        print('##################################################')
        print('Select isotope(s): ')
        print(iso_options)
        iso_input = str(input('Select an option: ')).split(',')
        if iso_input[0].lower() == 'all':
            iso_input = iso_options
        while not set(iso_input).issubset(iso_options):
            print('##################################################')
            print('Select isotope(s): ')
            print(iso_options)
            iso_input = str(input('Select an option: ')).split(',')
            if iso_input[0].lower() == 'all':
                iso_input = iso_options

        if plt_type.lower() == 'c':
            data = []
            labels = []
            total = 0
            for n_i, name in enumerate(names):
                if not (component_input == 'all' or component_input == name.lower()):
                    continue
                for i, iso in enumerate(isotopes[n_i]):
                    if iso not in iso_input:
                        continue
                    P = 0;
                    N = 0;
                    for x in range(len(components[n_i].IsoDecay[i])):
                        if components[n_i].IsoDecay[i][x] == 'Tl210':
                            P += (prompt_eff[n_i][i][x]*activity[n_i][i]*0.002)
                            N += (delayed_eff[n_i][i][x]*activity[n_i][i]*0.002)
                        else:
                            P += (activity[n_i][i]*prompt_eff[n_i][i][x])
                            N += (activity[n_i][i]*delayed_eff[n_i][i][x])
                    BG = P * N * dtCut * dsCut * (pow(60,2)*24)
                    total += BG
                    if BG != 0:
                        data.append(BG)
                        exponent = np.floor(np.log10(np.abs(BG))).astype(int)
                        BG_lab = BG/(1.*10.**exponent)
                        labels.append('%s:%s %1.2f$\\times10^{%i}$' % (name, isotopes[n_i][i], BG_lab, exponent))
            data = np.array(data)
            data = data/total

            fig1, ax1 = plt.subplots()
            wedges, texts, autotexts = ax1.pie(data, autopct='%1.1f%%', startangle=90)
            ax1.legend(wedges, labels, title='Rate (per day)', loc='center left')
            ax1.axis('equal')
            plt.tight_layout()
            plt.show()
            ans = menu()
            continue

        # Get the range of isotope concentration
        min_percent = -1
        while min_percent < 0:
            print('##################################################')
            print('Minimum value as percentage of default: ')
            try:
                min_percent = input('Value: ')
                min_percent = int(min_percent)
            except ValueError:
                min_percent = -1

        max_percent = -1
        while max_percent < 0:
            print('##################################################')
            print('Maximum value as percentage of default: ')
            try:
                max_percent = input('Value: ')
                max_percent = int(max_percent)
            except ValueError:
                max_percent = -1

        signal = 0.485 #0.387
        try:
            signal = literal_eval(input('Input signal rate: '))
            signal > 0
            print('Signal rate set to = %.5e' % signal)
        except:
            print('Signal set to default value = %.5e' % signal)

        days = 390 #1.92437e+02
        if plt_type == 'b':
            try:
                days = literal_eval(input('Input time to dection in days: '))
                days > 0
                print('Time to detection set to = %.5e' % days)
            except:
                print('Time to detection set to default value of = %.5e' % days)

        sigma = 4.65 #3
        try:
            sigma = literal_eval(input('Input significance of detection: '))
            sigma > 0
            print('Significance of detection set to = %.5e sigma' % sigma)
        except:
            print('Significance of detection set to default value of = %.5e sigma' % sigma)

        # Loop over the range and fill data for plotting
        x_data_list = []
        y_data_list = []
        y_label = 'Dwell time [days]'
        orig_activities = copy.deepcopy(activity)
        for isotope in iso_input:
            x_data = []
            y_data = []
            activity = copy.deepcopy(orig_activities)
            for percent in range(min_percent, max_percent, 1):
                x_data.append(percent)
                # Modify the activity
                if component_input != 'all':
                    comp_index = comp_options.index(component_input)
                    iso_index = iso_options.index(isotope)
                    activity[comp_index][iso_index] = orig_activities[comp_index][iso_index]*percent/100.
                else:
                    for comp in comp_options:
                        comp_index = comp_options.index(comp)
                        found = True
                        try:
                            iso_index = Iso.GetIsotopes(comp).index(isotope)
                            activity[comp_index][iso_index] = orig_activities[comp_index][iso_index]*percent/100.
                        except:
                            pass
                if plt_type == 'b':
                    y_data.append(BackgroundRatio(signal, days, sigma))
                    y_label = 'total/maximum accidental background'
                else:
                    y_data.append(DwellTime(signal, sigma))
            x_data_list.append(x_data)
            y_data_list.append(y_data)

        plt.cla()
        for i, isotope in enumerate(iso_input):
            plt.plot(np.array(x_data_list[i]), np.array(y_data_list[i]), label=isotope)
        plt.xlabel('% of default activity in ' + component_input)
        plt.ylabel(y_label)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        activity = orig_activities

        #reset
        ans = menu()
