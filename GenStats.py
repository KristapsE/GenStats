# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 15:01:02 2015

@author: ke291
"""

#import NewDP4data
#import Plot
import scipy
from numpy import arange
import pickle
import bisect
import subprocess
import os

#Constants for the standard dp4
meanC = 0.0
meanH = 0.0
stdevC = 2.269372270818724
stdevH = 0.18731058105269952

Conditions = [['', ""],                                             #  0
              ['_BF2', "-B '6-311g(d)' -F b3lyp"],                  #  1
              ['_BF3', "-B '6-31g(d,p)' -F mPW1PW91"],              #  2
              ['_BF4', "-B '6-311g(d)' -F mPW1PW91"],               #  3
              ['_BF5', "-B '6-31g(d,p)' -F wp04"],                  #  4
              ['_BF6', "-B '6-311g(d)' -F wp04"],                   #  5
              ['_BF7', "-B '6-31g(d,p)' -F m062x"],                 #  6
              ['_BF8', "-B '6-311g(d)' -F m062x"],                  #  7
              ['_BF9', "-B 'pcs-1' -F b3lyp"],                      #  8
              ['_BF10', "-B 'pcs-1' -F mPW1PW91"],                  #  9
              ['_BF11', "-B 'pcs-1' -F wp04"],                      #  10
              ['_BF12', "-B 'pcs-1' -F m062x"],                     #  11
              ['_BF1s', "--PM6Opt"],                                #  12
              ['_BF4s', "-B '6-311g(d)' -F mPW1PW91 --PM6Opt"],     #  13
              ['_BF7s', "-B '6-31g(d,p)' -F m062x --PM6Opt"],       #  14
              ['_BF1pm7', "--PM7Opt"],                              #  15
              ['_BF4pm7', "-B '6-311g(d)' -F mPW1PW91 --PM7Opt"],   #  16
              ['_BF7pm7', "-B '6-31g(d,p)' -F m062x --PM7Opt"],     #  17
              ['_BF1h', "--HFOpt"],                                 #  18
              ['_BF4h', "-B '6-311g(d)' -F mPW1PW91 --HFOpt"],      #  19
              ['_BF7h', "-B '6-31g(d,p)' -F m062x --HFOpt"],        #  20
              ['_BF1d', "--DFTOpt"],                                #  21
              ['_BF4d', "-B '6-311g(d)' -F mPW1PW91 --DFTOpt"],     #  22
              ['_BF7d', "-B '6-31g(d,p)' -F m062x --DFTOpt"],       #  23
              ['_BF1m', "--M06Opt"],                                #  24
              ['_BF4m', "-B '6-311g(d)' -F mPW1PW91 --M06Opt"],     #  25
              ['_BF7m', "-B '6-31g(d,p)' -F m062x --M06Opt"]]       #  26

ConTitles = ['b3lyp/6-31g(d,p)','b3lyp/6-311g(d)',
             'mPW1PW91/6-31g(d,p)','mPW1PW91/6-311g(d)',
             'wp04/6-31g(d,p)','wp04/6-311g(d)',
             'm062x/6-31g(d,p)','m062x/6-311g(d)',
             'b3lyp/pcS-1', 'mPW1PW91/pcS-1',
             'wp04/pcS-1', 'm062x/pcS-1']

ESuffixes = ["", "E7", "E5", "E4", "E3", "E2", "E1"]


def RunCalcSet(SetFile, StatsOpt = '', cond = 0, IgnoreIdx=-1, Print=True, Suffix = "", RealData=False):

    infile = open(SetFile, 'r')
    inp = infile.readlines()
    infile.close()
    cwd = os.getcwd()

    Titles = []
    Cdata = []
    Hdata = []
    FlatCdata = []
    FlatHdata = []
    Jdata = []
    CombRs = []
    CRs = []
    HRs = []
    JRs = []
    Positions = []
    TopConfs = []
    MaxCErrors = []
    MaxHErrors = []
    CCorrCoefs = []
    HCorrCoefs = []
    CMAEs = []
    HMAEs = []
    CCMAEs = []
    HCMAEs = []

    for i in range(len(inp)/4):
        if i != IgnoreIdx:
            [folder, correct] = inp[4*i].split(' ')
            cmd1 = inp[4*i+1][:-1]
            cmd2 = inp[4*i+2][:-1]
            Titles.append(cmd2.split(' ')[-1][:-3])
            os.chdir(folder+Conditions[cond][0]+Suffix)
            if Print == True:
                print os.getcwd()
            if StatsOpt == '':
                if Print ==True:
                    print cmd1 + ' ' + Conditions[cond][1] + ' ' + cmd2
                outp = subprocess.check_output(cmd1 + ' ' + Conditions[cond][1] +\
                    ' ' + cmd2, shell=True)
            else:
                if Print == True:
                    print cmd1 + ' ' + StatsOpt + ' ' + Conditions[cond][1] + ' ' + cmd2
                outp = subprocess.check_output(cmd1 + ' ' + StatsOpt + ' ' + \
                    Conditions[cond][1] + ' ' + cmd2, shell=True)
            C, H, J = DP4GetTables(outp, int(correct))
            CombRes, CRes, HRes, JRes = DP4GetRes(outp, int(correct))
            Pos, SigConfs, TopConf, MaxCError, MaxHError = DP4GetMisc(outp, int(correct))
            CCorrCoef, HCorrCoef, CMAE, HMAE, CCMAE, HCMAE = DP4GetSimpleStats(outp, int(correct))
            CCorrCoefs.append(CCorrCoef)
            HCorrCoefs.append(HCorrCoef)
            CMAEs.append(CMAE)
            HMAEs.append(HMAE)
            CCMAEs.append(CCMAE)
            HCMAEs.append(HCMAE)
            Cdata.append(C)
            Hdata.append(H)
            FlatCdata.extend(C)
            FlatHdata.extend(H)
            #Jdata.extend(J)
            Jdata.append(J)
            CombRs.append(CombRes)
            CRs.append(CRes)
            HRs.append(HRes)
            Positions.append(Pos)
            TopConfs.append(TopConf)
            MaxCErrors.append(MaxCError)
            MaxHErrors.append(MaxHError)

    os.chdir(cwd)
    
    #Plot.CScatter(FlatCdata, Hdata)
    print len(FlatCdata)
    #PrintSimpleStats(Titles, CCorrCoefs, HCorrCoefs, CMAEs, HMAEs, CCMAEs, HCMAEs)
    if not RealData:
        return FlatCdata, FlatHdata, Jdata, Titles, CombRs, CRs, HRs, JRs, Positions, TopConfs, MaxCErrors, MaxHErrors
    else:
        return FlatCdata, FlatHdata, Jdata, Titles, CombRs, CRs, HRs, JRs, Positions, TopConfs, MaxCErrors, MaxHErrors,\
               Cdata, Hdata, CCorrCoefs, HCorrCoefs, CMAEs, HMAEs, CCMAEs, HCMAEs


def RunSingleExp(SetFile, ExpIdx, StatsOpt = '', cond = 0, Econd=-1):

    infile = open(SetFile, 'r')
    inp = infile.readlines()
    infile.close()
    cwd = os.getcwd()

    i = ExpIdx

    [folder, correct] = inp[4*i].split(' ')
    cmd1 = inp[4*i+1][:-1]
    cmd2 = inp[4*i+2][:-1]
    Title = cmd2.split(' ')[-1][:-3]
    os.chdir(folder+Conditions[cond][0])
    print os.getcwd()
    if Econd==-1:
        if StatsOpt == '':
            print cmd1 + ' ' + Conditions[cond][1] + ' ' + cmd2
            outp = subprocess.check_output(cmd1 + ' ' + Conditions[cond][1] +\
                ' ' + cmd2, shell=True)
        else:
            print cmd1 + ' ' + StatsOpt + ' ' + Conditions[cond][1] + ' ' + cmd2
            outp = subprocess.check_output(cmd1 + ' ' + StatsOpt + ' ' + \
                Conditions[cond][1] + ' ' + cmd2, shell=True)
    else:
        Eopt = "-E '" + folder+Conditions[Econd][0] + "/'"
        if StatsOpt == '':
            print cmd1 + ' ' + Conditions[cond][1] + ' ' + Eopt + ' ' + cmd2
            outp = subprocess.check_output(cmd1 + ' ' + Conditions[cond][1] +\
                ' ' + Eopt + ' ' + cmd2, shell=True)
        else:
            print cmd1 + ' ' + StatsOpt + ' ' + Conditions[cond][1] + ' '+\
                Eopt + ' ' + cmd2
            outp = subprocess.check_output(cmd1 + ' ' + StatsOpt + ' ' + \
                Conditions[cond][1] + ' ' + Eopt + ' ' + cmd2, shell=True)

    #C, H, J, CombRes, CRes, HRes, JRes, Pos, SigConfs, _ = ReadDP4(outp, int(correct))
    C, H, J = DP4GetTables(outp, int(correct))
    CombRes, CRes, HRes, JRes = DP4GetRes(outp, int(correct))
    Pos, SigConfs, TopConf, MaxCError, MaxHError = DP4GetMisc(outp, int(correct))

    os.chdir(cwd)

    return Title, CombRes, CRes, HRes, JRes, Pos, C, H


#Out of sample testing of compound sets
def OOSTest(SetFile, Cregions = [50, 106, 148], Hregions = [3.0, 5.0, 7.0], conditions = 0, Econditions=-1):

    infile = open(SetFile, 'r')
    inp = infile.readlines()
    infile.close()

    nCompounds = len(inp)/4

    cwd = os.getcwd()
    AllCombRs = [[],[],[],[],[],[]]
    AllPos = [[],[],[],[],[],[]]
    Titles = []

    for i in range(nCompounds):

        #Run the set without the compound in question to get the data for stats
        if Econditions == -1:
            Cdata, Hdata, _, Titles, CombRs, _, _, _, Positions, TopConfs, MaxCErrors,\
            MaxHErrors = RunCalcSet(SetFile, '', conditions, IgnoreIdx=i)
        else:
            Cdata, Hdata,Titles, CombRs, CRs, HRs, Positions, TopConfs, \
                MaxCErrors, MaxHErrors = RunCalcSetEnergy(SetFile, '', conditions, Econditions, IgnoreIdx=i)

        print "Writing models excluding compound " + str(i) + " out of " + str(nCompounds)
        WriteStatsModels(SetFile + '_comp' + str(i), Cdata, Hdata, Cregions, Hregions, DrawPlot=False)
        #WriteKDEfiles(SetFile + '_comp' + str(i), Cdata, Hdata, Cregions, Hregions, DrawPlot=False)
        #GenRMultiGaus(SetFile + str(i) + 'r', 'r')
        GenMultiGaus(SetFile + '_comp' + str(i) + 'k', nGaus=2)

        cwd = os.getcwd()
        StatsOpts = ['',
                '-S g' + cwd + '/' + SetFile + '_comp' + str(i) + 'g',
                 '-S k' + cwd + '/' + SetFile + '_comp' + str(i) + 'k',
                 '-S r' + cwd + '/' + SetFile + '_comp' + str(i) + 'r',
                 '-S u' + cwd + '/' + SetFile + '_comp' + str(i) + 'u',
                 '-S m' + cwd + '/' + SetFile + '_comp' + str(i) + 'm2']
                 #'-S n' + cwd + '/' + SetFile + str(i) + 'm3',
                 #'-S n' + cwd + '/' + SetFile + str(i) + 'm4',
                 #'-S n' + cwd + '/' + SetFile + str(i) + 'm5']

        #Run a single experiment on the missing compound using the stats from the rest
        for j,opt in enumerate(StatsOpts):
            Title, CombRes, CRes, HRes, JRes, Position, C, H = \
                RunSingleExp(SetFile, i, opt, conditions, Econditions)
            AllCombRs[j].append(CombRes)
            AllPos[j].append(Position)

        Titles.append(Title)

    PrintAllStatsRes(Titles, AllCombRs, AllPos, None,
                 columns=['DP4', 'g', 'k', 'r', 'u'])


def OOSSingle(SetFile, Cregions, Hregions, conditions, Econditions, TestID):

    #Run the set without the compound in question to get the data for stats
    print "Initial run for compound " + str(TestID)

    if Econditions == -1:
        Cdata, Hdata, Jdata, _, CombRs, CRs, HRs, JRs, Positions = \
            RunCalcSet(SetFile, '', conditions, IgnoreIdx=TestID, Print=False)
    else:
        Cdata, Hdata,Titles, CombRs, CRs, HRs, Positions, TopConfs, \
            MaxCErrors, MaxHErrors = RunCalcSetEnergy(SetFile, '', conditions,\
                Econditions, IgnoreIdx=TestID, Print=False)

    print "Writing models excluding compound " + str(TestID)
    WriteStatsModels(SetFile + str(TestID), Cdata, Hdata, Cregions, Hregions, DrawPlot=False)
    GenRMultiGaus(SetFile + str(TestID) + 'r', 'r')
    #GenMultiGaus(SetFile + 'k')

    cwd = os.getcwd()
    StatsOpts = ['',
             '-S g' + cwd + '/' + SetFile + str(TestID) + 'g',
             '-S k' + cwd + '/' + SetFile + str(TestID) + 'k',
             '-S r' + cwd + '/' + SetFile + str(TestID) + 'r',
             '-S u' + cwd + '/' + SetFile + str(TestID) + 'u',
             '-S n' + cwd + '/' + SetFile + str(TestID) + 'm2',
             '-S n' + cwd + '/' + SetFile + str(TestID) + 'm3',
             '-S n' + cwd + '/' + SetFile + str(TestID) + 'm4',
             '-S n' + cwd + '/' + SetFile + str(TestID) + 'm5']

    #Run a single experiment on the missing compound using the stats from the rest
    print "Testing models for compounds " + str(TestID)
    CombRs = []
    AllPos = []
    for j,opt in enumerate(StatsOpts):
        Title, CombRes, CRes, HRes, JRes, Position = \
            RunSingleExp(SetFile, TestID, opt, conditions, Econditions)
        CombRs.append(CombRes)
        AllPos.append(Position)

    return CombRs, AllPos, Title


def WriteStatsModels(SetFile, Cdata, Hdata, Cregions, Hregions, DrawPlot=False):

    CShifts = [x[2] for x in Cdata]
    ScaledCErrors = [x[1] - x[2] for x in Cdata]
    UnscaledCErrors = [x[0] - x[2] for x in Cdata]

    HShifts = [x[2] for x in Hdata]
    ScaledHErrors = [x[1] - x[2] for x in Hdata]
    UnscaledHErrors = [x[0] - x[2] for x in Hdata]

    #Generate fitted gaussian for C and H
    ckde = scipy.stats.gaussian_kde(ScaledCErrors)
    Cmean, Cstdev = OptMultiGaus(ckde, [0.0], [1.0], 1, 'C')

    hkde = scipy.stats.gaussian_kde(ScaledHErrors)
    Hmean, Hstdev = OptMultiGaus(hkde, [0.0], [1.0], 1, 'H')
    WriteGausFile(SetFile + 'g', [Cmean[0],Cstdev[0]], [Hmean[0],Hstdev[0]])
    print "Gaussian model generated and written to " + SetFile + 'g'

    WriteKDEFile(SetFile + 'k', ScaledCErrors, ScaledHErrors)
    print "KDE model generated and written to " + SetFile + 'k'
    print "C model size: " + str(len(ScaledCErrors))
    print "H model size: " + str(len(ScaledHErrors))

    RCShifts, RCErrors = PickRegions(CShifts, ScaledCErrors, Cregions)
    RHShifts, RHErrors = PickRegions(HShifts, ScaledHErrors, Hregions)
    WriteRKDEFile(SetFile + 'r', 'r', RCErrors, Cregions, RHErrors, Hregions)
    print "RKDE model generated and written to " + SetFile + 'r'
    print "Regional C model sizes: "
    for r in RCErrors:
        print str(len(r))

    print "Regional H model sizes: "
    for r in RHErrors:
        print str(len(r))

    URCShifts, URCErrors = PickRegions(CShifts, UnscaledCErrors, Cregions)
    URHShifts, URHErrors = PickRegions(HShifts, UnscaledHErrors, Hregions)
    WriteRKDEFile(SetFile + 'u', 'u', URCErrors, Cregions, URHErrors, Hregions)
    print "URKDE model generated and written to " + SetFile + 'u'

    #Plot.CGKDE4Poster(Cmean, Cstdev, ScaledCErrors)
    #Plot.CURKDE4Poster(URCErrors)
    #Plot.HURKDE4Poster(URHErrors)
    #Plot.HGKDE4Poster(Hmean, Hstdev, ScaledHErrors)

    if DrawPlot:
        pass
        #Plot.CScatter(Cdata, Hdata)
        #Plot.AllStats(Cmean, Cstdev, Hmean, Hstdev, ScaledCErrors, ScaledHErrors,
        #      Cregions, RCErrors, URCErrors, Hregions, RHErrors, URHErrors,
        #      CShifts, UnscaledCErrors, HShifts, UnscaledHErrors)


def WriteKDEfiles(SetFile, Cdata, Hdata, Cregions, Hregions, DrawPlot=False):

    CShifts = [x[2] for x in Cdata]
    ScaledCErrors = [x[1] - x[2] for x in Cdata]
    UnscaledCErrors = [x[0] - x[2] for x in Cdata]

    HShifts = [x[2] for x in Hdata]
    ScaledHErrors = [x[1] - x[2] for x in Hdata]
    UnscaledHErrors = [x[0] - x[2] for x in Hdata]

    WriteKDEFile(SetFile + 'k', ScaledCErrors, ScaledHErrors)
    print "KDE model generated and written to " + SetFile + 'k'
    print "C model size: " + str(len(ScaledCErrors))
    print "H model size: " + str(len(ScaledHErrors))

    RCShifts, RCErrors = PickRegions(CShifts, ScaledCErrors, Cregions)
    RHShifts, RHErrors = PickRegions(HShifts, ScaledHErrors, Hregions)
    WriteRKDEFile(SetFile + 'r', 'r', RCErrors, Cregions, RHErrors, Hregions)
    print "RKDE model generated and written to " + SetFile + 'r'
    print "Regional C model sizes: "
    for r in RCErrors:
        print str(len(r))

    print "Regional H model sizes: "
    for r in RHErrors:
        print str(len(r))

    URCShifts, URCErrors = PickRegions(CShifts, UnscaledCErrors, Cregions)
    URHShifts, URHErrors = PickRegions(HShifts, UnscaledHErrors, Hregions)
    WriteRKDEFile(SetFile + 'u', 'u', URCErrors, Cregions, URHErrors, Hregions)
    print "URKDE model generated and written to " + SetFile + 'u'


def GenMultiGaus(kdefile, nGaus=5):

    Cerrors, Herrors = ReadKDEFile(kdefile)

    errors = Cerrors
    kde = scipy.stats.gaussian_kde(errors)
    AllCMeans = []
    AllCStDevs = []
    means = [0.0]
    stdevs = [1.0]
    for i in range(1,nGaus+1):
        means, stdevs = OptMultiGaus(kde, means, stdevs, i, 'C')
        AllCMeans.append(means)
        AllCStDevs.append(stdevs)
    #Plot.PlotMultiGaus(kde, AllCMeans, AllCStDevs, 'C')

    errors = Herrors
    kde = scipy.stats.gaussian_kde(errors)
    AllHMeans = []
    AllHStDevs = []
    means = [0.0]
    stdevs = [1.0]
    for i in range(1,nGaus+1):
        means, stdevs = OptMultiGaus(kde, means, stdevs, i, 'H')
        AllHMeans.append(means)
        AllHStDevs.append(stdevs)
    #Plot.PlotMultiGaus(kde, AllHMeans, AllHStDevs, 'H')

    for Cmeans, Cstdevs, Hmeans, Hstdevs in zip(AllCMeans, AllCStDevs,
                                                AllHMeans, AllHStDevs):
        WriteMultiGausFile(kdefile[:-1] + 'm' + str(len(Cmeans)), Cmeans,
                           Cstdevs, Hmeans, Hstdevs)

#Generates 5 different regional multigaus models - 1,2,3,4,5 gaussian versions
#Number of regions determined by the KDE file
def GenRMultiGaus(rkdefile, t, nGaus=5):

    Cregions, Cerrors, Hregions, Herrors = ReadRKDEFile(rkdefile)

    AllCmeans = [[] for x in range(nGaus)]
    AllCstdevs = [[] for x in range(nGaus)]
    for errors in Cerrors:
        kde = scipy.stats.gaussian_kde(errors)
        means = [0.0]
        stdevs = [1.0]
        for i in range(1,nGaus+1):
            means, stdevs = OptMultiGaus(kde, means, stdevs, i, 'C')
            AllCmeans[i-1].append(means)
            AllCstdevs[i-1].append(stdevs)
    #Plot.PlotMultiGaus(kde, AllCMeans, AllCStDevs, 'C')

    AllHmeans = [[] for x in range(nGaus)]
    AllHstdevs = [[] for x in range(nGaus)]
    for errors in Herrors:
        kde = scipy.stats.gaussian_kde(errors)
        means = [0.0]
        stdevs = [1.0]
        for i in range(1,nGaus+1):
            means, stdevs = OptMultiGaus(kde, means, stdevs, i, 'H')
            AllHmeans[i-1].append(means)
            AllHstdevs[i-1].append(stdevs)
    #Plot.PlotMultiGaus(kde, AllHMeans, AllHStDevs, 'H')

    for Cmeans, Cstdevs, Hmeans, Hstdevs in zip(AllCmeans, AllCstdevs,
                                                AllHmeans, AllHstdevs):
        WriteRMultiGausFile(rkdefile[:-1] + 'm' + str(len(Cmeans)) + str(len(Cmeans[0])), t,
                            Cregions, Cmeans, Cstdevs,
                            Hregions, Hmeans, Hstdevs)


def ReadKDEFile(f):

    infile = open(f, 'r')
    inp = infile.readlines()
    infile.close()

    if 'k' not in inp[0]:
        print "Wrong parameter file type, exiting..."
        quit()

    Cerrors = [float(x) for x in inp[1].split(',')]
    Herrors = [float(x) for x in inp[2].split(',')]
    return Cerrors, Herrors


def ReadRKDEFile(f):

    infile = open(f, 'r')
    inp = infile.readlines()
    infile.close()

    if ('r' not in inp[0]) and ('u' not in inp[0]):
        print "Wrong parameter file type, exiting..."
        quit()

    else:
        buf = inp[1].split(';')
        Cregions = [float(x) for x in buf[0].split(',')]
        Cerrors = [[float(x) for x in y.split(',')] for y in buf[1:]]

        buf = inp[2].split(';')
        Hregions = [float(x) for x in buf[0].split(',')]
        Herrors = [[float(x) for x in y.split(',')] for y in buf[1:]]

        return Cregions, Cerrors, Hregions, Herrors


def GenTestModels(SetFile, Cregions = [50, 106, 148], Hregions = [3.0, 5.0, 7.0],
                  cond = 0, Econd = -1, Suffix = ""):

    AllCombRs = []
    AllCRs = []
    AllHRs = []
    AllPos = []
    AllTopConfs = []

    if Econd == -1:
        Cdata, Hdata, _, Titles, CombRs, CRs, HRs, _, Positions, TopConfs, MaxCErrors,\
            MaxHErrors, Cdata2, Hdata2, CCorrCoefs, HCorrCoefs, CMAEs, HMAEs, CCMAEs, HCMAEs = \
            RunCalcSet(SetFile, '', cond, Suffix=Suffix, RealData=True)
    else:
        Cdata, Hdata,Titles, CombRs, CRs, HRs, Positions, TopConfs, MaxCErrors,\
            MaxHErrors, Cdata2, Hdata2, CCorrCoefs, HCorrCoefs, CMAEs, HMAEs, CCMAEs, HCMAEs = \
            RunCalcSetEnergy(SetFile, '', cond, Econd, Suffix=Suffix, RealData = True)

    AbsCErrors = [abs(x[1] - x[2]) for x in Cdata]
    AbsHErrors = [abs(x[1] - x[2]) for x in Hdata]
    CMAE = sum(AbsCErrors)/len(AbsCErrors)
    HMAE = sum(AbsHErrors)/len(AbsHErrors)
    print "CMAE: " + format(CMAE, "6.2f") + " ppm"
    print "HMAE: " + format(HMAE, "6.2f") + " ppm"
    AllCombRs.append(CombRs)
    AllCRs.append(CRs)
    AllHRs.append(HRs)
    AllPos.append(Positions)
    AllTopConfs.append(TopConfs)
    
    WriteStatsModels(SetFile, Cdata, Hdata, Cregions, Hregions, False)
    #GenRMultiGaus(SetFile + 'r', 'r', nGaus=3)
    #GenRMultiGaus(SetFile + 'u', 'r')   #Generate unscaled multigaus models
    GenMultiGaus(SetFile + 'k', nGaus=2)

    cwd = os.getcwd()
    #StatsOpts = ['-S g' + cwd + '/' + SetFile + 'g']
    """StatsOpts = ['-S g' + cwd + '/' + SetFile + 'g',
                 '-S k' + cwd + '/' + SetFile + 'k',
                 '-S r' + cwd + '/' + SetFile + 'r',
                 '-S u' + cwd + '/' + SetFile + 'u']"""
    """StatsOpts = ['-S g' + cwd + '/' + SetFile + 'g',
                 '-S k' + cwd + '/' + SetFile + 'k',
                 '-S r' + cwd + '/' + SetFile + 'r',
                 '-S u' + cwd + '/' + SetFile + 'u',
                 '-S n' + cwd + '/' + SetFile + 'm'+ str(len(Cregions)+1) + '1',
                 '-S n' + cwd + '/' + SetFile + 'm'+ str(len(Cregions)+1) + '2',
                 '-S n' + cwd + '/' + SetFile + 'm'+ str(len(Cregions)+1) + '3']
                 #'-S n' + cwd + '/' + SetFile + 'm'+ str(len(Cregions)+1) + '4',
                 #'-S n' + cwd + '/' + SetFile + 'm'+ str(len(Cregions)+1) + '5']
    """
    StatsOpts = ['-S g' + cwd + '/' + SetFile + 'g',
                 '-S k' + cwd + '/' + SetFile + 'k',
                 '-S r' + cwd + '/' + SetFile + 'r',
                 '-S u' + cwd + '/' + SetFile + 'u',
                 '-S m' + cwd + '/' + SetFile + 'm2']
                 #'-S m' + cwd + '/' + SetFile + 'm3',
                 #'-S m' + cwd + '/' + SetFile + 'm4',
                 #'-S m' + cwd + '/' + SetFile + 'm5']"""

    for opt in StatsOpts:
        if Econd == -1:
            _, _, _, Titles, CombRs, CRs, HRs, _, Positions, TopConfs, _, _ = \
                RunCalcSet(SetFile, opt, cond, Suffix = Suffix)
        else:
            _, _,Titles, CombRs, CRs, HRs, Positions, TopConfs, _, _ = \
                RunCalcSetEnergy(SetFile, opt, cond, Econd, Suffix=Suffix)

        AllCombRs.append(CombRs)
        AllCRs.append(CRs)
        AllHRs.append(HRs)
        AllPos.append(Positions)
        AllTopConfs.append(TopConfs)

    #res = PrintAllStatsRes(Titles, AllCombRs, AllPos, AllTopConfs,
    #                 columns=['DP4','g','k','r','u','41', '42', '43'])
    PrintAllStatsRes(Titles, AllCombRs, AllPos, AllTopConfs,
                     columns=['DP4', 'g', 'k', 'r','u'])
    print "CMAE: " + format(CMAE, "6.2f") + " ppm"
    print "HMAE: " + format(HMAE, "6.2f") + " ppm"
    #PrintMaxErrors(Titles, MaxCErrors, MaxHErrors)

    DumpCompoundDataToFile(cond, Titles, Cdata2, Hdata2, CCorrCoefs, HCorrCoefs, CMAEs, HMAEs, CCMAEs, HCMAEs,\
                           AllCombRs, AllCRs, AllHRs)
    #CSV format, if possible
    #    Title
    #    Cdata, Hdata
    #    order shifts ascending by experimental shifts - should enforce consistent ordering among different DFT conditions
    #    Rel. probabilities from the various stats models
    #return res


def DumpCompoundDataToFile(cond, Titles, Cdata, Hdata, CCorrCoefs, HCorrCoefs, CMAEs, HMAEs, CCMAEs, HCMAEs,\
                           AllCombRs, AllCRs, AllHRs):

    outfile = open('DP4CompDumpBF' + str(cond+1), 'w', 0)
    stats = ['DP4', 'Gaussian', 'KDE', 'RKDE', 'URKDE', '2gaus']
    for i in range(len(Titles)):
        title = Titles[i]
        C = Cdata[i]
        H = Hdata[i]
        CCorrCoef = CCorrCoefs[i]
        HCorrCoef = HCorrCoefs[i]
        CMAE = CMAEs[i]
        HMAE = HMAEs[i]
        CCMAE = CCMAEs[i]
        HCMAE = HCMAEs[i]
        CombRs = [x[i] for x in AllCombRs]
        CRs = [x[i] for x in AllCRs]
        HRs = [x[i] for x in AllHRs]

        outfile.write(title + ',,BF' + str(cond+1) + ';\n')
        outfile.write(',Exp,Calc,CalcCorrected;\n')
        #def PrintNMR(nucleus, labels, values, scaled, exp):
        #Print("\nAssigned " + nucleus +" shifts: (label, calc, corrected, exp, error)")
        C.sort(key=lambda x: x[2])
        H.sort(key=lambda x: x[2])
        outfile.write('Carbon shifts')
        for shift in C:
            outfile.write(',' + ','.join([format(shift[2], ".2f"),format(shift[0], ".2f"),format(shift[1], ".2f")]) + ';\n')
        outfile.write(';\n')
        outfile.write('Proton shifts')
        for shift in H:
            outfile.write(',' + ','.join([format(shift[2], ".2f"), format(shift[0], ".2f"), format(shift[1], ".2f")]) + ';\n')
        outfile.write(';\n')
        outfile.write('C MAE,,' + format(CMAE, ".2f") + ',' + format(CCMAE, ".2f") + ';\n')
        outfile.write('H MAE,,' + format(HMAE, ".2f") + ',' + format(HCMAE, ".2f") + ';\n')
        outfile.write('C CorrCoef,,' + format(CCorrCoef, ".4f") + ';\n')
        outfile.write('H CorrCoef,,' + format(HCorrCoef, ".4f") + ';\n')
        outfile.write(';\n')
        for CombR, CR, HR, s in zip(CombRs, CRs, HRs, stats):
            outfile.write(s + ' combined,,' + format(CombR, "4.1f") + '%;\n')
            outfile.write(s + ' C only,,' + CR + ';\n')
            outfile.write(s + ' H only,,' + HR + ';\n')
            outfile.write(';\n')
        outfile.write(';\n')
    outfile.close()
    print('NMR data written to DP4CompDumpBF'  + str(cond+1))


def RunCalcSetEnergy(SetFile, StatsOpt = '', cond = 0, Econd = 0, IgnoreIdx=-1, Print=True, Suffix="", RealData=False):

    infile = open(SetFile, 'r')
    inp = infile.readlines()
    infile.close()
    cwd = os.getcwd()

    Titles = []
    Cdata = []
    Hdata = []
    FlatCdata = []
    FlatHdata = []
    CombRs = []
    CRs = []
    HRs = []
    Positions = []
    TopConfs = []
    MaxCErrors = []
    MaxHErrors = []
    CCorrCoefs = []
    HCorrCoefs = []
    CMAEs = []
    HMAEs = []
    CCMAEs = []
    HCMAEs = []

    for i in range(len(inp)/4):
        if i != IgnoreIdx:
            [folder, correct] = inp[4*i].split(' ')
            cmd1 = inp[4*i+1][:-1]
            cmd2 = inp[4*i+2][:-1]
            Titles.append(cmd2.split(' ')[-1][:-3])
            os.chdir(folder+Conditions[cond][0]+Suffix)
            Eopt = "-E '" + folder+Conditions[Econd][0]+Suffix + "/'"
            if Print is True:
                print os.getcwd()
            if StatsOpt == '':
                if Print is True:
                    print cmd1 + ' ' + Conditions[cond][1] + ' '  + Eopt + ' ' + cmd2
                outp = subprocess.check_output(cmd1 + ' ' + Conditions[cond][1] +\
                    ' ' + Eopt + ' ' + cmd2, shell=True)
            else:
                if Print is True:
                    print cmd1 + ' ' + StatsOpt + ' ' + Conditions[cond][1] + ' ' + Eopt + ' ' +cmd2
                outp = subprocess.check_output(cmd1 + ' ' + StatsOpt + ' ' + \
                    Conditions[cond][1] + ' ' + Eopt + ' ' + cmd2, shell=True)
            C, H, J = DP4GetTables(outp, int(correct))
            CombRes, CRes, HRes, JRes = DP4GetRes(outp, int(correct))
            Pos, SigConfs, TopConf, MaxCError, MaxHError = DP4GetMisc(outp, int(correct))
            CCorrCoef, HCorrCoef, CMAE, HMAE, CCMAE, HCMAE = DP4GetSimpleStats(outp, int(correct))
            CCorrCoefs.append(CCorrCoef)
            HCorrCoefs.append(HCorrCoef)
            CMAEs.append(CMAE)
            HMAEs.append(HMAE)
            CCMAEs.append(CCMAE)
            HCMAEs.append(HCMAE)
            Cdata.append(C)
            Hdata.append(H)
            FlatCdata.extend(C)
            FlatHdata.extend(H)
            CombRs.append(CombRes)
            CRs.append(CRes)
            HRs.append(HRes)
            Positions.append(Pos)
            TopConfs.append(TopConf)
            MaxCErrors.append(MaxCError)
            MaxHErrors.append(MaxHError)

    os.chdir(cwd)
    #PrintSimpleStats(Titles, CCorrCoefs, HCorrCoefs, CMAEs, HMAEs, CCMAEs, HCMAEs)
    if not RealData:
        return FlatCdata, FlatHdata, Titles, CombRs, CRs, HRs, Positions, TopConfs, MaxCErrors, MaxHErrors
    else:
        return FlatCdata, FlatHdata, Titles, CombRs, CRs, HRs, Positions, TopConfs, MaxCErrors, MaxHErrors, \
               Cdata, Hdata, CCorrCoefs, HCorrCoefs, CMAEs, HMAEs, CCMAEs, HCMAEs


def GenTestModelsEnergies(SetFile, Cregions = [50, 106, 148], Hregions = [3.0, 5.0, 7.0], conditions = 0, Econditions = 0):

    AllCombRs = []
    AllPos = []

    #Cdata, Hdata = NewDP4data.GetNucleoData()
    Cdata, Hdata, Titles, CombRs, CRs, HRs, Positions, _, _, _ = \
        RunCalcSetEnergy(SetFile, '', conditions, Econditions)
    AllCombRs.append(CombRs)
    AllPos.append(Positions)

    CShifts = [x[2] for x in Cdata]
    ScaledCErrors = [x[1] - x[2] for x in Cdata]
    UnscaledCErrors = [x[0] - x[2] for x in Cdata]

    HShifts = [x[2] for x in Hdata]
    ScaledHErrors = [x[1] - x[2] for x in Hdata]
    UnscaledHErrors = [x[0] - x[2] for x in Hdata]

    #Generate fitted gaussian for C and H
    ckde = scipy.stats.gaussian_kde(ScaledCErrors)
    Cmean, Cstdev = OptMultiGaus(ckde, [0.0], [3.0], 1, 'C')
    #kde, InitMeans, InitStDevs, nGaus, nucleus
    hkde = scipy.stats.gaussian_kde(ScaledHErrors)
    Hmean, Hstdev = OptMultiGaus(hkde, [0.0], [0.3], 1, 'H')
    WriteGausFile(SetFile + 'g', [Cmean[0],Cstdev[0]], [Hmean[0],Hstdev[0]])
    print "Gaussian model generated and written to " + SetFile + 'g'

    WriteKDEFile(SetFile + 'k', ScaledCErrors, ScaledHErrors)
    print "KDE model generated and written to " + SetFile + 'k'

    RCShifts, RCErrors = PickRegions(CShifts, ScaledCErrors, Cregions)
    RHShifts, RHErrors = PickRegions(HShifts, ScaledHErrors, Hregions)
    WriteRKDEFile(SetFile + 'r', 'r', RCErrors, Cregions, RHErrors, Hregions)
    print "RKDE model generated and written to " + SetFile + 'r'

    URCShifts, URCErrors = PickRegions(CShifts, UnscaledCErrors, Cregions)
    URHShifts, URHErrors = PickRegions(HShifts, UnscaledHErrors, Hregions)
    WriteRKDEFile(SetFile + 'u', 'u', URCErrors, Cregions, URHErrors, Hregions)
    print "URKDE model generated and written to " + SetFile + 'u'

    """Plot.AllStats(Cmean, Cstdev, Hmean, Hstdev, ScaledCErrors, ScaledHErrors,
              Cregions, RCErrors, URCErrors, Hregions, RHErrors, URHErrors,
              CShifts, UnscaledCErrors, HShifts, UnscaledHErrors)
    """
    cwd = os.getcwd()
    StatsOpts = ['-S g' + cwd + '/' + SetFile + 'g',
                 '-S k' + cwd + '/' + SetFile + 'k',
                 '-S r' + cwd + '/' + SetFile + 'r',
                 '-S u' + cwd + '/' + SetFile + 'u']
    #StatsOpts = ['-S g' + SetFile + 'g',
    #             '-S k' + SetFile + 'k',
    #             '-S r' + SetFile + 'r',
    #             '-S u' + SetFile + 'u']

    for opt in StatsOpts:
        Cdata, Hdata, Titles, CombRs, CRs, HRs, Positions, _, _, _ = \
            RunCalcSetEnergy(SetFile, opt, conditions, Econditions)
        AllCombRs.append(CombRs)
        AllPos.append(Positions)

    top1s = []
    top2s = []
    for calc in AllPos:
        top1s.append(calc.count(0)/float(len(calc)))
        top2s.append((calc.count(0) + calc.count(1))/float(len(calc)))

    return top1s


def DP4GetTables(Input, CorrectIsomer):

    inp = Input.split('\n')
    CLines = []
    HLines = []
    JLines = []

    for i in range(len(inp)):
        if 'Assigned C shifts:' in inp[i]:
            CLines.append(i)
        if 'Assigned H shifts:' in inp[i]:
            HLines.append(i)
        if 'Direct J comparison table' in inp[i]:
            JLines.append(i)

    TableC = []
    TableH = []
    TableJ = []

    for line in inp[CLines[CorrectIsomer-1]+1:]:
        if ' error:' in line:
            break
        data = filter(None, line.split(' '))
        TableC.append([float(x) for x in data[1:4]])

    for line in inp[HLines[CorrectIsomer-1]+1:]:
        if ' error:' in line or 'J comparison' in line:
            break
        data = filter(None, line.split(' '))
        TableH.append([float(x) for x in data[1:4]])

    if len(JLines) > 1:
        for line in inp[JLines[CorrectIsomer-1]+1:]:
            if ' error:' in line:
                break
            data = filter(None, line.split(' '))
            print data
            TableJ.append([float(x) for x in data])
    else:
        TableJ = []

    return TableC, TableH, TableJ


def DP4GetSimpleStats(Input, CorrectIsomer):

    #Get exp, scaled, unscaled shifts for C and H
    #Calc correl coeff, MAE, CMAE for C and H

    inp = Input.split('\n')
    CLines = []
    HLines = []

    for i in range(len(inp)):
        if 'Assigned C shifts:' in inp[i]:
            CLines.append(i)
        if 'Assigned H shifts:' in inp[i]:
            HLines.append(i)

    TableC = []
    TableH = []

    for line in inp[CLines[CorrectIsomer-1]+1:]:
        if ' error:' in line:
            break
        data = filter(None, line.split(' '))
        TableC.append([float(x) for x in data[1:4]])

    for line in inp[HLines[CorrectIsomer-1]+1:]:
        if ' error:' in line or 'J comparison' in line:
            break
        data = filter(None, line.split(' '))
        TableH.append([float(x) for x in data[1:4]])

    print TableC, TableH
    CAbsoluteErrors = [abs(x[2] - x[0]) for x in TableC]
    CMAE = sum(CAbsoluteErrors)/len(CAbsoluteErrors)
    HAbsoluteErrors = [abs(x[2] - x[0]) for x in TableH]
    HMAE = sum(HAbsoluteErrors) / len(HAbsoluteErrors)

    CCAbsoluteErrors = [abs(x[2] - x[1]) for x in TableC]
    CCMAE = sum(CCAbsoluteErrors)/len(CCAbsoluteErrors)

    HCAbsoluteErrors = [abs(x[2] - x[1]) for x in TableH]
    HCMAE = sum(HCAbsoluteErrors) / len(HCAbsoluteErrors)

    from scipy.stats import pearsonr as pearsonr

    Cexp = [x[2] for x in TableC]
    Ccalc = [x[0] for x in TableC]

    Hexp = [x[2] for x in TableH]
    Hcalc = [x[0] for x in TableH]

    CCorrCoef = pearsonr(Cexp, Ccalc)[0]
    HCorrCoef = pearsonr(Hexp, Hcalc)[0]

    #print "Simple stats: "
    #print CCorrCoef, HCorrCoef, CMAE, HMAE, CCMAE, HCMAE
    return CCorrCoef, HCorrCoef, CMAE, HMAE, CCMAE, HCMAE


def DP4GetRes(Input, CorrectIsomer):

    inp = Input.split('\n')
    JLines = []

    for i in range(len(inp)):
        if 'Direct J comparison table' in inp[i]:
            JLines.append(i)
        if 'all available data' in inp[i]:
            CombIdx = i
        if 'carbon data only' in inp[i]:
            CIdx = i
        if 'proton data only' in inp[i]:
            HIdx = i
        if 'J value data only' in inp[i]:
            JIdx = i

    Results = []
    for line in inp[CombIdx+1:]:
        if 'Isomer ' in line:
            strval = filter(None, line.split(' '))[2]
            #print strval
            Results.append(float(strval.strip('%')))
            if 'Isomer ' + str(CorrectIsomer) + ':' in line:
                CombRes = Results[-1]
        else:
            break

    for line in inp[CIdx+1:]:
        if 'Isomer ' + str(CorrectIsomer) + ':' in line:
            CRes = filter(None, line.split(' '))[2]
            break

    for line in inp[HIdx+1:]:
        if 'Isomer ' + str(CorrectIsomer) + ':' in line:
            HRes = filter(None, line.split(' '))[2]
            break

    if len(JLines) > 1:
        for line in inp[JIdx+1:]:
            if 'Isomer ' + str(CorrectIsomer) + ':' in line:
                JRes = filter(None, line.split(' '))[2]
                break
    else:
        JRes = []

    return CombRes, CRes, HRes, JRes


def DP4GetMisc(Input, CorrectIsomer):

    inp = Input.split('\n')
    MaxCerrs = []
    MaxHerrs = []

    for i in range(len(inp)):
        if 'Max C error:' in inp[i]:
            MaxCerrs.append(i)
        if 'Max H error:' in inp[i]:
            MaxHerrs.append(i)
        if 'all available data' in inp[i]:
            CombIdx = i
        if 'conformers for isomer ' + str(CorrectIsomer) + ':' in inp[i]:
            SigConfs = int(inp[i].split(': ')[1])

    Results = []
    for line in inp[CombIdx+1:]:
        if 'Isomer ' in line:
            strval = filter(None, line.split(' '))[2]
            print strval
            Results.append(float(strval.strip('%')))
            if 'Isomer ' + str(CorrectIsomer) + ':' in line:
                CombRes = Results[-1]
        else:
            break

    print inp[MaxCerrs[CorrectIsomer-1]]
    print inp[MaxHerrs[CorrectIsomer-1]]
    MaxCError = float(inp[MaxCerrs[CorrectIsomer-1]].split(': ')[1])
    MaxHError = float(inp[MaxHerrs[CorrectIsomer-1]].split(': ')[1])
    Pos = sorted(Results, reverse = True).index(CombRes)
    TopConfidence = sorted(Results, reverse = True)[0]

    return Pos, SigConfs, TopConfidence, MaxCError, MaxHError


def PrintAllStatsRes(titles, AllCombs, AllPos, AllTopConfs, columns = ['DP4','g','k','r','u']):

    maxtitle = max([len(x) for x in titles])

    print "\nOverall probabilities given by various models:"
    print "Title".ljust(maxtitle) + " ".join([x.rjust(8) for x in columns])
    
    print len(titles), str(titles)
    print len(AllCombs)
    print ", ".join([str(len(x)) for x in AllCombs])
    print AllCombs
    
    for i, title in enumerate(titles):
        print title.ljust(maxtitle) + " " +\
            " ".join([format(AllCombs[x][i], "4.1f").rjust(8) for x in range(len(AllCombs))])

    print "\nCorrect structure positions given by various models:"
    print "Title".ljust(maxtitle) + " ".join([x.rjust(8) for x in columns])

    for i, title in enumerate(titles):
        print title.ljust(maxtitle) + " " +\
            " ".join([format(AllPos[x][i]+1, "4d").rjust(8) for x in range(len(AllPos))])

    print 40*'-'
    top1s = []
    top2s = []
    AvgConfs = []

    for calc in AllPos:
        top1s.append(calc.count(0)/float(len(calc)))
        top2s.append((calc.count(0) + calc.count(1))/float(len(calc)))
    if AllTopConfs is not None:
        for calc in AllTopConfs:
            AvgConfs.append(sum(calc)/float(len(calc)))

    print '% of 1 '.ljust(maxtitle) +\
        " ".join([format(x*100, "4.1f").rjust(8) for x in top1s])
    print '% of 1 or 2 '.ljust(maxtitle) +\
        " ".join([format(x*100, "4.1f").rjust(8) for x in top2s])
    if AllTopConfs is not None:    
        print 'Avg Conf '.ljust(maxtitle) +\
            " ".join([format(x, "4.1f").rjust(8) for x in AvgConfs])

    return " ".join([format(x*100, "4.1f").rjust(8) for x in top1s])


def PrintSimpleStats(titles, CCorrCoefs, HCorrCoefs, CMAEs, HMAEs, CCMAEs, HCMAEs):

    maxtitle = max([len(x) for x in titles])

    print "\nSimple stats of the NMR prediction:"
    print " ".ljust(maxtitle) + " C Corr".rjust(8) + "H Corr".rjust(8)\
          + "CMAE".rjust(8)+ "HMAE".rjust(8) + "CCMAE".rjust(8)  + "HCMAE".rjust(8)

    for title, CCorr, HCorr, CMAE, HMAE, CCMAE, HCMAE in\
            zip(titles, CCorrCoefs, HCorrCoefs, CMAEs, HMAEs, CCMAEs, HCMAEs):
        print title.ljust(maxtitle) + format(CCorr, "6.4f").rjust(8) + format(HCorr, "6.4f").rjust(8) \
              + format(CMAE, "4.2f").rjust(8) + format(HMAE, "4.2f").rjust(8)\
              + format(CCMAE, "4.2f").rjust(8) + format(HCMAE, "4.2f").rjust(8)


def PrintStatsRes(titles, relCombDP4s, relCs, relHs, Positions):

    maxtitle = max([len(x) for x in titles])

    print "\nOverall results of DP4:"
    print "       Overall     C    H"

    for title, relDP4, pos in zip(titles, relCombDP4s, Positions):
        print title.ljust(maxtitle) + " " + format(relDP4, "4.2f").rjust(8) +\
            str(pos+1).rjust(8)


def WriteMultiGausFile(f, Cmeans, Cstdevs, Hmeans, Hstdevs):
    f = file(f, 'w')
    f.write('m\n')
    f.write(','.join([format(x, "8.6f") for x in Cmeans]) + '\n')
    f.write(','.join([format(x, "8.6f") for x in Cstdevs]) + '\n')
    f.write(','.join([format(x, "8.6f") for x in Hmeans]) + '\n')
    f.write(','.join([format(x, "8.6f") for x in Hstdevs]) + '\n')
    f.close()


def WriteGausFile(f, Cdata, Hdata):
    f = file(f, 'w')
    f.write('g\n')
    f.write(format(Cdata[0], "8.6f") + ',' + format(Cdata[1], "8.6f") + '\n')
    f.write(format(Hdata[0], "8.6f") + ',' + format(Hdata[1], "8.6f") + '\n')
    f.close()


def WriteKDEFile(f, Cdata, Hdata):
    f = file(f, 'w')
    f.write('k\n')
    f.write(','.join([format(x, "4.4f") for x in Cdata]) + '\n')
    f.write(','.join([format(x, "4.4f") for x in Hdata]) + '\n')
    f.close()


def WriteRMultiGausFile(f, t, Cregions, Cmeans, Cstdevs, Hregions, Hmeans, Hstdevs):
    f = file(f, 'w')
    f.write(t+'\n')

    #write cregions
    f.write(','.join([format(x, "4.4f") for x in Cregions]) + '\n')
    #write cmeans
    buf = []
    for region in Cmeans:
        buf.append(','.join([format(x, "4.4f") for x in region]))
    f.write(';'.join(buf) + '\n')
    #write cstdevs
    buf = []
    for region in Cstdevs:
        buf.append(','.join([format(x, "4.4f") for x in region]))
    f.write(';'.join(buf) + '\n')

    #write hstdevs
    f.write(','.join([format(x, "4.4f") for x in Hregions]) + '\n')
    #write hmeans
    buf = []
    for region in Hmeans:
        buf.append(','.join([format(x, "4.4f") for x in region]))
    f.write(';'.join(buf) + '\n')
    #write hstdevs
    buf = []
    for region in Hstdevs:
        buf.append(','.join([format(x, "4.4f") for x in region]))
    f.write(';'.join(buf) + '\n')

    f.close()


def WriteRKDEFile(f, t, Cdata, Cregions, Hdata, Hregions):
    f = file(f, 'w')
    f.write(t+'\n')

    Cbuf = []
    Cbuf.append(','.join([format(x, "4.4f") for x in Cregions]))
    for region in Cdata:
        Cbuf.append(','.join([format(x, "4.4f") for x in region]))

    Hbuf = []
    Hbuf.append(','.join([format(x, "4.4f") for x in Hregions]))
    for region in Hdata:
        Hbuf.append(','.join([format(x, "4.4f") for x in region]))

    f.write(';'.join(Cbuf) + '\n')
    f.write(';'.join(Hbuf) + '\n')
    f.close()


def WriteURKDEFile(f, Cdata, Cregions, Hdata, Hregions):
    f = file(f, 'w')
    f.write('u\n')

    Cbuf = []
    Cbuf.append(','.join([format(x, "4.4f") for x in Cregions]))
    for region in Cdata:
        Cbuf.append(','.join([format(x, "4.4f") for x in region]))

    Hbuf = []
    Hbuf.append(','.join([format(x, "4.4f") for x in Hregions]))
    for region in Hdata:
        Hbuf.append(','.join([format(x, "4.4f") for x in region]))

    f.write(';'.join(Cbuf) + '\n')
    f.write(';'.join(Hbuf) + '\n')
    f.close()


def PickRegions(Shifts, Errors, regions):

    PartShifts = []
    PartErrors = []
    PartShifts.append([s for s in Shifts if s<regions[0]])
    PartErrors.append([e for s,e in zip(Shifts, Errors) if s<regions[0]])

    for i in range(len(regions)-1):
        PartShifts.append([s for s in Shifts
                           if s>regions[i] and s<regions[i+1]])
        PartErrors.append([e for s,e in zip(Shifts, Errors)
                           if s>regions[i] and s<regions[i+1]])

    PartShifts.append([s for s in Shifts if s>regions[-1]])
    PartErrors.append([e for s,e in zip(Shifts, Errors) if s>regions[-1]])

    return PartShifts, PartErrors


def OptMultiGaus(kde, InitMeans, InitStDevs, nGaus, nucleus):

    MeansStdevs = InitMeans + [0.0]*(nGaus - len(InitMeans))
    MeansStdevs += InitStDevs + [1.0]*(nGaus - len(InitMeans))

    #MeansStdevs = [0.0 for x in range(nGaus)] + [1.0 for x in range(nGaus)]
    #MeansStdevs = [-0.0051,0.00194,0.0023, 0.35,0.17625,0.15]
    print MeansStdevs
    #f = lambda w: TautError(Cvalues, ExpCvalues, Hvalues, ExpHvalues, w)
    f = lambda w: kdeMultiGausError(kde, w, nGaus, nucleus)
    res = scipy.optimize.minimize(f, MeansStdevs, method='Nelder-Mead')

    print res.x
    return list(res.x)[:nGaus], list(res.x)[nGaus:]


def kdeMultiGausError(kde, w, nGaus, nucleus):

    error = 0

    if nucleus=='C':
        if any([abs(x)>10 for x in w[:nGaus]]):
            return 1000
        if any([(x>10 or x<0.4) for x in w[nGaus:]]):
            return 1000
        X = arange(-20, 20, 0.1)
    elif nucleus=='H':
        if any([abs(x)>2 for x in w[:nGaus]]):
            return 1000
        if any([(x>2 or x<0.02) for x in w[nGaus:]]):
            return 1000
        X = arange(-2, 2, 0.01)
    elif nucleus=='J':
        if any([abs(x)>2 for x in w[:nGaus]]):
            return 1000
        if any([(x>4 or x<0.02) for x in w[nGaus:]]):
            return 1000
        X = arange(-3, 3, 0.01)

    for x in X:
        error += abs(kde(x)-MultiGaus(w[:nGaus], w[nGaus:],x))

    return error/len(X)


def MultiGaus(means, stdevs, x):

    res = 0
    for i in range(len(means)):
        res += scipy.stats.norm(means[i], stdevs[i]).pdf(x)

    return res/len(means)
