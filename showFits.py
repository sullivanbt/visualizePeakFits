# This box deals with importing and defines a few functions.  Do not change.
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from mantid.simpleapi import *
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import ICCFitTools as ICCFT
import BVGFitTools as BVGFT
import pySlice
import pickle

# These are parameters for fitting.
qLow = -25  # Lowest value of q for ConvertToMD 
qHigh = 25; # Highest value of q for ConvertToMD
Q3DFrame='Q_lab' # Either 'Q_lab' or 'Q_sample'; Q_lab recommended if using a strong peaks
                 # profile library from a different sample
eventFileName = '/SNS/TOPAZ/IPTS-18474/data/TOPAZ_26751_event.nxs' #Full path to the event nexus file
peaksFile = '/SNS/TOPAZ/IPTS-18474/shared/SC100K_useDetCal/26751_Niggli.integrate' #Full path to the ISAW peaks file
UBFile = '/SNS/TOPAZ/IPTS-18474/shared/SC100K_useDetCal/26751_Niggli.mat' #Full path to the ISAW UB file
#strongPeakParamsFile = '/SNS/MANDI/shared/ProfileFitting/strongPeakParams_beta_lac_mut_mbvg.pkl' #Full path to pkl file
strongPeakParamsFile = None
moderatorCoefficientsFile = '/SNS/MANDI/shared/ProfileFitting/franz_coefficients_2017.dat' #Full path to pkl file
IntensityCutoff = -500 # Minimum number of counts to not force a profile
EdgeCutoff = 10 # Pixels within EdgeCutoff from a detector edge will be have a profile forced. Currently for Anger cameras only.
FracStop = 0.05 # Fraction of max counts to include in peak selection.
MinpplFrac = 0.9 # Min fraction of predicted background level to check
MaxpplFrac = 1.1 # Max fraction of predicted background level to check
DQMax = 0.25 # Largest total side length (in Angstrom-1) to consider for profile fitting.
plotResults = True #Show BVG and ICC Fits separately.

#=================================================================================
# Do not edit below this line
#=================================================================================

def addInstrumentParameters(peaks_ws):
    """
    This function adds parameters to instrument files.  This is only done as a TEMPORARY workaround 
    until the instrument files with these parameters are included in the stable release of mantid 
    which is available on the analysis jupyter server.
    """
    instrumentName = peaks_ws.getInstrument().getName()
    if instrumentName == 'MANDI':
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='fitConvolvedPeak', ParameterType='Bool', Value='False')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='sigX0Scale', ParameterType='Number', Value='1.0')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='sigY0Scale', ParameterType='Number', Value='1.0')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='numDetRows', ParameterType='Number', Value='255')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='numDetCols', ParameterType='Number', Value='255')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='numBinsTheta', ParameterType='Number', Value='50')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='numBinsPhi', ParameterType='Number', Value='50')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='fracHKL', ParameterType='Number', Value='0.4')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='dQPixel', ParameterType='Number', Value='0.003')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='mindtBinWidth', ParameterType='Number', Value='15.0')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='maxdtBinWidth', ParameterType='Number', Value='50.0')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='peakMaskSize', ParameterType='Number', Value='5')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='iccKConv', ParameterType='String', Value='100.0 140.0 120.0')

    elif instrumentName == 'TOPAZ':
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='fitConvolvedPeak', ParameterType='Bool', Value='False')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='sigX0Scale', ParameterType='Number', Value='3.0')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='sigY0Scale', ParameterType='Number', Value='3.0')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='numDetRows', ParameterType='Number', Value='255')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='numDetCols', ParameterType='Number', Value='255')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='numBinsTheta', ParameterType='Number', Value='50')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='numBinsPhi', ParameterType='Number', Value='50')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='fracHKL', ParameterType='Number', Value='0.4')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='dQPixel', ParameterType='Number', Value='0.01')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='mindtBinWidth', ParameterType='Number', Value='2.0')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='maxdtBinWidth', ParameterType='Number', Value='15.0')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='peakMaskSize', ParameterType='Number', Value='15')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='iccB', ParameterType='String', Value='0.001 0.3 0.005')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='iccKConv', ParameterType='String', Value='10.0 1000.0 500.0')

    elif instrumentName == 'CORELLI':
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='fitConvolvedPeak', ParameterType='Bool', Value='True')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='sigX0Scale', ParameterType='Number', Value='2.0')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='sigY0Scale', ParameterType='Number', Value='2.0')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='numDetRows', ParameterType='Number', Value='255')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='numDetCols', ParameterType='Number', Value='16')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='numBinsTheta', ParameterType='Number', Value='50')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='numBinsPhi', ParameterType='Number', Value='50')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='fracHKL', ParameterType='Number', Value='0.4')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='dQPixel', ParameterType='Number', Value='0.007')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='mindtBinWidth', ParameterType='Number', Value='2.0')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='maxdtBinWidth', ParameterType='Number', Value='60.0')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='peakMaskSize', ParameterType='Number', Value='10')
        #SetInstrumentParameter(Workspace='peaks_ws', ParameterName='iccA', ParameterType='String', Value='0.25 0.75 0.5')
        #SetInstrumentParameter(Workspace='peaks_ws', ParameterName='iccB', ParameterType='String', Value='0.001 0.3 0.005')
        #SetInstrumentParameter(Workspace='peaks_ws', ParameterName='iccR', ParameterType='String', Value='0.05 1. 0.1')
        #SetInstrumentParameter(Workspace='peaks_ws', ParameterName='iccKConv', ParameterType='String', Value='10.0 800.0 100.0')



print("Which peak?")
peakNumber = input()

#Load peaks file if not already loaded
importPeaks = True
for ws in mtd.getObjectNames():
    if mtd[ws].getComment() == '%s'%peaksFile:
        print('    using already loaded peaks file')
        importPeaks = False
        peaks_ws = mtd[ws]
if importPeaks:
    peaks_ws = LoadIsawPeaks(Filename = peaksFile)
    peaks_ws.setComment(peaksFile)

# Load MDdata if not already loaded
importFlag=True
for ws in mtd.getObjectNames():
    if mtd[ws].getComment() == 'MD_%i'%peaks_ws.getPeak(peakNumber).getRunNumber():
        print('    using already loaded MDdata')
        MDdata = mtd[ws]
        importFlag = False
        break
if importFlag:
    event_ws = Load(Filename=eventFileName, OutputWorkspace='event_ws')
    MDdata = ConvertToMD(InputWorkspace='event_ws',  QDimensions='Q3D', dEAnalysisMode='Elastic',
                             Q3DFrames=Q3DFrame, QConversionScales='Q in A^-1',
                             MinValues='%f, %f, %f' % (qLow, qLow, qLow), Maxvalues='%f, %f, %f' % (qHigh, qHigh, qHigh), MaxRecursionDepth=10,
                             LorentzCorrection=False)
    MDdata.setComment('MD_%i'%peaks_ws.getPeak(peakNumber).getRunNumber())
addInstrumentParameters(peaks_ws)

LoadIsawUB(InputWorkspace=peaks_ws, Filename=UBFile)
UBMatrix = peaks_ws.sample().getOrientedLattice().getUB()
dQ = np.abs(ICCFT.getDQFracHKL(UBMatrix, frac=0.5))
dQ[dQ>DQMax]=DQMax

if strongPeakParamsFile is None:
    strongPeakParams = None
else:
    strongPeakParams = pickle.load(open(strongPeakParamsFile, 'rb'))

padeCoefficients = ICCFT.getModeratorCoefficients(moderatorCoefficientsFile)
NTheta = peaks_ws.getInstrument().getIntParameter("numBinsTheta")[0]
NPhi = peaks_ws.getInstrument().getIntParameter("numBinsPhi")[0]
MindtBinWidth = peaks_ws.getInstrument().getNumberParameter("minDTBinWidth")[0]
MaxdtBinWidth = peaks_ws.getInstrument().getNumberParameter("maxDTBinWidth")[0]
FracHKL = 0.4 # Fraction of HKL to consider for profile fitting.
DQPixel = peaks_ws.getInstrument().getNumberParameter("DQPixel")[0]
np.warnings.filterwarnings('ignore') # There can be a lot of warnings for bad solutions.

q_frame = 'lab'
#Get some peak variables
peak = peaks_ws.getPeak(peakNumber)
Box = ICCFT.getBoxFracHKL(peak, peaks_ws, MDdata, UBMatrix, peakNumber, dQ, fracHKL=0.5, dQPixel=DQPixel,  q_frame=q_frame)
box = Box
#Set up our filters
qMask = ICCFT.getHKLMask(UBMatrix, frac=FracHKL, dQPixel=DQPixel, dQ=dQ)
n_events = Box.getNumEventsArray()
qMask = ICCFT.getHKLMask(UBMatrix, frac=0.4, dQPixel=DQPixel, dQ=dQ)

iccFitDict = ICCFT.parseConstraints(peaks_ws)
Y3D, gIDX2, pp_lambda2, params2 = BVGFT.get3DPeak(peak, peaks_ws, box, padeCoefficients,qMask,nTheta=NTheta, nPhi=NPhi,
                                               plotResults=plotResults,
                                               zBG=1.96,fracBoxToHistogram=1.0,bgPolyOrder=1, strongPeakParams=strongPeakParams,
                                               q_frame=q_frame, mindtBinWidth=MindtBinWidth, maxdtBinWidth=MaxdtBinWidth,
                                               pplmin_frac=0.9, pplmax_frac=1.1,forceCutoff=IntensityCutoff,
                                               edgeCutoff=EdgeCutoff, peakMaskSize = 5, figureNumber=2, iccFitDict=iccFitDict)

#Do some annotation
plt.figure(1)
ax = plt.gca()
paramNames = ['alpha', 'beta', 'R', 'T0', 'scale', 'HatWidth', 'KConv', 'bgOffset', 'bgSlope', 'chiSq']
annotation = ''.join(['%s %4.4f\n'%(a,b) for a,b in zip(paramNames, mtd['fit_Parameters'].column(1))])
anchored_text = AnchoredText(annotation[:-1],loc=2)
ax.add_artist(anchored_text)
plt.title('%s d=%4.4f wl=%4.4f'%(str(peak.getHKL()),peak.getDSpacing(), peak.getWavelength()))

pySlice.simpleSlices(n_events, Y3D)




