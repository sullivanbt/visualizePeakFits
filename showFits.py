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
from scipy.ndimage.filters import convolve


# These are parameters for fitting.
qLow = -25  # Lowest value of q for ConvertToMD 
qHigh = 25; # Highest value of q for ConvertToMD
Q3DFrame='Q_lab' # Either 'Q_lab' or 'Q_sample'; Q_lab recommended if using a strong peaks
                 # profile library from a different sample
eventFileName = '/SNS/TOPAZ/IPTS-18474/data/TOPAZ_26751_event.nxs' #Full path to the event nexus file
peaksFile = '/SNS/TOPAZ/IPTS-18474/shared/SC100K_useDetCal/26751_Niggli.integrate' #Full path to the ISAW peaks file
UBFile = '/SNS/TOPAZ/IPTS-18474/shared/SC100K_useDetCal/26751_Niggli.mat' #Full path to the ISAW UB file
#strongPeakParamsFile = '/SNS/MANDI/shared/ProfileFitting/strongPeakParams_beta_lac_mut_mbvg.pkl' #Full path to pkl file
strongPeakParamsFile = '/SNS/users/ntv/integrate/strongPeakParams_sc2018.pkl'
strongPeakParamsFile = None
#moderatorCoefficientsFile = '/SNS/MANDI/shared/ProfileFitting/franz_coefficients_2017.dat' #Full path to pkl file
moderatorCoefficientsFile = '/home/ntv/integrate/bl11_moderatorCoefficients_2018.dat' #Full path to pkl file
IntensityCutoff = 200 # Profile paramters from a nearby strong Bragg peak will be forced for peaks with counts below IntensityCutoff.
EdgeCutoff = 10 # Profile parameters from a nearby strong Bragg peak will be forced for peaks within EdgeCutoff of a detector edge.
FracStop = 0.05 # Fraction of max counts to include in peak selection.
MinpplFrac = 0.9 # Min fraction of predicted background level to check
MaxpplFrac = 1.2 # Max fraction of predicted background level to check
DQMax = 0.2 # Largest total side length (in Angstrom-1) to consider for profile fitting.
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
        #SetInstrumentParameter(Workspace='peaks_ws', ParameterName='iccB', ParameterType='String', Value='0.001 0.3 0.005')
        SetInstrumentParameter(Workspace='peaks_ws', ParameterName='iccKConv', ParameterType='String', Value='10.0 10000.0 400.0')
        #SetInstrumentParameter(Workspace='peaks_ws', ParameterName='iccR', ParameterType='String', Value='0.0 1.0 0.05')

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


try:
    print('Current peak is %i'%peakNumber)
except:
    pass
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
peakMaskSize = peaks_ws.getInstrument().getIntParameter("peakMaskSize")[0]
np.warnings.filterwarnings('ignore') # There can be a lot of warnings for bad solutions.

q_frame = 'lab'
#Get some peak variables
peak = peaks_ws.getPeak(peakNumber)
Box = ICCFT.getBoxFracHKL(peak, peaks_ws, MDdata, UBMatrix, peakNumber, dQ, fracHKL=0.5, dQPixel=DQPixel,  q_frame=q_frame)
box = Box
#Set up our filters
qMask = ICCFT.getHKLMask(UBMatrix, frac=FracHKL, dQPixel=DQPixel, dQ=dQ)
n_events = Box.getNumEventsArray()


iccFitDict = ICCFT.parseConstraints(peaks_ws)
Y3D, goodIDX, pp_lambda2, params2 = BVGFT.get3DPeak(peak, peaks_ws, box, padeCoefficients,qMask,nTheta=NTheta, nPhi=NPhi,
                                               plotResults=plotResults,
                                               zBG=1.96,fracBoxToHistogram=1.0,bgPolyOrder=1, strongPeakParams=strongPeakParams,
                                               q_frame=q_frame, mindtBinWidth=MindtBinWidth, maxdtBinWidth=MaxdtBinWidth,
                                               pplmin_frac=MinpplFrac, pplmax_frac=MaxpplFrac,forceCutoff=IntensityCutoff,
                                               edgeCutoff=EdgeCutoff, peakMaskSize = peakMaskSize, figureNumber=2, iccFitDict=iccFitDict)

#Calcualte I and Sigma I
neigh_length_m=3
peakIDX = Y3D/Y3D.max() > 0.05
intensity = np.sum(Y3D[peakIDX])
convBox = 1.0*np.ones([neigh_length_m, neigh_length_m,neigh_length_m]) / neigh_length_m**3
conv_n_events = convolve(n_events,convBox)
bgIDX = reduce(np.logical_and,[~goodIDX, qMask, conv_n_events>0])
bgEvents = np.mean(n_events[bgIDX])*np.sum(peakIDX)
w_events = n_events.copy()
w_events[w_events==0] = 1
varFit = np.average((n_events[peakIDX]-Y3D[peakIDX])*(n_events[peakIDX]-Y3D[peakIDX]), weights=(w_events[peakIDX]))
sigma = np.sqrt(intensity + bgEvents + varFit)

intensityData = np.sum(n_events[peakIDX])
sigmaData = np.sqrt(intensityData + bgEvents + varFit)

peakI = peak.getIntensity()
peakS = peak.getSigmaIntensity()
print('Elliptical          : %4.2f +- %4.2f (%4.2f)'%(peakI, peakS, 1.*peakI/peakS))
print('Profile fitted model: %4.2f +- %4.2f (%4.2f)'%(intensity, sigma, 1.*intensity/sigma))
print('Profile fitted data : %4.2f +- %4.2f (%4.2f)'%(intensityData, sigmaData, 1.*intensityData/sigmaData))
#Do some annotation
if plotResults:
    plt.figure(1)
    ax = plt.gca()
    paramNames = ['alpha', 'beta', 'R', 'T0', 'scale', 'HatWidth', 'KConv', 'bgOffset', 'bgSlope', 'chiSq']
    annotation = ''.join(['%s %4.4f\n'%(a,b) for a,b in zip(paramNames, mtd['fit_Parameters'].column(1))])
    annotation += 'E (meV) %4.4f'%peak.getInitialEnergy()
    anchored_text = AnchoredText(annotation[:-1],loc=2)
    ax.add_artist(anchored_text)
    plt.title('%s d=%4.4f wl=%4.4f'%(str(peak.getHKL()),peak.getDSpacing(), peak.getWavelength()))

#Show interactive slices
pySlice.simpleSlices(n_events, Y3D)




