import os
import sys
import argparse
import yaml
from itertools import product
from ROOT import TCanvas, TLegend, TText, TFile, gROOT, TImage, TH1F, TMath
from ROOT import kBlack, kWhite, kGray, kRed, kBlue, kGreen # pylint: disable=import-error,no-name-in-module
from ROOT import kFullCircle, kFullSquare, kFullDiamond, kFullCross, kFullTriangleUp, kFullTriangleDown # pylint: disable=import-error,no-name-in-module
from ROOT import kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross, kOpenTriangleUp, kOpenTriangleDown # pylint: disable=import-error,no-name-in-module
from pathlib import Path
sys.path.append('/home/wuct/ALICE/local/RTools/cfAnRes/')
from utils.get_method import get_paths
from cfanres.doc_loader import hytask_thn_loader
from cfanres.thn_operator import multi_proj, get_task_histo

from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

def get_labels(rawLabels, thns=None, *sublabels):
    if thns is not None:
        thnsKeys = [thn[1].split('_')[0] for thn in thns]
    else:
        thnsKeys = []
    if len(sublabels) > 0:
        subComb = [list(sublabel) for sublabel in sublabels] # convert to list for combination
    else:
        subComb = []
    comb = [rawLabels] + [thnsKeys] + subComb # convert to list for combination
    labels = []
    for label in product(*comb):
        if isinstance(label, tuple):
            label = ' '.join(label)
        labels.append(label)
    return labels

def ComputeEfficiency(recoCounts, genCounts, recoCountsError, genCountsError):
    '''
    Method to compute efficiency

    Parameters
    ----------
    - recoCounts: number of reconstructed D
    - genCounts: number of genertated D
    - recoCountsError: error on number of reconstructed D
    - genCountsError: error on number of generated D

    Returns
    ----------
    - efficiency, error on efficiency
    '''
    hTmpNum = TH1F('hTmpNum', '', 1, 0, 1)
    hTmpDen = TH1F('hTmpDen', '', 1, 0, 1)
    hTmpNum.SetBinContent(1, recoCounts)
    hTmpDen.SetBinContent(1, genCounts)
    hTmpNum.SetBinError(1, recoCountsError)
    hTmpDen.SetBinError(1, genCountsError)
    hTmpNum.Divide(hTmpNum, hTmpDen, 1., 1, 'B')

    return hTmpNum.GetBinContent(1), hTmpNum.GetBinError(1)

def divide_full_corr(num, den, num_err, den_err, mode="rho1"):
    """
    mode:
      - 'rho1' : full_corr rho=1 → cov = num_err * den_err
      - 'thin' : den - Poisson   → cov = num_err**2
      - float  : rho ∈ [0,1]     → cov = rho * num_err * den_err
    """
    if den <= 0:
        return 0.0, 0.0

    r = num / den
    dRdNum = 1.0 / den
    dRdDen = -num / (den**2)

    if mode == "rho1":
        cov = num_err * den_err # rho = 1.0
    elif mode == "thin":
        cov = num_err**2  # rho = num_err / den_err
    elif isinstance(mode, (float, int)):
        rho = max(0.0, min(1.0, float(mode)))
        cov = rho * num_err * den_err
    else:
        cov = 0.0

    sigma2 = (dRdNum**2) * (num_err**2) + (dRdDen**2) * (den_err**2) + 2*dRdNum*dRdDen * cov

    if sigma2 < 0 and abs(sigma2) < 1e-18:
        sigma2 = 0.0
    sigma = TMath.Sqrt(sigma2) if sigma2 > 0 else 0.0
    return r, sigma

def compute_ratio(num_hist, den_hist, opt='eff'):
    """
    calculate the ratio of two histograms, assuming the errors are correlated
    """
    ratio_hist = den_hist.Clone(f'ratio_{num_hist.GetName()}_{den_hist.GetName()}')
    ratio_hist.Reset()

    for iBin in range(1, den_hist.GetNbinsX() + 1):
        num = num_hist.GetBinContent(iBin)
        num_err = num_hist.GetBinError(iBin)
        den = den_hist.GetBinContent(iBin)
        den_err = den_hist.GetBinError(iBin)

        if den != 0:
            if opt == 'eff':
                r, re = ComputeEfficiency(num, den, num_err, den_err)
            elif opt == 'full_corr':
                r, re = divide_full_corr(num, den, num_err, den_err, mode="rho1")
        else:
            r, re = 0, 0

        ratio_hist.SetBinContent(iBin, r)
        ratio_hist.SetBinError(iBin, re)

    return ratio_hist

# def proj_rebin(histo2D, proj_axis, binning):
#     hProjRebin = TH1F('hProjRebin', '', len(binning) - 1, np.array(binning, dtype="float64"))
#     binMins = binning[:-1]
#     binMaxs = binning[1:]
    
#     histo1D = histo2D.ProjectionY('histo1D')
#     for iBin, (binMin, binMax) in enumerate(zip(binMins, binMaxs)):
#         oBinMin = histo1D.GetXaxis().FindBin(binMin*1.00001)
#         oBinMax = histo1D.GetXaxis().FindBin(binMax*0.00009)
#         histo2Dtemp = histo2D.Clone(f'temp_2D_{iBin}')
#         histo1Dtemp = histo2Dtemp.ProjectionY(f'temp_proj_{iBin}', oBinMin, oBinMax, 'e')
        

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description='Basic script to process THnSparse')
    argparser.add_argument('config', type=str, metavar="text",
						  default="config_basic.yml", help='Configuration file for processing')
    args = argparser.parse_args()

    # Load configuration
    with open(args.config, 'r') as cfg:
        config = yaml.load(cfg, Loader=yaml.FullLoader)

    inputs = config['inputs']
    taskThns = config['taskThns']
    taskHistos = config['taskHistos']
    rawLabels = config['labels']
    axisProj = config['axisProj']
    axisFilters = config['axisFilters']
    outfile = config['outfile']
    dataset = config['dataset']
    outputFile = outfile + dataset
    inputFiles = [get_paths(input, fileType='.root') for input in inputs]
    marks = ['ITS', 'FT0C']
    # occuBinning = sorted({edge for filters in axisFilters['occ'] for edge in filters})
    # print(occuBinning)
    # exit()
    occuBinnings = [
        [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 20000],
        [0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 200000]
    ]

    # process THnSparse
    listHistos = []
    ratios = []
    for iFile, inputFile in enumerate(inputFiles):
        listHistos.append([])
        # for thn in taskThns:
        #     histos = multi_proj(inputFile, thn, axisProj, axisFilters)
        #     listHistos[-1].extend(histos)
            # histos = multi_proj(inputFile, thn, 'score_FD', axisFilters)
            # listHistos[-1].extend(histos)
        ry = TFile(f'./raw_yields_{marks[iFile]}.root', 'READ')
        hRawYields = ry.Get('hRawYield')
        hRawYields.SetDirectory(0)  # Detach from file
        listHistos[-1].append(hRawYields)
    
        histos = [get_task_histo(inputFile, taskHisto).ProjectionY(f'{iFile}_{iHisto}') for iHisto, taskHisto in enumerate(taskHistos)]
        
        for iHisto, taskHisto in enumerate(histos):
            histoRebined = taskHisto.Rebin(len(occuBinnings[iFile]) - 1, f'{taskHisto.GetName()}_rebined', np.array(occuBinnings[iFile], dtype="float64"))

        listHistos[-1].append(histoRebined)
        
        ratios.append(compute_ratio(listHistos[-1][0], listHistos[-1][1], opt='full_corr'))
    
    gROOT.SetBatch(True)
    canvas = TCanvas('canvas', '', 2400, 1200)
    canvas.Divide(2, 1)
    canvas.SetGridy()
    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.10)
    canvas.SetRightMargin(0.10)
    canvas.SetTopMargin(0.08)
    
    logX = True
    logY = True
    gridX = True
    gridY = True
    
    if logX:
        canvas.cd(1).SetLogx()
        canvas.cd(2).SetLogx()
    if logY:
        canvas.cd(1).SetLogy()
        canvas.cd(2).SetLogy()
    if gridX:
        canvas.cd(1).SetGridx()
        canvas.cd(2).SetGridx()
    if gridY:
        canvas.cd(1).SetGridy()
        canvas.cd(2).SetGridy()

    legend = TLegend(0.18, 0.18, 0.5, 0.5)
    legend.SetTextSize(0.04)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetLineColor(0)
    legend.SetFillColor(0)
    legend.SetFillStyle(1)
    
    legRatio = legend.Clone()

    text = TText(0.70, 0.90, dataset)
    text.SetTextSize(0.035)
    text.SetTextFont(42)
    text.SetTextAlign(22)
    text.SetNDC()
    text.SetTextColor(kBlack)

    markerStyles = [kFullCircle, kFullCross, kOpenCircle, kOpenCross]
    colors = [kBlack, kBlack, kBlue, kBlue]

    maxX, maxY = 0, 0
    minX = 0 if not logX else 0.1
    minY = 0 if not logY else 0.1
    for iFile, listHisto in enumerate(listHistos):
        if iFile == 0:
            e = legend.AddEntry('', 'ITS 452638', '')
            e.SetTextFont(62)
            er = legRatio.AddEntry('', 'ITS 452638', '')
            er.SetTextFont(62)
        if iFile == 1:
            e = legend.AddEntry('', 'FT0C 468349', '')
            e.SetTextFont(62)
            er = legRatio.AddEntry('', 'FT0C 468349', '')
            er.SetTextFont(62)

        for iHisto, histo in enumerate(listHisto):
            histo.SetMarkerStyle(markerStyles[iHisto + iFile * len(listHisto)])
            histo.SetMarkerColor(colors[iHisto + iFile * len(listHisto)])
            histo.SetMarkerSize(1.2)
            histo.SetLineColor(colors[iHisto + iFile * len(listHisto)])
            histo.SetLineWidth(2)
            histo.SetStats(0)  # Disable statistics box
            legend.AddEntry(histo, f"{rawLabels[iHisto + iFile * len(listHisto)]}")
            maxY = max(maxY, histo.GetMaximum())
            minY = min(minY, histo.GetMinimum()) if not logY else min(minY, histo.GetMinimum() if histo.GetMinimum() > 0 else 0.1)
            maxX = max(maxX, histo.GetXaxis().GetXmax())
            minX = min(minX, histo.GetXaxis().GetXmin()) if not logX else min(minX, histo.GetXaxis().GetXmin() if histo.GetXaxis().GetXmin() > 0 else 0.1)
            
            if iHisto != 0:
                ratios[iHisto - 1 + iFile * (len(listHisto) - 1)].SetMarkerStyle(markerStyles[iHisto + iFile * len(listHisto)])
                ratios[iHisto - 1 + iFile * (len(listHisto) - 1)].SetMarkerColor(colors[iHisto + iFile * len(listHisto)])
                ratios[iHisto - 1 + iFile * (len(listHisto) - 1)].SetMarkerSize(1.2)
                ratios[iHisto - 1 + iFile * (len(listHisto) - 1)].SetLineColor(colors[iHisto + iFile * len(listHisto)])
                ratios[iHisto - 1 + iFile * (len(listHisto) - 1)].SetLineWidth(2)
                ratios[iHisto - 1 + iFile * (len(listHisto) - 1)].SetStats(0)  # Disable statistics box
                legRatio.AddEntry(ratios[iHisto - 1 + iFile * (len(listHisto) - 1)], f"raw yields / {rawLabels[iHisto + iFile * len(listHisto)]}")

    canvas.cd(1).DrawFrame(500, minY, maxX, maxY*1.4, "Occupancy distribution of D0 in 30-50% Pb-Pb collisions")
    canvas.cd(2).DrawFrame(500, 0.0001, maxX, 1.2, "")

    # save results
    canvas.cd(1)
    for iFile, listHisto in enumerate(listHistos):
        for iHisto, histo in enumerate(listHisto):
            histo.Draw('same lp')

    legend.Draw()
    text.Draw()
    canvas.Modified()
    canvas.Update()

    canvas.cd(2)
    for iFile, ratio in enumerate(ratios):
        ratio.Draw('same lp')
    legRatio.Draw()
    canvas.Modified()
    canvas.Update()
    
    # input("Press Enter to save the canvas...")  # Pause for user input before saving

    # canvas.RedrawAxis()
    # canvas.Modified()
    # canvas.Update()

    # canvas.SaveAs(outputFile + '.png')
    img = TImage.Create()
    img.FromPad(canvas)
    img.WriteImage(outputFile + ".png")
    canvas.Print(outputFile + ".pdf")

    output = TFile(outputFile + '.root', 'RECREATE')
    canvas.Write()
    for iFile, listHisto in enumerate(listHistos):
        for iHisto, histo in enumerate(listHisto):
            histo.Write()
    for ratio in ratios:
        ratio.Write()
    output.Close()