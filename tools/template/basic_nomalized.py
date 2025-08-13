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

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description='Basic script to process THnSparse')
    argparser.add_argument('config', type=str, metavar="text",
						  default="config_basic.yml", help='Configuration file for processing')
    args = argparser.parse_args()

    inputFiles = [
        '/data/meta/AnRes/flow/pass4/T03_occupancy_ITS_FT0c/occu_ratio_rebin_pickedMC_pass4_medium.root',
        '/data/meta/AnRes/flow/pass4/T03_occupancy_ITS_FT0c/occu_ratiodata_pass4_medium.root',
    ]

    mc = TFile.Open(inputFiles[0], 'read')
    mc_rec_occu_its = mc.Get('proj_11_axisscore_bkg_0_0.004_axisscore_FD_0_0.96_rebined_0')
    mc_rec_occu_its.SetDirectory(0)
    mc_evt_occu_its = mc.Get('0_0_rebined_0')
    mc_evt_occu_its.SetDirectory(0)
    for key in mc.GetListOfKeys():
        print(key.GetName())
    # ft0 = TFile.Open(inputFiles[0], 'read')
    mc_rec_occu_ft0 = mc.Get('proj_11_axisscore_bkg_0_0.004_axisscore_FD_0_0.96_rebined_1')
    mc_rec_occu_ft0.SetDirectory(0)
    mc_evt_occu_ft0 = mc.Get('1_0_rebined_1')
    mc_evt_occu_ft0.SetDirectory(0)

    data = TFile.Open(inputFiles[1], 'read')
    data_occu_its = data.Get('hRawYield;1')
    data_occu_its.SetDirectory(0)
    data_occu_ft0 = data.Get('hRawYield;2')
    data_occu_ft0.SetDirectory(0)
    
    data_evt_occu_its = data.Get('0_0_rebined')
    data_evt_occu_its.SetDirectory(0)
    data_evt_occu_ft0 = data.Get('1_0_rebined')
    data_evt_occu_ft0.SetDirectory(0)

    listHisto = [
        mc_rec_occu_its, mc_evt_occu_its,
        data_occu_its, data_evt_occu_its,
        mc_rec_occu_ft0, mc_evt_occu_ft0,
        data_occu_ft0, data_evt_occu_ft0
    ]

    label = [
        'MC rec', 'MC evt',
        'Data raw yields', 'Data evt', 
        'MC rec', 'MC evt',
        'Data raw yields', 'Data evt'
    ]

    ratio_labels = [
        'MC rec / Data raw yields', 'MC / Data evt',
        'MC rec / Data raw yields', 'MC / Data evt'
    ]

    for histo in listHisto:
        histo.Scale(1.0 / histo.Integral(), "width")
        histo.Sumw2()

    ratios = []
    for iHisto in [0,1,4,5]:
        ratios.append(compute_ratio(listHisto[iHisto], listHisto[iHisto+2]))

    # prepare the output
    gROOT.SetBatch(True)
    canvas = TCanvas('canvas', 'Occupancy distribution of D0 in 30-50% Pb-Pb collisions', 2400, 1200)
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
    legend.SetTextSize(0.035)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetLineColor(0)
    legend.SetFillColor(0)
    legend.SetFillStyle(1)
    
    legRatio = legend.Clone()

    text = TText(0.70, 0.90, "")
    text.SetTextSize(0.035)
    text.SetTextFont(42)
    text.SetTextAlign(22)
    text.SetNDC()
    text.SetTextColor(kBlack)

    markerStyles = [kFullCircle, kFullSquare, kFullDiamond, kFullCross, kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross]
    colors = [kBlack, kRed, kBlue, kGreen, kBlack, kRed, kBlue, kGreen]

    maxX, maxY = 0, 0
    minX = 0 if not logX else 0.1
    minY = 0 if not logY else 0.1

    iRatio = 0
    for iHisto, histo in enumerate(listHisto):
        if iHisto == 0:
            e = legend.AddEntry('', 'ITS', '')
            e.SetTextFont(62)
        elif iHisto == 4:
            e = legend.AddEntry('', 'FT0C', '')
            e.SetTextFont(62)

        histo.SetMarkerStyle(markerStyles[iHisto])
        histo.SetMarkerColor(colors[iHisto])
        histo.SetMarkerSize(1.2)
        histo.SetLineColor(colors[iHisto])
        histo.SetLineWidth(2)
        histo.SetStats(0)
        legend.AddEntry(histo, f"{label[iHisto]}")
        
        maxY = max(maxY, histo.GetMaximum())
        maxX = max(maxX, histo.GetXaxis().GetXmax())
        minX = min(minX, histo.GetXaxis().GetXmin()) if not logX else min(minX, histo.GetXaxis().GetXmin() if histo.GetXaxis().GetXmin() > 0 else 0.1)
        minY = min(minY, histo.GetMinimum()) if not logY else min(minY, histo.GetMinimum() if histo.GetMinimum() > 0 else 0.1)
        if iHisto in [0, 1, 4, 5]:
            ratios[iRatio].SetMarkerStyle(markerStyles[iHisto])
            ratios[iRatio].SetMarkerColor(colors[iHisto])
            ratios[iRatio].SetMarkerSize(1.2)
            ratios[iRatio].SetLineColor(colors[iHisto])
            ratios[iRatio].SetLineWidth(2)
            ratios[iRatio].SetStats(0)
            e = legRatio.AddEntry(histo, f"{ratio_labels[iRatio]}")
            iRatio += 1

    canvas.cd(1).DrawFrame(500, minY*0.8, maxX, maxY*1.4, "Occupancy distribution of D0 in 30-50% Pb-Pb collisions;Occupancy;Normalized")
    canvas.cd(2).DrawFrame(500, 0.2, maxX, 20, "Ratio of occupancy distributions;Occupancy;Ratio of mc to data")
    canvas.cd(1)
    for iHisto, histo in enumerate(listHisto):
        histo.Draw("same")

    legend.Draw()
    # text.Draw()
    canvas.Modified()
    canvas.Update()

    canvas.cd(2)
    for iFile, ratio in enumerate(ratios):
        ratio.Draw('same lp')
    legRatio.Draw()
    canvas.Modified()
    canvas.Update()

    canvas.SaveAs("occupancy_distribution_picke.root")
    

    # # prepare the output
    # gROOT.SetBatch(True)
    # canvas = TCanvas('canvas', 'Occupancy distribution of D0 in 30-50% Pb-Pb collisions', 1200, 1200)
    # # canvas.Divide(2, 1)
    # canvas.SetGridy()
    # canvas.SetLeftMargin(0.12)
    # canvas.SetBottomMargin(0.10)
    # canvas.SetRightMargin(0.10)
    # canvas.SetTopMargin(0.08)
    
    # logX = True
    # logY = True
    # gridX = True
    # gridY = True
    
    # if logX:
    #     canvas.cd(1).SetLogx()
    #     canvas.cd(2).SetLogx()
    # if logY:
    #     canvas.cd(1).SetLogy()
    #     canvas.cd(2).SetLogy()
    # if gridX:
    #     canvas.cd(1).SetGridx()
    #     canvas.cd(2).SetGridx()
    # if gridY:
    #     canvas.cd(1).SetGridy()
    #     canvas.cd(2).SetGridy()

    # legend = TLegend(0.18, 0.18, 0.5, 0.5)
    # legend.SetTextSize(0.04)
    # legend.SetTextFont(42)
    # legend.SetBorderSize(0)
    # legend.SetLineColor(0)
    # legend.SetFillColor(0)
    # legend.SetFillStyle(1)
    
    # legRatio = legend.Clone()

    # text = TText(0.70, 0.90, "")
    # text.SetTextSize(0.035)
    # text.SetTextFont(42)
    # text.SetTextAlign(22)
    # text.SetNDC()
    # text.SetTextColor(kBlack)

    # markerStyles = [kFullCircle, kFullSquare, kFullCross, kOpenCircle, kOpenSquare, kOpenCross]
    # colors = [kBlack, kRed, kBlue, kGreen, kBlack, kRed, kBlue, kGreen]



    # maxX, maxY = 0, 0
    # minX = 0 if not logX else 0.1
    # minY = 0 if not logY else 0.1
    # for iFile, listHisto in enumerate(listHistos):
    #     if iFile == 0:
    #         e = legend.AddEntry('', 'ITS 452638', '')
    #         e.SetTextFont(62)
    #         er = legRatio.AddEntry('', 'ITS 452638', '')
    #         er.SetTextFont(62)
    #     if iFile == 1:
    #         e = legend.AddEntry('', 'FT0C 468349', '')
    #         e.SetTextFont(62)
    #         er = legRatio.AddEntry('', 'FT0C 468349', '')
    #         er.SetTextFont(62)

    #     for iHisto, histo in enumerate(listHisto):
    #         histo.SetMarkerStyle(markerStyles[iHisto + iFile * len(listHisto)])
    #         histo.SetMarkerColor(colors[iHisto + iFile * len(listHisto)])
    #         histo.SetMarkerSize(1.2)
    #         histo.SetLineColor(colors[iHisto + iFile * len(listHisto)])
    #         histo.SetLineWidth(2)
    #         histo.SetStats(0)  # Disable statistics box
    #         legend.AddEntry(histo, f"{rawLabels[iHisto + iFile * len(listHisto)]}")
    #         maxY = max(maxY, histo.GetMaximum())
    #         minY = min(minY, histo.GetMinimum()) if not logY else min(minY, histo.GetMinimum() if histo.GetMinimum() > 0 else 0.1)
    #         maxX = max(maxX, histo.GetXaxis().GetXmax())
    #         minX = min(minX, histo.GetXaxis().GetXmin()) if not logX else min(minX, histo.GetXaxis().GetXmin() if histo.GetXaxis().GetXmin() > 0 else 0.1)
            
    #         if iHisto != 0:
    #             ratios[iHisto - 1 + iFile * (len(listHisto) - 1)].SetMarkerStyle(markerStyles[iHisto + iFile * len(listHisto)])
    #             ratios[iHisto - 1 + iFile * (len(listHisto) - 1)].SetMarkerColor(colors[iHisto + iFile * len(listHisto)])
    #             ratios[iHisto - 1 + iFile * (len(listHisto) - 1)].SetMarkerSize(1.2)
    #             ratios[iHisto - 1 + iFile * (len(listHisto) - 1)].SetLineColor(colors[iHisto + iFile * len(listHisto)])
    #             ratios[iHisto - 1 + iFile * (len(listHisto) - 1)].SetLineWidth(2)
    #             ratios[iHisto - 1 + iFile * (len(listHisto) - 1)].SetStats(0)  # Disable statistics box
    #             legRatio.AddEntry(ratios[iHisto - 1 + iFile * (len(listHisto) - 1)], f"candidates / {rawLabels[iHisto + iFile * len(listHisto)]}")

    # canvas.cd(1).DrawFrame(500, minY, maxX, maxY*1.4, "Occupancy distribution of D0 in 30-50% Pb-Pb collisions")
    # canvas.cd(2).DrawFrame(500, 0.0001, maxX, 1.8, "")

    # # save results
    # canvas.cd(1)
    # for iFile, listHisto in enumerate(listHistos):
    #     for iHisto, histo in enumerate(listHisto):
    #         histo.Draw('same lp')

    # legend.Draw()
    # text.Draw()
    # canvas.Modified()
    # canvas.Update()

    # canvas.cd(2)
    # for iFile, ratio in enumerate(ratios):
    #     ratio.Draw('same lp')
    # legRatio.Draw()
    # canvas.Modified()
    # canvas.Update()
    
    # # input("Press Enter to save the canvas...")  # Pause for user input before saving

    # # canvas.RedrawAxis()
    # # canvas.Modified()
    # # canvas.Update()

    # # canvas.SaveAs(outputFile + '.png')
    # img = TImage.Create()
    # img.FromPad(canvas)
    # img.WriteImage(outputFile + ".png")
    # canvas.Print(outputFile + ".pdf")

    # output = TFile(outputFile + '.root', 'RECREATE')
    # canvas.Write()
    # for listHisto in listHistos:
    #     for histo in listHisto:
    #         histo.Write()
    # for ratio in ratios:
    #     ratio.Write()
    # output.Close()