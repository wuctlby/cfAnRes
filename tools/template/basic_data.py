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

import matplotlib.pyplot as plt
import math
import io


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

def combine_figures(figs, title="", cols=3):
    """将多个 Figure 对象组合到一个 16:9 大图中"""
    # 计算行数
    rows = int(np.ceil(len(figs) / cols))

    # 创建 16:9 大画布
    fig_combined, axes = plt.subplots(rows, cols, figsize=(16, 9))
    fig_combined.suptitle(title)

    # axes 可能是二维或一维数组
    axes = np.atleast_2d(axes)

    for idx, fig in enumerate(figs):
        r, c = divmod(idx, cols)
        ax = axes[r][c]

        # 从原 fig 拿到 Axes 对象，并绘制到目标 ax 上
        for original_ax in fig.axes:
            for artist in original_ax.get_children():
                try:
                    artist.remove()  # 防止重复引用
                except:
                    pass
            for line in original_ax.get_lines():
                ax.plot(line.get_xdata(), line.get_ydata(), label=line.get_label())

        # 复制标题和标签
        ax.set_title(original_ax.get_title())
        ax.set_xlabel(original_ax.get_xlabel())
        ax.set_ylabel(original_ax.get_ylabel())
        ax.legend()

    # 删除多余空白子图
    for j in range(len(figs), rows * cols):
        fig_combined.delaxes(axes.flatten()[j])

    plt.tight_layout()
    return fig_combined

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
    mark = 'FT0C'

    # process THnSparse
    listHistos = []
    ratios = []
    for iFile, inputFile in enumerate(inputFiles):
        listHistos.append([])
        for thn in taskThns:
            histos = multi_proj(inputFile, thn, axisProj, axisFilters)
            listHistos[-1].extend(histos)
            # histos = multi_proj(inputFile, thn, 'score_FD', axisFilters)
            # listHistos[-1].extend(histos)

    histoNames = []
    projFile = TFile(f'./proj_{mark}.root', 'RECREATE')
    for listHisto in listHistos:
        for iHisto, histo in enumerate(listHisto):
            histoNames.append(histo.GetName())
            histo.Write()
    projFile.Close()

    fitters = []
    plotMassFits, plotRawResiduals = [], []
    means, sigmas, rawYields, roverss, significances = [], [], [], [], []
    meanErrors, sigmaErrors, rawYieldErrors, roversErrors, significancesErrors = [], [], [], [], []
    for histName in histoNames:

        massDis = DataHandler(
            data=f'proj_{mark}.root',
            histoname=histName,
            limits=[1.72, 2.02],
            rebin=2
        )

        
        fitters.append(
            F2MassFitter(
                data_handler=massDis,
                name_signal_pdf=['gaussian'],
                name_background_pdf=['expo'],
                name=f"{histName}_fit"
            )
        )

        fitters[-1].set_particle_mass(0, pdg_id=421, limits=[1.850, 1.890])

        fitters[-1].mass_zfit()

        plotMassFit, _ = fitters[-1].plot_mass_fit(style="ATLAS", show_extra_info=False, extra_info_loc=['upper left', 'lower right'], axis_title=r"$M_{K\pi} (\mathrm{GeV}/c^2)$")
        plotRawResidual = fitters[-1].plot_raw_residuals(style="ATLAS", axis_title=r"$M_{pK\pi}$")
        mean, meanError = fitters[-1].get_mass()
        sigma, sigmaError = fitters[-1].get_sigma()
        rawyield, rawyieldError = fitters[-1].get_raw_yield(0)
        rovers, roversError = fitters[-1].get_signal_over_background()
        significance, signifError = fitters[-1].get_significance()

        plotMassFits.append(plotMassFit)
        plotRawResiduals.append(plotRawResidual)
        means.append(mean)
        meanErrors.append(meanError)
        sigmas.append(sigma)
        sigmaErrors.append(sigmaError)
        rawYields.append(rawyield)
        rawYieldErrors.append(rawyieldError)
        roverss.append(rovers)
        significances.append(significance)
        roversErrors.append(roversError)
        significancesErrors.append(signifError)

    # prepare the output
    occuBinning = sorted({edge for filters in axisFilters['occ'] for edge in filters})
    occuLows = occuBinning[:-1]
    occuHighs = occuBinning[1:]
    hMean = TH1F('hMean', 'Mean vs occupancy;Occupancy;Mean (GeV/c^2)', len(occuBinning) - 1, np.array(occuBinning, dtype="float64"))
    hSigma = TH1F('hSigma', 'Sigma vs occupancy;Occupancy;Sigma (GeV/c^2)', len(occuBinning) - 1, np.array(occuBinning, dtype="float64"))
    hRawYield = TH1F('hRawYield', 'Raw Yield vs occupancy;Occupancy;Raw Yield', len(occuBinning) - 1, np.array(occuBinning, dtype="float64"))

    for iBin, (mean, meanError, sigma, sigmaError, rawyield, rawyieldError, rovers, significances, significancesError, occuLows, occuHighs) in enumerate(zip(
        means, meanErrors, sigmas, sigmaErrors, rawYields, rawYieldErrors, roverss, significances, significancesErrors, occuLows, occuHighs
    )):
        hMean.SetBinContent(iBin + 1, mean)
        hMean.SetBinError(iBin + 1, meanErrors[iBin])
        hSigma.SetBinContent(iBin + 1, sigma)
        hSigma.SetBinError(iBin + 1, sigmaErrors[iBin])
        hRawYield.SetBinContent(iBin + 1, rawyield)
        hRawYield.SetBinError(iBin + 1, rawYieldErrors[iBin])
        
        latex_expr_mass = f"$\mu=$" + str(round(mean, 4)) + f"$\pm$" + str(round(meanError, 4)) + "$\;\mathrm{GeV}/c^2$"
        latex_expr_width = f"$\sigma=$" + str(round(sigma, 3)) + f"$\pm$" + str(round(sigmaError, 3)) + "$\;\mathrm{GeV}/c^2$"
        latex_expr_rawyield = r"$N_{\mathrm{D}^{\mathrm{0}}}\:=\:$" + str(int(rawyield)) + f"$\:\pm\:$" + str(int(rawyieldError))
        latex_expr_rovers = r"$r/s\:=\:$" + str(round(rovers, 3))
        latex_expr_significance = r"$s/\sqrt{s+b}\:(3\sigma)\:=\:$" + str(round(significance, 1)) + f"$\:\pm\:$" + str(round(significancesError, 1))
        latex_expr_occupancy = f"{occuLows} <" + r"$\mathrm{Occupancy}$" + f'({mark})' + f" < {occuHighs}"
        
        plotMassFits[iBin].text(0.195, 0.83, r'$\mathrm{D^0} \rightarrow \mathrm{K}^- \mathrm{\pi}^+ + \mathrm{c.c.}$')
        plotMassFits[iBin].text(0.195, 0.79, rf'{30-50}$\% ~\mathrm{{Pb-Pb}}, \sqrt{{\it{{s}}_\mathrm{{NN}}}} = 5.36 \: \mathrm{{TeV}}$')
        plotMassFits[iBin].text(0.195, 0.75, f'0 < $p_{{\mathrm{{T}}}} < 24 \: \mathrm{{GeV/c}}$')
        plotMassFits[iBin].text(0.195, 0.71, latex_expr_occupancy)
        plotMassFits[iBin].text(0.62, 0.67, latex_expr_mass, fontsize=13)
        plotMassFits[iBin].text(0.62, 0.62, latex_expr_width, fontsize=13)
        plotMassFits[iBin].text(0.195, 0.40, latex_expr_rawyield, fontsize=15)
        plotMassFits[iBin].text(0.195, 0.35, latex_expr_rovers, fontsize=15)
        
        plotRawResiduals[iBin].text(0.195, 0.83, r'$\mathrm{D^0} \rightarrow \mathrm{K}^- \mathrm{\pi}^+ + \mathrm{c.c.}$')
        plotRawResiduals[iBin].text(0.195, 0.75, f'0 < $p_{{\mathrm{{T}}}} < 24 \: \mathrm{{GeV/c}}$')
        plotRawResiduals[iBin].text(0.195, 0.71, latex_expr_occupancy)

    ry = TFile(f'./raw_yields_{mark}.root', 'RECREATE')
    hMean.Write()
    hSigma.Write()
    hRawYield.Write()
    ry.Close()

    # with PdfPages(f'./fitSpectrum_D0_PbPb_{mark}.pdf') as pdf:
    #     combined_fig = combine_figures(plotMassFits, title="Mass Fits")
    #     if combined_fig:
    #         pdf.savefig(combined_fig)

    # with PdfPages(f'./fitResiduals_D0_PbPb_{mark}.pdf') as pdf:
    #     combined_fig = combine_figures(plotRawResiduals, title="Raw Residuals")
    #     if combined_fig:
    #         pdf.savefig(combined_fig)
    
    with PdfPages(f'./fitSpectrum_D0_PbPb_{mark}.pdf') as pdf:
        for fig in plotMassFits:
            pdf.savefig(fig)

    with PdfPages(f'./fitResiduals_D0_PbPb_{mark}.pdf') as pdf:
        for fig in plotRawResiduals:
            pdf.savefig(fig)

    exit(0)
    
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
                legRatio.AddEntry(ratios[iHisto - 1 + iFile * (len(listHisto) - 1)], f"candidates / {rawLabels[iHisto + iFile * len(listHisto)]}")

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
    output.Close()