#!/usr/bin/env python3
"""
data_driven_fraction.py — Compute prompt/feed-down fractions using data-driven method.

Takes cutVar.root + efficiencies → produces fraction ROOT files per cut set.

Usage:
    python3 data_driven_fraction.py cutVar.root effs_dir -o output_dir -b
"""
import argparse
import os
import sys
import ROOT

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from utils import logger, list_files, load_eff_histos, load_cutVar_histos, GetPromptFDFractionCutSet

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)


def data_driven_frac(output_dir, i_file, hEffPrompt, hEffFD,
                     hPromptFrac, hFDFrac, hPromptFracCorr, hFDFracCorr,
                     hCorrYieldPrompt, hCorrYieldFD,
                     hCovPP, hCovPF, hCovFF):
    """Compute prompt/FD fractions for one cut set and save to file."""
    for iPt in range(hEffPrompt.GetNbinsX()):
        pt_bin = iPt + 1
        
        # Check that the cutVar histogram has this bin
        if pt_bin > hCorrYieldPrompt.GetNbinsX():
            logger(f"  Skipping pt bin {pt_bin} (cutVar has only {hCorrYieldPrompt.GetNbinsX()} bins)", "WARNING")
            continue
            
        eff_p = hEffPrompt.GetBinContent(pt_bin)
        eff_f = hEffFD.GetBinContent(pt_bin)
        
        corr_p = hCorrYieldPrompt.GetBinContent(pt_bin)
        corr_f = hCorrYieldFD.GetBinContent(pt_bin)
        
        # Skip bins with no signal
        if corr_p == 0 and corr_f == 0:
            logger(f"  Skipping pt bin {pt_bin} (no corrected yield)", "WARNING")
            continue
            
        cov_pp = hCovPP.GetBinContent(pt_bin)
        cov_pf = hCovPF.GetBinContent(pt_bin)
        cov_ff = hCovFF.GetBinContent(pt_bin)

        frac, unc_frac = GetPromptFDFractionCutSet(eff_p, eff_f, corr_p, corr_f, cov_pp, cov_ff, cov_pf)
        frac_corr, unc_frac_corr = GetPromptFDFractionCutSet(1., 1., corr_p, corr_f, cov_pp, cov_ff, cov_pf)

        hPromptFrac.SetBinContent(pt_bin, frac[0])
        hPromptFrac.SetBinError(pt_bin, unc_frac[0])
        hFDFrac.SetBinContent(pt_bin, frac[1])
        hFDFrac.SetBinError(pt_bin, unc_frac[1])
        hPromptFracCorr.SetBinContent(pt_bin, frac_corr[0])
        hPromptFracCorr.SetBinError(pt_bin, unc_frac_corr[0])
        hFDFracCorr.SetBinContent(pt_bin, frac_corr[1])
        hFDFracCorr.SetBinError(pt_bin, unc_frac_corr[1])

    out_path = os.path.join(output_dir, f"frac_{i_file:02d}.root")
    outfile = ROOT.TFile(out_path, "RECREATE")
    hEffPrompt.Write()
    hEffFD.Write()
    hPromptFrac.Write()
    hFDFrac.Write()
    hPromptFracCorr.Write()
    hFDFracCorr.Write()
    outfile.Close()

    logger(f"  Saved {out_path}")


def main_data_driven_frac(cutVar_file, eff_path, output_dir="", batch=True):
    """Main entry point for data-driven fraction computation."""
    if batch:
        ROOT.gROOT.SetBatch(True)

    hEffPrompts, hEffFDs, hPromptFracs, hFDFracs, hPromptFracCorrs, hFDFracCorrs = \
        load_eff_histos(list_files(eff_path, prefix="eff_"))
    hCorrYieldPrompt, hCorrYieldFD, hCovPP, hCovPF, hCovFF = load_cutVar_histos(cutVar_file)

    out_dir = output_dir if output_dir else os.path.join(os.path.dirname(eff_path), "frac")
    os.makedirs(out_dir, exist_ok=True)
    logger(f"Output directory: {out_dir}")

    for i_file in range(len(hEffPrompts)):
        data_driven_frac(out_dir, i_file,
                         hEffPrompts[i_file], hEffFDs[i_file],
                         hPromptFracs[i_file], hFDFracs[i_file],
                         hPromptFracCorrs[i_file], hFDFracCorrs[i_file],
                         hCorrYieldPrompt, hCorrYieldFD,
                         hCovPP, hCovPF, hCovFF)

    logger(f"All fractions saved to {out_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Data-driven fraction computation")
    parser.add_argument("cutVar_file", help="Path to cutVar.root")
    parser.add_argument("eff_path", help="Path to efficiencies directory")
    parser.add_argument("-o", "--output-dir", default="", help="Output directory")
    parser.add_argument("-b", "--batch", action="store_true", default=True, help="Batch mode")
    args = parser.parse_args()

    main_data_driven_frac(args.cutVar_file, args.eff_path, args.output_dir, args.batch)
