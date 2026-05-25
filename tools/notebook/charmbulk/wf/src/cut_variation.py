#!/usr/bin/env python3
"""
cut_variation.py — Perform cut variation analysis.

Takes raw yields and efficiencies → produces corrected yields with the chi2 minimisation.
Output: {outdir}/cutVar/cutVar.root

Usage:
    python3 cut_variation.py flow_config.yml raw_yields_dir effs_dir -b
"""
import argparse
import os
import sys
import numpy as np
import yaml
import ROOT
from itertools import product

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from utils import logger, list_files, GetMinimisation

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetPaintTextFormat("4.2f")


def compute_cut_var(config_file, input_path_ry, input_path_eff, batch=True):
    """Compute cut variation: minimisation over cut sets to get corrected yields."""
    if batch:
        ROOT.gROOT.SetBatch(True)

    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    ptmins = config["ptbins"][:-1]
    ptmaxs = config["ptbins"][1:]
    nPtBins = len(ptmins)

    # Load raw yield files
    raw_yield_files = list_files(input_path_ry, prefix="raw_yields_")
    eff_files = list_files(input_path_eff, prefix="eff_")

    if len(eff_files) != len(raw_yield_files):
        logger(f"Number of eff files ({len(eff_files)}) != raw yield files ({len(raw_yield_files)})", "ERROR")
        raise ValueError("Mismatch between eff and raw yield files")

    eff_files.sort()
    raw_yield_files.sort()

    # Load histograms
    hRawYields = []
    hEffPrompt, hEffFD = [], []

    for ry_file, eff_file in zip(raw_yield_files, eff_files):
        fin_ry = ROOT.TFile.Open(ry_file)
        if not fin_ry.GetListOfKeys().Contains("hRawYieldsSimFit"):
            logger(f"File {ry_file} has no hRawYieldsSimFit, skipping", "WARNING")
            continue
        h = fin_ry.Get("hRawYieldsSimFit")
        h.SetDirectory(0)
        hRawYields.append(h)
        fin_ry.Close()

        fin_eff = ROOT.TFile.Open(eff_file)
        hEp = fin_eff.Get("hEffPrompt")
        hEf = fin_eff.Get("hEffFD")
        hEp.SetDirectory(0)
        hEf.SetDirectory(0)
        hEffPrompt.append(hEp)
        hEffFD.append(hEf)
        fin_eff.Close()

    logger(f"Loaded {len(hRawYields)} raw yield files, {len(hEffPrompt)} eff files")

    # Cutset files — try multiple locations
    raw_yields_dir = input_path_ry.rstrip("/")
    base_dir = os.path.dirname(raw_yields_dir)
    possible_dirs = [
        raw_yields_dir,
        os.path.join(base_dir, "cutsets"),
        os.path.join(base_dir, "..", "cutsets"),
    ]
    # Also check cutvar_* subdirectories in the parent
    if os.path.exists(base_dir):
        for d in os.listdir(base_dir):
            full = os.path.join(base_dir, d, "cutsets")
            if os.path.isdir(full):
                possible_dirs.append(full)
    cutset_dir = None
    for d in possible_dirs:
        if os.path.exists(d):
            for f in os.listdir(d):
                if f.endswith(".yml"):
                    cutset_dir = d
                    break
        if cutset_dir:
            break

    if cutset_dir is None:
        logger(f"Could not find cutset dir. Tried: {possible_dirs}", "ERROR")
        raise FileNotFoundError("No cutset directory found")

    cutset_files = sorted(f for f in os.listdir(cutset_dir) if f.endswith(".yml"))
    cutset_paths = [os.path.join(cutset_dir, f) for f in cutset_files]
    logger(f"Found {len(cutset_paths)} cutset files in {cutset_dir}", "INFO")

    # Create output histograms
    pt_arr = np.array(config["ptbins"], dtype="d")
    hCorrYieldPrompt = ROOT.TH1F("hCorrYieldPrompt", ";p_{T} (GeV/c);N_{prompt}", nPtBins, pt_arr)
    hCorrYieldFD = ROOT.TH1F("hCorrYieldFD", ";p_{T} (GeV/c);N_{non-prompt}", nPtBins, pt_arr)
    hCovPP = ROOT.TH1F("hCovPromptPrompt", ";p_{T} (GeV/c);#sigma(N_{p},N_{p})", nPtBins, pt_arr)
    hCovPF = ROOT.TH1F("hCovPromptFD", ";p_{T} (GeV/c);#sigma(N_{p},N_{FD})", nPtBins, pt_arr)
    hCovFP = ROOT.TH1F("hCovFDPrompt", ";p_{T} (GeV/c);#sigma(N_{FD},N_{p})", nPtBins, pt_arr)
    hCovFF = ROOT.TH1F("hCovFDFD", ";p_{T} (GeV/c);#sigma(N_{FD},N_{FD})", nPtBins, pt_arr)

    # Per-pT bin minimisation
    for iPt in range(nPtBins):
        pt_min, pt_max = ptmins[iPt], ptmaxs[iPt]
        logger(f"Processing pT bin {pt_min:.1f}–{pt_max:.1f} GeV/c (bin {iPt + 1}/{nPtBins})")

        list_ry, list_ry_unc = [], []
        list_ep, list_ep_unc = [], []
        list_ef, list_ef_unc = [], []

        last_fd_min = None
        for iCut, (hRy, hEp, hEf) in enumerate(zip(hRawYields, hEffPrompt, hEffFD)):
            # Skip duplicate FD cuts (correlated case — same fd_min repeats)
            cutset_path = cutset_paths[iCut]
            with open(cutset_path) as cf:
                cs = yaml.safe_load(cf)
            fd_min_this = cs["ScoreFD"]["min"][iPt]
            if last_fd_min is not None and fd_min_this == last_fd_min:
                continue
            last_fd_min = fd_min_this

            list_ry.append(hRy.GetBinContent(iPt + 1))
            list_ry_unc.append(hRy.GetBinError(iPt + 1))
            list_ep.append(hEp.GetBinContent(iPt + 1))
            list_ep_unc.append(hEp.GetBinError(iPt + 1))
            list_ef.append(hEf.GetBinContent(iPt + 1))
            list_ef_unc.append(hEf.GetBinError(iPt + 1))

        if len(list_ry) < 2:
            logger(f"  Not enough cut sets ({len(list_ry)}) for pT bin {iPt}, skipping", "ERROR")
            continue

        # Apply minimisation
        corrYields, covMatrix, chi2, matrices = GetMinimisation(
            list_ep, list_ef, list_ry,
            list_ep_unc, list_ef_unc, list_ry_unc
        )

        hCorrYieldPrompt.SetBinContent(iPt + 1, float(corrYields[0]))
        hCorrYieldPrompt.SetBinError(iPt + 1, np.sqrt(abs(float(covMatrix[0, 0]))))
        hCorrYieldFD.SetBinContent(iPt + 1, float(corrYields[1]))
        hCorrYieldFD.SetBinError(iPt + 1, np.sqrt(abs(float(covMatrix[1, 1]))))
        hCovPP.SetBinContent(iPt + 1, float(covMatrix[0, 0]))
        hCovPF.SetBinContent(iPt + 1, float(covMatrix[0, 1]))
        hCovFP.SetBinContent(iPt + 1, float(covMatrix[1, 0]))
        hCovFF.SetBinContent(iPt + 1, float(covMatrix[1, 1]))

        logger(f"  N_prompt={float(corrYields[0]):.1f}±{np.sqrt(abs(float(covMatrix[0, 0]))):.1f}, "
               f"N_FD={float(corrYields[1]):.1f}±{np.sqrt(abs(float(covMatrix[1, 1]))):.1f}, χ²={chi2:.3f}")

    # Save
    out_dir = os.path.join(os.path.dirname(input_path_ry), "cutVar")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "cutVar.root")

    outfile = ROOT.TFile(out_path, "RECREATE")
    hCorrYieldPrompt.Write()
    hCorrYieldFD.Write()
    hCovPP.Write()
    hCovPF.Write()
    hCovFP.Write()
    hCovFF.Write()
    outfile.Close()

    logger(f"Cut variation results saved to {out_path}")
    return out_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cut variation analysis")
    parser.add_argument("config", help="Flow configuration YAML")
    parser.add_argument("ry_path", help="Path to raw yields directory")
    parser.add_argument("eff_path", help="Path to efficiencies directory")
    parser.add_argument("-b", "--batch", action="store_true", default=True, help="Batch mode")
    args = parser.parse_args()

    compute_cut_var(args.config, args.ry_path, args.eff_path, args.batch)
