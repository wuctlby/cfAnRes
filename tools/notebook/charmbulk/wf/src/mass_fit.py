#!/usr/bin/env python3
"""
mass_fit.py — Fit invariant mass distributions from projection files
using flarefly's F2MassFitter.

Takes projection ROOT file → fits hMassData per pT bin → outputs raw yields histograms.

Usage:
    CUDA_VISIBLE_DEVICES=-1 python3 mass_fit.py flow_config.yml proj_file.root -o output_dir
"""
import argparse
import os
import sys
import yaml
import numpy as np
import ROOT

# Add flarefly
flarefly_dir = os.environ.get("FLAREFLY_DIR", "/home/wuct/ALICE/reps/flarefly")
if flarefly_dir not in sys.path:
    sys.path.insert(0, flarefly_dir)

from flarefly import F2MassFitter
from flarefly.data_handler import DataHandler

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "."))
from utils import logger

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)


# Map from our naming to flarefly naming
SGN_MAP = {
    "kGaus": "gaussian",
    "k2Gaus": "doublegaus",
    "k2GausSigmaRatioPar": "doublegaus",
    "crystalball": "crystalball",
    "doublecb": "doublecb",
}

BKG_MAP = {
    "kExpo": "expo",
    "kLin": "chebpol",
    "kPol2": "chebpol",
    "kPol3": "chebpol",
}


def get_bkg_config(bkg_str):
    """Return (pdf_name, kwargs_for_bkg) from our string."""
    if bkg_str == "kExpo":
        return "expo", {"expo": None}
    elif bkg_str == "kLin":
        return "chebpol", {"chebpol": 1}
    elif bkg_str == "kPol2":
        return "chebpol", {"chebpol": 2}
    elif bkg_str == "kPol3":
        return "chebpol", {"chebpol": 3}
    else:
        logger(f"Unknown bkg {bkg_str}, falling back to expo", "WARNING")
        return "expo", {"expo": None}


def try_fit(h_mass, mass_min, mass_max, rebin, sgn_name, bkg_name,
            bkg_pdf_arg, sig0):
    """Attempt a mass fit; returns (fitter, True) or (None, False)."""
    h = h_mass.Clone()
    h.SetDirectory(0)
    h.GetXaxis().SetRangeUser(mass_min, mass_max)
    if rebin > 1:
        h.Rebin(rebin)

    dh = DataHandler(h, var_name="mass", limits=[mass_min, mass_max])
    fitter = F2MassFitter(
        dh,
        name_signal_pdf=[sgn_name],
        name_background_pdf=bkg_pdf_arg,
    )

    # Set initial parameters
    fitter.set_signal_initpar(0, "mu", 1.86)
    fitter.set_signal_initpar(0, "sigma", sig0)
    fitter.set_background_initpar(0, "lamb", -0.5)

    fitter.mass_zfit(do_prefit=False)
    return fitter, True


def run_mass_fit(config_file, proj_file, output_dir, batch=True):
    """Fit invariant mass distributions from projection file for each pT bin."""
    if batch:
        ROOT.gROOT.SetBatch(True)

    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    ptbins = config["ptbins"]
    ptmins = ptbins[:-1]
    ptmaxs = ptbins[1:]
    nPtBins = len(ptmins)

    # Fit configuration
    fit_cfg = config.get("v2extraction", config.get("fitConfig", {}))
    sgn_func = fit_cfg.get("SgnFunc", "kGaus")
    bkg_func = fit_cfg.get("BkgFunc", "kExpo")
    mass_ranges = fit_cfg.get("MassFitRanges", [[1.7, 2.05]] * nPtBins)
    sigma_init = fit_cfg.get("Sigma", [0.02] * nPtBins)
    rebin_factors = fit_cfg.get("Rebin", 1)
    if not isinstance(rebin_factors, list):
        rebin_factors = [rebin_factors] * nPtBins
    if not isinstance(sigma_init, list):
        sigma_init = [sigma_init] * nPtBins
    if not isinstance(mass_ranges[0], list):
        mass_ranges = [mass_ranges] * nPtBins

    sgn_name = SGN_MAP.get(sgn_func, "gaussian")
    bkg_name, bkg_pdf_arg = get_bkg_config(bkg_func)

    # Open projection file
    infile = ROOT.TFile.Open(proj_file)
    if not infile or not infile.IsOpen():
        logger(f"Cannot open {proj_file}", "ERROR")
        sys.exit(1)

    # Create output histograms
    pt_arr = np.array(ptbins, dtype="d")
    hRawYields = ROOT.TH1F("hRawYieldsSimFit", ";p_{T} (GeV/c);raw yield", nPtBins, pt_arr)
    hMean = ROOT.TH1F("hMeanSimFit", ";p_{T} (GeV/c);mean", nPtBins, pt_arr)
    hSigma = ROOT.TH1F("hSigmaSimFit", ";p_{T} (GeV/c);#sigma", nPtBins, pt_arr)
    hChi2 = ROOT.TH1F("hRedChi2SimFit", ";p_{T} (GeV/c);#chi^{2}/ndf", nPtBins, pt_arr)
    hSignif = ROOT.TH1F("hRawYieldsSignificanceSimFit", ";p_{T} (GeV/c);significance", nPtBins, pt_arr)
    hSoverB = ROOT.TH1F("hRawYieldsSoverBSimFit", ";p_{T} (GeV/c);S/B", nPtBins, pt_arr)

    for iPt, (pt_min, pt_max) in enumerate(zip(ptmins, ptmaxs)):
        pt_label = f"pt_{int(pt_min*10)}_{int(pt_max*10)}"
        logger(f"Fitting {pt_label}: {pt_min:.1f}–{pt_max:.1f} GeV/c")

        mass_min, mass_max = mass_ranges[min(iPt, len(mass_ranges) - 1)]
        rebin = rebin_factors[min(iPt, len(rebin_factors) - 1)]
        sig0 = sigma_init[min(iPt, len(sigma_init) - 1)]

        # Get mass histogram
        h_mass = infile.Get(f"{pt_label}/hMassData")
        if not h_mass:
            logger(f"  No hMassData for {pt_label}, skipping", "WARNING")
            continue
        h_mass.SetDirectory(0)

        # Fit
        try:
            fitter, ok = try_fit(h_mass, mass_min, mass_max, rebin,
                                 sgn_name, bkg_name, bkg_pdf_arg, sig0)
            if not ok:
                continue
        except Exception as e:
            logger(f"  Fit failed: {e}", "WARNING")
            continue

        # Extract results
        try:
            # Raw yield: tuple (value, error)
            ry_raw = fitter.get_raw_yield(0)
            ry_val, ry_err = float(ry_raw[0]), float(ry_raw[1])

            # Signal pars: list of dicts
            sgn_pars = fitter.get_signal_pars()[0]
            sgn_pars_uncs = fitter.get_signal_pars_uncs()[0]

            mean_val = float(sgn_pars.get("mu", 1.86))
            mean_unc = float(sgn_pars_uncs.get("mu", 0.001))
            sigma_val = float(sgn_pars.get("sigma", sig0))
            sigma_unc = float(sgn_pars_uncs.get("sigma", 0.001))

            chi2 = fitter.get_chi2()
            ndf = fitter.get_ndf()
            chi2_ndf = float(chi2) / max(float(ndf), 1)

            # Significance and S/B: tuples (value, error)
            signif_val = float(fitter.get_significance(0)[0])
            sob_val = float(fitter.get_signal_over_background(0)[0])

        except Exception as e:
            logger(f"  Error extracting results: {e}", "WARNING")
            continue

        bin_idx = iPt + 1
        hRawYields.SetBinContent(bin_idx, ry_val)
        hRawYields.SetBinError(bin_idx, ry_err)
        hMean.SetBinContent(bin_idx, mean_val)
        hMean.SetBinError(bin_idx, mean_unc)
        hSigma.SetBinContent(bin_idx, sigma_val)
        hSigma.SetBinError(bin_idx, sigma_unc)
        hChi2.SetBinContent(bin_idx, max(chi2_ndf, 0))
        hSignif.SetBinContent(bin_idx, signif_val)
        hSoverB.SetBinContent(bin_idx, sob_val)

        logger(f"  → yield={ry_val:.1f}±{ry_err:.1f}, µ={mean_val:.4f}±{mean_unc:.4f}, "
               f"σ={sigma_val:.4f}±{sigma_unc:.4f}, χ²/ndf={chi2_ndf:.2f}")

    infile.Close()

    # Save output
    os.makedirs(output_dir, exist_ok=True)
    proj_basename = os.path.basename(proj_file).replace("proj_", "raw_yields_")
    out_path = os.path.join(output_dir, proj_basename)
    outfile = ROOT.TFile(out_path, "RECREATE")

    hRawYields.Write()
    hMean.Write()
    hSigma.Write()
    hChi2.Write()
    hSignif.Write()
    hSoverB.Write()

    outfile.Close()
    logger(f"Results saved to {out_path}")
    return out_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mass fit from projection files")
    parser.add_argument("config", help="Flow configuration YAML")
    parser.add_argument("proj_file", help="Projection ROOT file")
    parser.add_argument("-o", "--output-dir", default="", help="Output directory")
    parser.add_argument("-b", "--batch", action="store_true", default=True, help="Batch mode")
    args = parser.parse_args()

    out_dir = args.output_dir
    if not out_dir:
        out_dir = os.path.join(os.path.dirname(os.path.dirname(args.proj_file)), "raw_yields")

    run_mass_fit(args.config, args.proj_file, out_dir, args.batch)
