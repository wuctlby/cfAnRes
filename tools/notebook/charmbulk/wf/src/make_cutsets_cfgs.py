#!/usr/bin/env python3
"""
make_cutsets_cfgs.py — Generate cutset YAML files from the flow configuration.

Correlated mode: all pT bins share the same BDT cut index.
Output: {outdir}/cutsets/cutset_00.yml, cutset_01.yml, ...

Usage:
    python3 make_cutsets_cfgs.py flow_config.yml -o output_dir [--correlated]
"""
import argparse
import os
import yaml
import numpy as np
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from utils import logger


def pad_to_length(lst, target_len):
    """Pad a list to target length by repeating the last element."""
    lst_len = len(lst)
    if lst_len >= target_len:
        return lst
    return lst + [lst[-1]] * (target_len - lst_len)


def make_cutsets(config_file, output_dir, correlated=True):
    """Generate cutset YAML files for cut variation analysis."""
    with open(config_file, "r") as f:
        cfg = yaml.safe_load(f)

    ptmins = cfg["ptbins"][:-1]
    ptmaxs = cfg["ptbins"][1:]
    nPtBins = len(ptmins)

    cfg_cutvar = cfg["cut_variation"]

    if correlated:
        sig = cfg_cutvar["corr_bdt_cut"]["sig"]
        sig_cuts_lower = [list(np.arange(sig["min"][i], sig["max"][i], sig["step"][i])) for i in range(min(nPtBins, len(sig["min"])))]
        sig_cuts_upper = [[1.0] * len(cuts) for cuts in sig_cuts_lower]
    else:
        sig = cfg_cutvar["uncorr_bdt_cut"]["sig"]
        sig_cuts_lower = [sig[i][:-1] for i in range(min(nPtBins, len(sig)))]
        sig_cuts_upper = [sig[i][1:] for i in range(min(nPtBins, len(sig)))]

    # Pad all to uniform length
    maxCutSets = max(len(cuts) for cuts in sig_cuts_lower) if sig_cuts_lower else 0
    if maxCutSets == 0:
        logger("No cuts found!", "ERROR")
        return

    sig_cuts_lower = [pad_to_length(cuts, maxCutSets) for cuts in sig_cuts_lower]
    sig_cuts_upper = [pad_to_length(cuts, maxCutSets) for cuts in sig_cuts_upper]

    # Transpose: [iPt][iCut] → [iCut][iPt]
    sig_cuts_lower = list(map(list, zip(*sig_cuts_lower)))
    sig_cuts_upper = list(map(list, zip(*sig_cuts_upper)))

    if correlated:
        bkg_cuts_upper = [cfg_cutvar["corr_bdt_cut"]["bkg_max"][:nPtBins]] * maxCutSets
    else:
        bkg_data = cfg_cutvar["uncorr_bdt_cut"]["bkg_max"]
        bkg_cuts_upper = [pad_to_length(b, maxCutSets) for b in bkg_data[:nPtBins]]
        bkg_cuts_upper = list(map(list, zip(*bkg_cuts_upper)))

    os.makedirs(f"{output_dir}/cutsets", exist_ok=True)
    for iCut, (bkg_maxs, fd_mins, fd_maxs) in enumerate(zip(bkg_cuts_upper, sig_cuts_lower, sig_cuts_upper)):
        bkg_max = list(map(float, bkg_maxs[:nPtBins]))
        fd_min = list(map(float, fd_mins[:nPtBins]))
        fd_max = list(map(float, fd_maxs[:nPtBins]))

        combinations = {
            "icutset": iCut,
            "Pt": {"min": list(ptmins), "max": list(ptmaxs)},
            "ScoreBkg": {"min": [0.0] * nPtBins, "max": bkg_max},
            "ScoreFD": {"min": fd_min, "max": fd_max},
        }

        out_path = f"{output_dir}/cutsets/cutset_{iCut:02d}.yml"
        with open(out_path, "w") as f:
            yaml.dump(combinations, f, default_flow_style=False, sort_keys=False)

    logger(f"Created {maxCutSets} cutset files in {output_dir}/cutsets")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate cutset YAML files")
    parser.add_argument("config", help="Flow configuration YAML")
    parser.add_argument("-o", "--output-dir", default=".", help="Output directory")
    parser.add_argument("--correlated", "-c", action="store_true", help="Correlated cuts")
    args = parser.parse_args()

    make_cutsets(args.config, args.output_dir, args.correlated)
