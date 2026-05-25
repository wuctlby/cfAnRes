"""Minimal utilities for the charm bulk workflow."""

import os
import shutil
import sys
import ROOT
import numpy as np
from ROOT import TH1F, TFile

ROOT.TH1.AddDirectory(False)

work_dir = os.path.dirname(os.path.realpath(__file__))

_SCRIPTS = {
	"Preprocess": 'pre_process.py',
	"YamlCuts": 'pre_process.py',
	"CutVariation": 'pre_process.py',
	"DataDrivenFraction": 'data_driven_fraction.py',
}

def _get_script_name(script_path):
    script_path = os.path.basename(script_path)  # Get the filename from the path
    if script_path not in _SCRIPTS.values():
        raise ValueError(f"Script path '{script_path}' not found. Available keys: {list(_SCRIPTS.values())}")
    return f"[{next(key for key, value in _SCRIPTS.items() if value == script_path)}]"

def make_dir_root_file(directory, file, verbose=True):
    if not file.GetDirectory(directory):
        file.mkdir(directory)
        if verbose:
            logger(f"Created directory {directory} in file {file.GetName()}", level='WARNING')
    else:
        if verbose:
            logger(f"Directory {directory} already exists in file {file.GetName()}", level='WARNING')

def check_dir(dir_path, clean=True, script=None):
    if os.path.exists(dir_path):
        if clean:
            logger(f"{dir_path} already exists, it will be overwritten", level="WARNING", script=script)
            shutil.rmtree(dir_path)
            os.makedirs(dir_path)
        return
    os.makedirs(dir_path)

def logger(message, level="INFO", script=None):
    msg = f"[{level}] {message}"
    if script:
        message = f"{_get_script_name(script)} {msg}"
    if level == "ERROR":
        print(f"\033[31m{msg}\033[0m")
    elif level == "WARNING":
        print(f"\033[33m{msg}\033[0m")
    elif level == "COMMAND":
        print(f"\033[34m{msg}\033[0m")
    else:
        print(f"\033[32m{msg}\033[0m")

def list_files(input_path, prefix="", suffix="", exfix="",
               subdir="", exdir="", script=None):
    """List files under *directory* with flexible filtering.

    Parameters
    ----------
    input_path : str
        Base input_path, could be a directory, a text file containing paths, or a list of file paths.
    prefix : str
        If non empty, only files whose name **starts with** *prefix*.
    suffix : str
        If non empty, only files whose name **ends with** *suffix*.
    exfix : str
        If non-empty, **exclude** files whose name **contains** *exfix*.
    subdir : str
        If non-empty, only descend into sub-directories whose
        **full path** contains *subdir*.
    exdir : str
        If non-empty, **skip** sub-directories whose **full path**
        contains *exdir*.
    script : str
        Optional script path for logging context.

    Returns
    -------
    list[str]
        Sorted list of matching absolute file paths.

    Examples
    --------
    >>> list_files("/data/projs", prefix="proj_")
    >>> list_files("/data/effs",  suffix=".root")
    >>> list_files("/data",       prefix="eff_", subdir="cutvar")
    >>> list_files("/data",       prefix="raw_yields_", exdir="uncorr")
    """
    from natsort import natsorted
    from pathlib import Path
    
    def _walk_directory(directory, prefix, suffix, exfix, subdir, exdir):
        if not os.path.isdir(directory):
            raise NotADirectoryError(f"Not a directory: {directory}")

        matches: list[str] = []
        for root, dirs, filenames in os.walk(directory):
            if subdir and subdir not in root:
                continue
            if exdir and exdir in root:
                continue
            for filename in filenames:
                if prefix and not filename.startswith(prefix):
                    continue
                if suffix and not filename.endswith(suffix):
                    continue
                if exfix and exfix in filename:
                    continue
                matches.append(os.path.join(root, filename))
        return natsorted(matches)

    file_paths: list[str] = []
    if (isinstance(input_path, str) or isinstance(input_path, Path)) and not str(input_path).endswith(".txt"):
    # path
        file_paths = _walk_directory(input_path, prefix, suffix, exfix, subdir, exdir)
    elif (isinstance(input_path, str) or isinstance(input_path, Path)) and str(input_path).endswith(".txt"):
    # txt file with paths
        with open(input_path, 'r') as f:
            for line in f:
                path = line.strip()
                if path and not path.startswith("#"):  # skip empty lines and comments
                    if os.path.isfile(path):
                        file_paths.append(path)
                    else:
                        file_paths.extend(_walk_directory(path, prefix, suffix, exfix, subdir, exdir))
    elif isinstance(input_path, list):
    # list of paths
        file_paths = input_path
    else:
        logger(f"Invalid type for {input_path} in configuration. Must be a string (directory or text file) or list of file paths.", "ERROR", script=script)
        sys.exit(1)
    return natsorted(file_paths)


def load_root_files(inputPath, prefix, suffix=".root"):
    """Load root files from a directory matching prefix and suffix."""
    if not os.path.exists(inputPath):
        raise ValueError(f"No folder found in {inputPath}")
    files = sorted(
        f for f in os.listdir(inputPath)
        if f.startswith(prefix) and f.endswith(suffix)
    )
    return [os.path.join(inputPath, f) for f in files]


def load_eff_histos(effFiles):
    """Load efficiency histograms from a file or list of files."""
    def _load_single(f):
        fin = TFile.Open(f)
        hEp = fin.Get("hEffPrompt").Clone()
        hEf = fin.Get("hEffFD").Clone()
        hEp.SetDirectory(0)
        hEf.SetDirectory(0)
        hPF = hEp.Clone("hPromptFrac")
        hFF = hEf.Clone("hFDFrac")
        hPFc = hEp.Clone("hPromptFracCorr")
        hFFc = hEf.Clone("hFDFracCorr")
        fin.Close()
        return hEp, hEf, hPF, hFF, hPFc, hFFc

    if isinstance(effFiles, str):
        return _load_single(effFiles)
    lists = [[] for _ in range(6)]
    for f in effFiles:
        for i, h in enumerate(_load_single(f)):
            lists[i].append(h)
    return tuple(lists)


def load_cutVar_histos(cutVarFile):
    """Load histograms from cutVar.root."""
    fin = TFile.Open(cutVarFile)
    hCP = fin.Get("hCorrYieldPrompt")
    hCF = fin.Get("hCorrYieldFD")
    hPP = fin.Get("hCovPromptPrompt")
    hPF = fin.Get("hCovPromptFD")
    hFF = fin.Get("hCovFDFD")
    for h in [hCP, hCF, hPP, hPF, hFF]:
        h.SetDirectory(0)
    fin.Close()
    return hCP, hCF, hPP, hPF, hFF


def GetMinimisation(effPromptList, effFDList, rawYieldList,
                    effPromptUncList, effFDUncList, rawYieldUncList,
                    corr=True, precision=1e-8, nMaxIter=100):
    """Analytic system minimisation to retrieve prompt & FD corrected yields."""
    nCutSets = len(effPromptList)
    mRawYield = np.zeros((nCutSets, 1))
    mEff = np.zeros((nCutSets, 2))
    for i, (ry, ep, ef) in enumerate(zip(rawYieldList, effPromptList, effFDList)):
        mRawYield[i] = ry
        mEff[i, 0] = ep
        mEff[i, 1] = ef
    mRawYield = np.matrix(mRawYield)
    mEff = np.matrix(mEff)

    mCorrYield = np.matrix(np.zeros((2, 1)))
    mCorrYieldOld = np.matrix(np.zeros((2, 1)))

    for iIter in range(nMaxIter):
        if iIter == 0:
            mCorrYield[0] = 0
            mCorrYield[1] = 0
        mCovSets = np.matrix(np.zeros((nCutSets, nCutSets)))
        mCorrSets = np.matrix(np.zeros((nCutSets, nCutSets)))
        for iR, (ruR, epR, efR) in enumerate(zip(rawYieldUncList, effPromptUncList, effFDUncList)):
            for iC, (ruC, epC, efC) in enumerate(zip(rawYieldUncList, effPromptUncList, effFDUncList)):
                uncR = np.sqrt(ruR**2 + epR**2 * float(mCorrYield[0])**2 + efR**2 * float(mCorrYield[1])**2)
                uncC = np.sqrt(ruC**2 + epC**2 * float(mCorrYield[0])**2 + efC**2 * float(mCorrYield[1])**2)
                # Diagonal covariance only for numerical stability
                rho = 1.0 if iR == iC else 0.0
                mCovSets[iR, iC] = rho * uncR * uncC
                mCorrSets[iR, iC] = rho

        # Diagonal — simple inverse
        mWeights = np.linalg.inv(mCovSets)
        mEffT = mEff.T
        mCovariance = np.linalg.inv(mEffT * mWeights * mEff)
        mCorrYield = mCovariance * mEffT * mWeights * mRawYield
        mRes = mEff * mCorrYield - mRawYield
        if (abs(mCorrYield[0] - mCorrYieldOld[0]) / mCorrYield[0] < precision and
                abs(mCorrYield[1] - mCorrYieldOld[1]) / mCorrYield[1] < precision):
            break
        mCorrYieldOld = np.copy(mCorrYield)

    mResT = np.transpose(mRes)
    redChiSquare = float(mResT * mWeights * mRes / (nCutSets - 2))
    dic = {"covMatrix": mCovSets, "weightMatrix": mWeights, "corrMatrix": mCorrSets}
    return mCorrYield, mCovariance, redChiSquare, dic


def GetPromptFDFractionCutSet(accEffPrompt, accEffFD, corrYieldPrompt, corrYieldFD,
                              covPromptPrompt, covFDFD, covPromptFD):
    """Get prompt and FD fractions with uncertainties for a cutset."""
    denom = accEffPrompt * corrYieldPrompt + accEffFD * corrYieldFD
    if denom == 0:
        return [0., 0.], [0., 0.]

    fP = accEffPrompt * corrYieldPrompt / denom
    dPdeNP = (accEffPrompt * (accEffPrompt * corrYieldPrompt + accEffFD * corrYieldFD)
              - accEffPrompt**2 * corrYieldPrompt) / (accEffPrompt * corrYieldPrompt + accEffFD * corrYieldFD)**2
    dPdeNF = -accEffFD * accEffPrompt * corrYieldPrompt / (accEffPrompt * corrYieldPrompt + accEffFD * corrYieldFD)**2
    fPUnc = np.sqrt(dPdeNP**2 * covPromptPrompt + dPdeNF**2 * covFDFD + 2 * dPdeNP * dPdeNF * covPromptFD)

    fFD = accEffFD * corrYieldFD / (accEffPrompt * corrYieldPrompt + accEffFD * corrYieldFD)
    dFdeNF = (accEffFD * (accEffFD * corrYieldFD + accEffPrompt * corrYieldPrompt)
              - accEffFD**2 * corrYieldFD) / (accEffPrompt * corrYieldPrompt + accEffFD * corrYieldFD)**2
    dFdeNP = -accEffFD * accEffPrompt * corrYieldFD / (accEffPrompt * corrYieldPrompt + accEffFD * corrYieldFD)**2
    fFDUnc = np.sqrt(dFdeNF**2 * covFDFD + dFdeNP**2 * covPromptPrompt + 2 * dFdeNF * dFdeNP * covPromptFD)

    return [fP, fFD], [fPUnc, fFDUnc]


def get_centrality_bins(centrality):
    '''
    Get centrality bins

    Input:
        - centrality:
            str, centrality class (e.g. 'k3050')

    Output:
        - cent_bins:
            list of floats, centrality bins
        - cent_label:
            str, centrality label
    '''
    if centrality == 'k05':
        return '0_5', [0, 5]
    if centrality == 'k510':
        return '5_10', [5, 10]
    if centrality == 'k010':
        return '0_10', [0, 10]
    if centrality == 'k1015':
        return '10_15', [10, 15]
    if centrality == 'k1520':
        return '15_20', [15, 20]
    if centrality == 'k1020':
        return '10_20', [10, 20]
    if centrality == 'k1030':
        return '10_30', [10, 30]
    if centrality == 'k020':
        return '0_20', [0, 20]
    if centrality == 'k1030':
        return '10_30', [10, 30]
    if centrality == 'k2030':
        return '20_30', [20, 30]
    elif centrality == 'k2050':
        return '20_50', [20, 50]
    elif centrality == 'k3040':
        return '30_40', [30, 40]
    elif centrality == 'k3050':
        return '30_50', [30, 50]
    elif centrality == 'k4050':
        return '40_50', [40, 50]
    elif centrality == 'k2060':
        return '20_60', [20, 60]
    elif centrality == 'k4060':
        return '40_60', [40, 60]
    elif centrality == 'k4080':
        return '40_80', [40, 80]
    elif centrality == 'k5060':
        return '50_60', [50, 60]
    elif centrality == 'k5080':
        return '50_80', [50, 80]
    elif centrality == 'k50100':
        return '50_100', [50, 100]
    elif centrality == 'k6070':
        return '60_70', [60, 70]
    elif centrality == 'k6080':
        return '60_80', [60, 80]
    elif centrality == 'k7080':
        return '70_80', [70, 80]
    elif centrality == 'k0100':
        return '0_100', [0, 100]
    else:
        print(f"ERROR: cent class \\'{centrality}\\' is not supported! Exit")
    sys.exit()


# ── histogram reweighting helpers ─────────────────────────────────────────────

def reweight_histo_1D(hist, spline, binned=False):
    """Reweight a 1‑D histogram using a spline or binned weights."""
    if binned:
        for i in range(1, hist.GetNbinsX() + 1):
            hist.SetBinContent(i, hist.GetBinContent(i) * spline.GetBinContent(i))
        return hist
    for i in range(1, hist.GetNbinsX() + 1):
        pt = hist.GetBinCenter(i)
        w = spline(pt)
        hist.SetBinContent(i, hist.GetBinContent(i) * w)
    return hist


def reweight_histo_2D(hist, spline, binned=False):
    """Reweight a 2‑D histogram using a spline or binned weights on the X axis."""
    for i in range(1, hist.GetNbinsX() + 1):
        pt = hist.GetXaxis().GetBinCenter(i)
        w = spline(pt) if not binned else spline.GetBinContent(i)
        for j in range(1, hist.GetNbinsY() + 1):
            hist.SetBinContent(i, j, hist.GetBinContent(i, j) * w)
    return hist


def reweight_histo_3D(hist, spline, species_weights):
    """Reweight a 3‑D histogram with a spline and per‑species weights."""
    for i in range(1, hist.GetNbinsX() + 1):
        pt = hist.GetXaxis().GetBinCenter(i)
        w_pt = spline(pt)
        for j in range(1, hist.GetNbinsY() + 1):
            for k in range(1, hist.GetNbinsZ() + 1):
                w_species = species_weights[j - 1]  # simplified; adjust as needed
                hist.SetBinContent(i, j, k, hist.GetBinContent(i, j, k) * w_pt * w_species)
    return hist


# ── vn vs mass helpers ────────────────────────────────────────────────────────

def get_vn_versus_mass(sparse, mass_bins, mass_axis, sp_axis):
    """Compute vn vs invariant mass from a THnSparse (scalar‑product method)."""
    import ROOT, array
    hist = ROOT.TH1F("hVnVsMass", ";#it{m}_{inv} (GeV/#it{c}^{2});#it{v}_{2}",
                     len(mass_bins) - 1, array.array("d", mass_bins))
    for i_bin in range(1, hist.GetNbinsX() + 1):
        m_low = mass_bins[i_bin - 1]
        m_high = mass_bins[i_bin]
        sparse.GetAxis(mass_axis).SetRangeUser(m_low, m_high)
        hsp = sparse.Projection(sp_axis)
        if hsp.Integral() > 0:
            vn = hsp.GetMean()
            err = hsp.GetMeanError()
        else:
            vn = err = 0
        hist.SetBinContent(i_bin, vn)
        hist.SetBinError(i_bin, err)
    return hist


def profile_mass_sp(sparse, mass_bins, resolution=1.0):
    """Profile (Sp) as a function of mass, used for multitrial systematics."""
    import ROOT, array
    hist = ROOT.TProfile("hVnVsMass", ";#it{m}_{inv} (GeV/#it{c}^{2});#it{v}_{2}",
                         len(mass_bins) - 1, array.array("d", mass_bins))
    for i_bin in range(1, hist.GetNbinsX() + 1):
        m_low = mass_bins[i_bin - 1]
        m_high = mass_bins[i_bin]
        sparse.GetAxis(0).SetRangeUser(m_low, m_high)
        hsp = sparse.Projection(1)
        if hsp.Integral() > 0:
            mean = hsp.GetMean() / resolution
            err = hsp.GetMeanError() / resolution
        else:
            mean = err = 0
        hist.Fill((m_low + m_high) / 2, mean)
    return hist