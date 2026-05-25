"""
Projection utilities for the charm bulk workflow.

Extracted from :file:`proj_thn.py` — contains all ``proj_xx`` functions and
the unified :func:`get_pt_preprocessed_sparses` which handles **pre=True**
(preprocessed files) with **pre=False** fallback to original files.
"""
import os
import glob
import sys
import yaml
import numpy as np
from pathlib import Path
from functools import partial
from concurrent.futures import ThreadPoolExecutor
from scipy.interpolate import make_interp_spline

import ROOT
from ROOT import TFile, TObject, TH1F

from utils_thn import GetTHnInfo
from utils import (
    reweight_histo_1D, reweight_histo_2D, reweight_histo_3D,
    get_vn_versus_mass, profile_mass_sp,
    make_dir_root_file, logger, get_centrality_bins, list_files
)
sys.path.append(f"{os.path.dirname(os.path.abspath(__file__))}/..")
from pre_process import get_inputs_sparse

ROOT.TH1.AddDirectory(False)

# ── pt weight helper ─────────────────────────────────────────────────────────
def get_pt_weights(cfgProj):
    """Get pt weights and return weights flags with spline

    Args:
        cfgProj (dict): Configuration dictionary for projections

    Outputs:
        sPtWeights (spline): Spline for ptWeights interpolation
        sPtWeightsB (spline): Spline for ptWeightsB weights interpolation
        Bspeciesweights (str): B species weights # TODO
    """

    # REVIEW: the ptWeights inputed is a list, but the ptWeights outputed is a TH1D object
    # and actually ptweights is used as a flag
        # compute info for pt weights
    if not cfgProj.get('PtWeightsFile'):
        logger('No pt weights for D and B mesons provided in the config file!', level='WARNING')
        return None, None, None
        
    ptWeightsFile = TFile.Open(cfgProj["PtWeightsFile"], 'r')

    if cfgProj.get('ApplyPtWeightsD'):
        hPtWeightsD = ptWeightsFile.Get('hPtWeightsFONLLtimesTAMUDcent')
        ptBinCentersD = [ (hPtWeightsD.GetBinLowEdge(i)+hPtWeightsD.GetBinLowEdge(i+1))/2 for i in range(1, hPtWeightsD.GetNbinsX()+1)]
        ptBinContentsD = [hPtWeightsD.GetBinContent(i) for i in range(1, hPtWeightsD.GetNbinsX()+1)]
        sPtWeights = make_interp_spline(ptBinCentersD, ptBinContentsD)
    else:
        logger('pt weights for D mesons will not be provided!', level='WARNING')
        sPtWeights = None

    if cfgProj.get('ApplyPtWeightsB'):
        hPtWeightsB = ptWeightsFile.Get('hPtWeightsFONLLtimesTAMUBcent')
        ptBinCentersB = [ (hPtWeightsB.GetBinLowEdge(i)+hPtWeightsB.GetBinLowEdge(i+1))/2 for i in range(1, hPtWeightsB.GetNbinsX()+1)]
        ptBinContentsB = [hPtWeightsB.GetBinContent(i) for i in range(1, hPtWeightsB.GetNbinsX()+1)]
        sPtWeightsB = make_interp_spline(ptBinCentersB, ptBinContentsB)
    else:
        logger('pt weights for B mesons will not be provided!', level='WARNING')
        sPtWeightsB = None

    if cfgProj.get('ApplyBSpeciesWeights'):
        Bspeciesweights = cfgProj['Bspeciesweights']
    else:
        logger('B species weights will not be provided!', level='WARNING')
        Bspeciesweights = None
    
    return sPtWeights, sPtWeightsB, Bspeciesweights

# ── resolution helper ─────────────────────────────────────────────────────────
def get_resolution(config, operations, multitrial_folder=False):
    if operations.get("proj_data") or multitrial_folder:
        reso_file = TFile.Open(config["projections"]["Resolution"], 'r')
        det_A = config["projections"].get('detA', 'FT0c')
        det_B = config["projections"].get('detB', 'FV0a')
        det_C = config["projections"].get('detC', 'TPCtot')
        logger(f"Getting resolution histogram from file {config['projections']['Resolution']} for triplet {det_A}_{det_B}_{det_C}",  "WARNING")
        reso_hist = reso_file.Get(f'{det_A}_{det_B}_{det_C}/histo_reso_delta_cent')
        resolution = reso_hist.GetBinContent(1)
        reso_hist.SetDirectory(0)
        reso_file.Close()
    return reso_hist,resolution

# ── multitrial projection ────────────────────────────────────────────────────
def proj_multitrial(config, multitrial_folder, workers, operations):

    if operations.get("proj_data", False):
         reso_hist, resolution = get_resolution(config, operations, multitrial_folder=True)

    pt_bin_label = Path(multitrial_folder).name
    logger(f"Processing multitrial projections for pt bin {pt_bin_label} ...", level='INFO')

    # Load default cutsets
    cutset_dir = f"{config['outdir']}/cutvar_{config['suffix']}_combined/cutsets"
    default_cutsets = [f"{cutset_dir}/{f}" for f in os.listdir(cutset_dir) if f.endswith('.yml')]
    # Load Mass and MassSp histos from the default cases
    default_histos = {}
    for default_cutset in default_cutsets:
        suffix = os.path.basename(default_cutset).replace(".yml", "").replace("cutset_", "")
        default_proj = TFile.Open(default_cutset.replace(".yml", ".root").replace("cutset", "proj"), "READ")
        default_histos[suffix] = {}
        default_histos[suffix]['Mass'] = default_proj.Get(f"{pt_bin_label}/hMassData")
        default_histos[suffix]['MassSp'] = default_proj.Get(f"{pt_bin_label}/hMassSpData")
        default_proj.Close()

    def process_cutset(multitrial_dir, default_histos):
        trial_number = Path(multitrial_dir).name.replace("trial_", "")
        try:
            with open(f"{multitrial_dir}/config_trial_{trial_number}.yml", 'r') as ymlCutSetFile:
                config_trial = yaml.safe_load(ymlCutSetFile)
        except Exception as e:
            logger(f"Error opening or reading config file for trial {trial_number}: {e}", level='ERROR')
            return

        for suffix, histo in default_histos.items():
            output_dir = f"{multitrial_dir}/projs"
            os.makedirs(output_dir, exist_ok=True)
            output_path = f"{output_dir}/proj_{suffix}.root"
            output_file = TFile.Open(output_path, "RECREATE")
            output_file.mkdir(pt_bin_label)
            output_file.cd(pt_bin_label)
            default_histos[suffix]['Mass'].Write("hMassData")
            hist_vn_vs_mass = profile_mass_sp(default_histos[suffix]['MassSp'], config_trial['projections']['inv_mass_bins'][0], resolution)
            hist_vn_vs_mass.Write("hVnVsMassData")
            output_file.Close()

        logger(f"[{trial_number}] Completed projections!", level='INFO')

    # Parallel execution
    multitrial_dirs = [f for f in glob.glob(f"{multitrial_folder}/trials/*") if os.path.isdir(f)]
    with ThreadPoolExecutor(max_workers=workers) as executor:
        executor.map(partial(process_cutset, default_histos=default_histos), multitrial_dirs)


# ── loading sparses ───────────────────────────────────────────────────────
def _load_preprocessed_sparse(pre_dir, pt_label, thn_info):
    """Try to load a THnSparse from the preprocessed file."""
    pre_file_path = f"{pre_dir}/preprocess/{pt_label}/{thn_info.datatype}/AnalysisResults_{pt_label}.root"
    if not os.path.exists(pre_file_path):
        return None, None, False
    pre_file = TFile.Open(pre_file_path, "READ")
    if not pre_file or pre_file.IsZombie():
        return None, None, False
    axes = thn_info.axis_name_id_map
    sparse = pre_file.Get(f"{thn_info.pre_thpath}/{thn_info.pre_thname}")
    pre_file.Close()
    return sparse, axes, not isinstance(sparse, ROOT.TObject)

def _load_original_sparse(origin_files, thn_info, Dmeson=None, cent=None, pt_label=None):
    """Load a THnSparse from the original file."""
    file_paths = list_files(origin_files, prefix="AnalysisResults", suffix=".root")

    def _process_sparse(i_file, infile_path, Dmeson, cent, pt_label):
        infile = TFile.Open(infile_path, 'read')
        logger(f'Processing sparse {thn_info.thname} for file {i_file}', level='INFO')
        sparse = get_inputs_sparse(infile, Dmeson, thn_info, thn_info.thname, debug=False)
        if sparse is None:
            logger(f"Sparse {thn_info.thname} not found in file {infile_path} at path {thn_info.thurl}", level='ERROR')
            infile.Close()
            return
        else:
            logger(f"Sparse {thn_info.thname} loaded from {thn_info.thurl}", level='INFO')

        # Apply centrality cut if centrality axis exists
        if thn_info.axis_id('cent') is not None:
            cent_min, cent_max = get_centrality_bins(cent)[1] # Default to 0-100% if not specified
            logger(f"Applying cent cut to sparse {thn_info.thname} with value {cent_min} -- {cent_max}", "INFO")
            sparse.GetAxis(thn_info.axis_id('cent')).SetRangeUser(cent_min, cent_max)
        # Apply pt cut if pt axis exists
        if thn_info.axis_id('pt') is not None and pt_label is not None:
            pt_min, pt_max = map(float, pt_label.replace("pt_", "").split("_"))
            logger(f"Applying pt cut to sparse {thn_info.thname} with value {pt_min} -- {pt_max}", "INFO")
            sparse.GetAxis(thn_info.axis_id('pt')).SetRangeUser(pt_min, pt_max)
        infile.Close()

        return sparse

    import concurrent
    sparses = []
    with concurrent.futures.ThreadPoolExecutor(4) as executor:
        tasks_sparses = [executor.submit(_process_sparse, i_file, file, Dmeson, cent, pt_label) for i_file, file in enumerate(file_paths)]
        for task in tasks_sparses:
            sparses.append(task.result())

    merged = sparses[0].Clone()
    for sparse in sparses[1:]:
        merged.Add(sparse)
    return merged, thn_info.axis_name_id_map, not isinstance(sparse, ROOT.TObject)

def get_pt_preprocessed_sparses(config, pt_label):
    """ Load preprocessed sparses for a given pt bin label from the preprocess output files.
    Args:
        config (dict): Configuration dictionary.
        pt_label (str): Pt bin label to load sparses for.
    Returns:
        sparses_data (dict): Dictionary of data sparses.
        sparses_reco (dict): Dictionary of reconstructed MC sparses.
        sparses_gen (dict): Dictionary of generated MC sparses.
        axes (dict): Dictionary of axes for each sparse.
    """
    from utils_thn import GetTHnInfo
    sparses_data, sparses_reco, sparses_gen, axes = {}, {}, {}, {}
    pre_cfg = config['preprocess']

    # Find preprocess config of sparse
    thns_proj_data_info = []
    thns_proj_mc_info = []
    for input_cfg in pre_cfg['inputs']:
        # print(f"input_cfg: {input_cfg}")
        if not 'sparses' in input_cfg:
            continue
        for sparse_cfg in input_cfg['sparses']:
            thn_info = GetTHnInfo.thn(sparse_cfg['name'], pre=True, cand=config.get("Dmeson", "D0"))
            if thn_info.datatype == "Data" and sparse_cfg.get("proj_data", False):
                thns_proj_data_info.append(thn_info)
            elif thn_info.datatype == "MC" and sparse_cfg.get("proj_mc", False):
                thns_proj_mc_info.append(thn_info)

    prep_dir = config.get("preprocess", {}).get("outdir") or config["outdir"] + "/preprocess"
    # input_files = {sparse_cfg['name']: input_cfg['files'] for input_cfg in pre_cfg['inputs'] if 'sparses' in input_cfg for sparse_cfg in input_cfg['sparses']}
    input_files = {}
    for input_cfg in pre_cfg['inputs']:
        for sparse_cfg in input_cfg.get('sparses', []):
            input_files[sparse_cfg['name']] = input_cfg['files']

    if config["operations"].get("proj_data") or config["operations"].get("perform_fits_syst_reso"):
        for thn_proj_data_info in thns_proj_data_info:
            sparse, axes, success = _load_preprocessed_sparse(prep_dir, pt_label, thn_proj_data_info)
            if not success:
                original_files = input_files.get(thn_proj_data_info.name, [])
                if not original_files:
                    logger(f"No input files specified for {thn_proj_data_info.name} in config!", level='ERROR')
                    continue
                sparse, axes, success = _load_original_sparse(original_files, thn_proj_data_info, Dmeson=config.get("Dmeson"), cent=config.get("centrality", "k0100"), pt_label=pt_label)
            sparses_data[thn_proj_data_info.name] = sparse
            axes[thn_proj_data_info.name] = axes

    if config["operations"].get("proj_mc"):
        for thn_proj_mc_info in thns_proj_mc_info:
            sparse, axes, success = _load_preprocessed_sparse(prep_dir, pt_label, thn_proj_mc_info)
            if not success:
                original_files = input_files.get(thn_proj_mc_info.name, [])
                if not original_files:
                    logger(f"No input files specified for {thn_proj_mc_info.name} in config!", level='ERROR')
                    continue
                sparse, axes, success = _load_original_sparse(original_files, thn_proj_mc_info, Dmeson=config.get("Dmeson"), cent=config.get("centrality", "k0100"), pt_label=pt_label)
            if thn_proj_mc_info.name.startswith("Reco"):
                sparses_reco[thn_proj_mc_info.name] = sparse
            elif thn_proj_mc_info.name.startswith("Gen"):
                sparses_gen[thn_proj_mc_info.name] = sparse
            axes[thn_proj_mc_info.name] = axes

    return sparses_data, sparses_reco, sparses_gen, axes

# ── mc projections ─────────────────────────────────────────────────────────
def proj_mc_reco(sparses_reco, sPtWeightsD, sPtWeightsB, Bspeciesweights, writeopt, axes, pt_min, pt_max, save_centrality=False):

    for key, sparse in sparses_reco.items():
        if key != 'RecoPrompt' and key != 'RecoFD':
            sparse.Projection(axes[key]['Mass']).Write(f'h{key}Mass')
            sparse.Projection(axes[key]['Pt']).Write(f'h{key}Pt')

    hMassPrompt = sparses_reco['RecoPrompt'].Projection(axes['RecoPrompt']['Mass'])
    hMassPrompt.SetName(f'hPromptMass_{pt_min}_{pt_max}')
    hMassFD = sparses_reco['RecoFD'].Projection(axes['RecoFD']['Mass'])
    hMassFD.SetName(f'hFDMass_{pt_min}_{pt_max}')

    ### project pt prompt
    hPtPrompt = sparses_reco['RecoPrompt'].Projection(axes['RecoPrompt']['Pt'])
    if sPtWeightsD:
        hPtPrompt = reweight_histo_1D(hPtPrompt, sPtWeightsD, binned=False)

    ### project pt FD
    if sPtWeightsD:
        hPtFD = reweight_histo_1D(sparses_reco['RecoFD'].Projection(axes['RecoFD']['Pt']), sPtWeightsD, binned=False)
    elif sPtWeightsB:
        if Bspeciesweights:
            hPtFD = reweight_histo_3D(
                sparses_reco['RecoFD'].Projection(axes['RecoFD']['Pt'], axes['RecoFD']['PtBMoth'], axes['RecoFD']['FlagBHad']),
                sPtWeightsB, Bspeciesweights
            )
        else:
            hPtFD = reweight_histo_2D(
                sparses_reco['RecoFD'].Projection(axes['RecoFD']['PtBMoth'], axes['RecoFD']['Pt']),          # 2D projection: Projection(ydim, xdim)
                sPtWeightsB, binned=False
            )
    elif Bspeciesweights:
        hPtFD = reweight_histo_2D(
            sparses_reco['RecoFD'].Projection(axes['RecoFD']['FlagBHad'], axes['RecoFD']['Pt']),             # 2D projection: Projection(ydim, xdim)
            Bspeciesweights, binned=True
        )
    else:
        hPtFD = sparses_reco['RecoFD'].Projection(axes['RecoFD']['Pt'])

    ## write the output
    hMassPrompt.Write('hPromptMass', writeopt)
    hMassFD.Write('hFDMass', writeopt)
    hPtPrompt.Write('hPromptPt', writeopt)
    hPtFD.Write('hFDPt', writeopt)

    # Store also centrality of prompt, if available
    if save_centrality and 'Cent' in axes['RecoPrompt']:
        hRecoCentPrompt = sparses_reco['RecoPrompt'].Projection(axes['RecoPrompt']['Cent'])
        hRecoCentPrompt.Write('hPromptRecoCent', writeopt)

    return hPtPrompt, hPtFD

def proj_mc_gen(sparses_gen, sPtWeightsD, sPtWeightsB, Bspeciesweights, writeopt, axes, pt_min, pt_max, save_centrality=False):

    for key, sparse in sparses_gen.items():
        if key != 'GenPrompt' and key != 'GenFD':
            sparse.Projection(axes[key]['Pt']).Write(f'h{key}Pt')

    ### prompt
    hGenPtPrompt = sparses_gen['GenPrompt'].Projection(axes['GenPrompt']['Pt'])
    if sPtWeightsD:
        hGenPtPrompt = reweight_histo_1D(hGenPtPrompt, sPtWeightsD, binned=False)

    ### FD
    if sPtWeightsD:
        hGenPtFD = reweight_histo_1D(sparses_gen['GenFD'].Projection(axes['GenFD']['Pt']), sPtWeightsD, binned=False)
    elif sPtWeightsB:
        if Bspeciesweights:
            hGenPtFD = reweight_histo_3D(
                sparses_gen['GenFD'].Projection(axes['GenFD']['Pt'], axes['GenFD']['PtBMoth'], axes['GenFD']['FlagBHad']),
                sPtWeightsB, Bspeciesweights
            )
        else:
            hGenPtFD = reweight_histo_2D(
                sparses_gen['GenFD'].Projection(axes['GenFD']['PtBMoth'], axes['GenFD']['Pt']),         # 2D projection: Projection(ydim, xdim)
                sPtWeightsB, binned=False
            )
    elif Bspeciesweights:
        hGenPtFD = reweight_histo_2D(
            sparses_gen['GenFD'].Projection(axes['GenFD']['FlagBHad'], axes['GenFD']['Pt']),            # 2D projection: Projection(ydim, xdim)
            Bspeciesweights, binned=True
        )
    else:
        hGenPtFD = sparses_gen['GenFD'].Projection(axes['GenFD']['Pt'])

    ## write the output
    hGenPtPrompt.Write('hPromptGenPt', writeopt)
    hGenPtFD.Write('hFDGenPt', writeopt)

    # Store also centrality of prompt, if available
    if save_centrality and 'Cent' in axes['GenPrompt']:
        hGenCentPrompt = sparses_gen['GenPrompt'].Projection(axes['GenPrompt']['Cent'])
        hGenCentPrompt.Write('hPromptGenCent', writeopt)

    return hGenPtPrompt, hGenPtFD

# ── data projections ─────────────────────────────────────────────────────────
def proj_data(i_bin, process, sparse, axes, proj_cfg, writeopt):

    if process == "proj_data":
         proj_vars = proj_cfg.get('ProjVars', [])
         proj_vars += ['Mass']
         proj_axes = [axes['Data'][var] for var in proj_vars]
    elif process == "proj_data_charmbulk":
        proj_vars = proj_cfg.get('ProjVars', [])
        proj_vars += ['Mass']
        proj_axes = [axes['CharmBulk'][var] for var in proj_vars]
        if proj_cfg.get('y_bins'):
            y_low, y_high = proj_cfg['y_bins']
            sparse.GetAxis(axes['CharmBulk']['y']).SetRangeUser(y_low, y_high)
    elif process == "proj_data_sp":
        proj_vars = proj_cfg.get('ProjVars', [])
        proj_vars += ['Mass', 'Sp']
        proj_axes = [axes['FlowSP'][var] for var in proj_vars]

    for var, ax in zip(proj_vars, proj_axes):
        hist_var = sparse.Projection(ax)
        hist_var.Write(f'h{var.capitalize()}Data', writeopt)

    if process == "proj_data_charmbulk":
        pass
    elif process == "proj_data_sp":
        hist_vn_sp = get_vn_versus_mass(sparse, proj_cfg['inv_mass_bins'][i_bin], axes['FlowSP']['Mass'], axes['FlowSP']['Sp'])
        _, resolution = get_resolution(proj_cfg, {"proj_data": True})
        hist_vn_sp.Scale(1/resolution) # Correct for resolution
        hist_vn_sp.Write('hVnVsMassData', writeopt)

        # Save a TH2 of (Mass, Sp) for the multitrial systematic
        mass_lowest_bin, mass_highest_bin, sp_lowest_bin, sp_highest_bin = -1, 1e10, -1, 1e10
        hist_mass_sp = sparse.Projection(axes['FlowSP']['Sp'], axes['FlowSP']['Mass'])
        sp_lowest_bin = hist_mass_sp.ProjectionY().FindFirstBinAbove(0)
        sp_highest_bin = hist_mass_sp.ProjectionY().FindLastBinAbove(0)
        mass_lowest_bin = hist_mass_sp.ProjectionX().FindFirstBinAbove(0)
        mass_highest_bin = hist_mass_sp.ProjectionX().FindLastBinAbove(0)
        hist_mass_sp.GetXaxis().SetRange(mass_lowest_bin, mass_highest_bin)
        hist_mass_sp.GetYaxis().SetRange(sp_lowest_bin, sp_highest_bin)
        hist_mass_sp.Write('hMassSpData', writeopt)