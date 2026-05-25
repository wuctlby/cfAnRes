'''
This script is used to pre-process a/multi large AnRes.root for the BDT training:
    - split the input by pT
    - obtain the sigma from prompt enhance sample
python3 pre_process.py config.yml AnRes_1.root AnRes_2.root --pre --sigma  
'''
import os
import sys
import yaml
import numpy as np
import array
from ROOT import TFile, TObject
import argparse
import gc
import concurrent.futures
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
SCRIPT = os.path.basename(__file__)
sys.path.append(f"{SCRIPT_DIR}/utils/")
from utils import get_centrality_bins, logger, list_files
from utils_thn import GetTHnInfo
import uproot
import pandas as pd
import awkward as ak

def check_existing_outputs(file_path):
    """
    Check if output file already exists, and open it accordingly.

    Args:
        file_path (str): path to the output ROOT file.

    Returns:
        tuple: (out_file, write_opt) where out_file is the opened TFile and write_opt is the write option.
    """
    if os.path.exists(file_path):
        logger(f"    Updating file: {file_path}", level='WARNING')
        out_file = TFile(file_path, 'update')
        write_opt = TObject.kOverwrite
    else:
        logger(f"    Creating file: {file_path}", level='WARNING')
        out_file = TFile.Open(file_path, 'recreate')
        write_opt = TObject.kSingleKey # Standard

    return out_file, write_opt

def make_dir_root_file(directory, file, verbose=True):
    if not file.GetDirectory(directory):
        file.mkdir(directory)
        if verbose:
            logger(f"Created directory {directory} in file {file.GetName()}", level='WARNING', script=SCRIPT)
    else:
        if verbose:
            logger(f"Directory {directory} already exists in file {file.GetName()}", level='WARNING', script=SCRIPT)

def get_inputs_sparse(file, full_cfg, thn_info, sparse_name, debug=False):
    """Load a single sparse and axes info
    
    Args:
        file (TFile): input ROOT file
        full_cfg (dict): full configuration dictionary
        thn_info (THnSparseInfo): THnSparseInfo object
        sparse_name (str): name of the sparse to load
        debug (bool, optional): print debug info. Defaults to False.
    
    Returns:
        tuple: (sparse, axes) where sparse is the loaded sparse histogram and axes is the dictionary of axes information
    """
    sparse = file.Get(thn_info.thurl)
    if sparse is None:
        logger(f"Sparse {thn_info.thname} not found in file {file.GetName()} at path {thn_info.thpath}", level='ERROR', script=SCRIPT)
    else:
        logger(f"Sparse {sparse} loaded from {thn_info.thpath}", level='INFO', script=SCRIPT)

    if full_cfg['Dmeson'] == 'Dzero' and thn_info.datatype == 'MC':
        # TODO: safety checks for Dmeson reflecton and secondary peak
        if sparse_name == "RecoPrompt":
            sparse.GetAxis(thn_info.axis_id('Origin')).SetRange(2, 2)       # select prompt
            sparse.GetAxis(thn_info.axis_id('CandType')).SetRange(1, 2)     # select signal
        elif sparse_name == "RecoFD":
            sparse.GetAxis(thn_info.axis_id('Origin')).SetRange(3, 3)       # select non-prompt
            sparse.GetAxis(thn_info.axis_id('CandType')).SetRange(1, 2)     # select signal
        elif sparse_name == "RecoRefl":
            sparse.GetAxis(thn_info.axis_id('CandType')).SetRange(3, 4)     # select reflection
        elif sparse_name == "RecoReflPrompt":
            sparse.GetAxis(thn_info.axis_id('CandType')).SetRange(3, 4)     # select reflection
            sparse.GetAxis(thn_info.axis_id('Origin')).SetRange(2, 2)       # select prompt
        elif sparse_name == "RecoReflFD":
            sparse.GetAxis(thn_info.axis_id('CandType')).SetRange(3, 4)     # select reflection
            sparse.GetAxis(thn_info.axis_id('Origin')).SetRange(3, 3)       # select FD
        elif sparse_name == "GenPrompt":
            sparse.GetAxis(thn_info.axis_id('Origin')).SetRange(2, 2)       # select prompt
        elif sparse_name == "GenFD":
            sparse.GetAxis(thn_info.axis_id('Origin')).SetRange(3, 3)       # select non-prompt
        else:
            logger(f"Unknown sparse type for Dzero {thn_info.thname}", level='ERROR')

    if debug:
        print('\n')
        logger('###############################################################', level='DEBUG')
        for key, value in thn_info.axis_name_id_map.items():
            logger(f"{key}: {value}", level='DEBUG')
        logger('###############################################################\n', level='DEBUG')

    return sparse

def process_sparse(i_file, infile_path, full_cfg, sparse_name, prep_out_dir):
    """
    Process a single sparse from an input file for all pt bins according to the configuration.
    
    Args:
        i_file (int): index of the input file
        infile (TFile): input ROOT file
        full_cfg (dict): full configuration dictionary
        sparse_name (str): name of the sparse to process
        prep_out_dir (str): output directory for pre-processed files
    """
    infile = TFile.Open(infile_path, 'read')
    thn_info = GetTHnInfo.thn(sparse_name, cand=full_cfg['Dmeson'])
    logger(f'Processing sparse {sparse_name} for file {i_file}', level='INFO', script=SCRIPT)
    sparse = get_inputs_sparse(infile, full_cfg, thn_info, sparse_name, debug=False)
    if sparse is None:
        logger(f"Sparse {sparse_name} not found in file {infile_path} at path {thn_info.thurl}", level='ERROR', script=SCRIPT)
        infile.Close()
        return
    else:
        logger(f"Sparse {sparse} loaded from {thn_info.thurl}", level='INFO', script=SCRIPT)

    # Apply centrality cut if centrality axis exists
    if thn_info.axis_id('cent') is not None:
        cent_min, cent_max = get_centrality_bins(full_cfg.get('centrality', 'k0100'))[1] # Default to 0-100% if not specified
        logger(f"Applying cent cut to sparse {sparse} with value {cent_min} -- {cent_max}", "INFO", script=SCRIPT)
        sparse.GetAxis(thn_info.axis_id('cent')).SetRangeUser(cent_min, cent_max)

    pt_mins, pt_maxs = full_cfg['ptbins'][:-1], full_cfg['ptbins'][1:]
    bkg_maxs = full_cfg['preprocess']['bkg_cuts']
    axes_to_keep = [thn_info.axis_name(i) for i in thn_info.axisids_kept]
    type_name = thn_info.datatype

    logger(f"Projecting sparse {sparse} for file {i_file} with axes to keep {axes_to_keep}", level='INFO', script=SCRIPT)
    for pt_min, pt_max, bkg_max in zip(pt_mins, pt_maxs, bkg_maxs):
        logger(f"Processing pT bin {pt_min} - {pt_max} with bkg max {bkg_max}", level='INFO', script=SCRIPT)
        # Create output file
        out_file_dir = f"{prep_out_dir}/preprocess/pt_{int(pt_min*10)}_{int(pt_max*10)}/{type_name}/jobs"
        os.makedirs(out_file_dir, exist_ok=True)
        out_file_path = f'{out_file_dir}/AnalysisResults_{i_file}.root'
        out_file, write_opt = check_existing_outputs(out_file_path)

        sparse.GetAxis(thn_info.axis_id('pt')).SetRangeUser(pt_min, pt_max) # PtTrig for correlations, Pt for SP flow
        if thn_info.axis_id('bkg score') is not None: # Skip sparses for generated info
            sparse.GetAxis(thn_info.axis_id('bkg score')).SetRangeUser(0, bkg_max)
        proj_axes = thn_info.axisids_kept
        proj_sparse = sparse.Projection(len(proj_axes), array.array('i', proj_axes), 'O')
        proj_sparse.SetName(thn_info.pre_thname)
        # proj_sparse = proj_sparse.Rebin(array.array('i', rebin))
        make_dir_root_file(thn_info.pre_thpath, out_file)
        out_file.cd(thn_info.pre_thpath)
        proj_sparse.Write(thn_info.pre_thname, write_opt)
        out_file.Delete(thn_info.pre_thname + ";*")
        proj_sparse.Delete()
        del proj_sparse

        gc.collect()

        out_file.Close()
        logger(f"Saved projected sparse to {out_file_path}", level='INFO', script=SCRIPT)

    del sparse
    gc.collect()
    infile.Close()

def process_tree(i_file, infile_path, full_cfg, trees_cfg, prep_out_dir, input_out_dir):
    logger(f'[Data] Processing tree file {i_file}, {infile_path}', "INFO")
    cols_dict, columns_to_keep = {}, []

    # Build columns dictionary and columns to keep
    for tree_cfg in trees_cfg['trees']:
        cols_dict.update(get_tree_dict(tree_cfg['table']))
        columns_to_keep.extend([cols_dict[col] for col in tree_cfg['cols']])

    try:
        def collect_trees(uproot_dir, tree_names, trees):
            for key in uproot_dir.keys():
                obj = uproot_dir[key]
                if '/' in key: # Else one gets duplicate entries
                    continue
                if isinstance(obj, uproot.TTree) and obj.name in tree_names:
                    # extend existing awkward array with new data
                    trees[obj.name] = ak.concatenate([trees[obj.name], obj.arrays(library="ak")])
                elif isinstance(obj, uproot.reading.ReadOnlyDirectory):
                    collect_trees(obj, tree_names, trees)  # recurse into subdirectory

        trees = {tree_cfg['table']: ak.Array([]) for tree_cfg in trees_cfg['trees']}
        tree_names = [t['table'] for t in trees_cfg['trees']]
        with uproot.open(infile_path) as f:
            collect_trees(f, tree_names, trees)

    except Exception as e:
        print(f"Failed to open file {infile_path}: {e}")
        return

    try:
        merged_tree = ak.zip(
            {field: trees[tree][field]
            for tree in trees 
            for field in trees[tree].fields}
        )
    except Exception as e:
        import traceback
        traceback.print_exc()
        logger(f"Caught Python exception while merging trees: {e}", "ERROR")
        return

    # Loop over pt bins
    pt_mins, pt_maxs = full_cfg['ptbins'][:-1], full_cfg['ptbins'][1:]
    bkg_maxs = full_cfg['preprocess']['bkg_cuts']

    try:
        for pt_min, pt_max, bkg_max in zip(pt_mins, pt_maxs, bkg_maxs):
            out_file_dir = f"{prep_out_dir}/preprocess/pt_{int(pt_min*10)}_{int(pt_max*10)}/{input_out_dir}"
            os.makedirs(out_file_dir, exist_ok=True)
            out_file_path = f'{out_file_dir}/AO2D_{i_file}.root'

            # Selection (centrality is optional)
            if 'Cent' in cols_dict:
                cent_min, cent_max = get_centrality_bins(full_cfg['centrality'])[1]
                mask = (
                    (merged_tree[cols_dict['Pt']] >= pt_min) &
                    (merged_tree[cols_dict['Pt']] <= pt_max) &
                    (merged_tree[cols_dict['Cent']] >= cent_min) &
                    (merged_tree[cols_dict['Cent']] <= cent_max) &
                    (merged_tree[cols_dict['ScoreBkg']] >= 0) &
                    (merged_tree[cols_dict['ScoreBkg']] <= bkg_max)
                )
            else:
                mask = (
                    (merged_tree[cols_dict['Pt']] >= pt_min) &
                    (merged_tree[cols_dict['Pt']] <= pt_max) &
                    (merged_tree[cols_dict['ScoreBkg']] >= 0) &
                    (merged_tree[cols_dict['ScoreBkg']] <= bkg_max)
            )
            merged_tree_sel = merged_tree[mask]

            # Save snapshot
            uproot.recreate(out_file_path)[trees_cfg['final_table_name']] = {
                col: merged_tree_sel[col] for col in columns_to_keep
            }
            logger(f"Saved snapshot to {out_file_path}", "INFO")
    except Exception as e:
        import traceback
        traceback.print_exc()
        logger(f"Caught Python exception while saving snapshot: {e}", "ERROR")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument('config', metavar='text', default='config.yml', help='configuration file')
    parser.add_argument("--workers", "-w", type=int, default=1, help="number of workers")
    args = parser.parse_args()

    with open(args.config, 'r') as cfg_pre:
        full_cfg = yaml.safe_load(cfg_pre)

    output_dir = full_cfg['preprocess']['outdir'] if full_cfg.get('preprocess').get('outdir') else full_cfg['outdir']
    for input_cfg in full_cfg['preprocess']['inputs']:
        files = input_cfg['files']

        thn_info = GetTHnInfo.thn(input_cfg['sparses'][0])
        type_name = thn_info.datatype

        if input_cfg.get('sparses'):
            file_paths = list_files(files, prefix="AnalysisResults", suffix=".root", script=SCRIPT_DIR)
            for sparse_name in input_cfg['sparses']:
                logger(f"##### Skimming {sparse_name} #####", "INFO")
                with concurrent.futures.ThreadPoolExecutor(args.workers) as executor:
                    tasks_sparses = [executor.submit(process_sparse, i_file, file, full_cfg, sparse_name, output_dir) for i_file, file in enumerate(file_paths)]
                    for task in tasks_sparses:
                        task.result()

        elif input_cfg.get('trees'):
            logger(f"##### Skimming {input_cfg} trees #####", "INFO")
            file_paths = list_files(files, prefix="AO2D", suffix=".root", script=SCRIPT_DIR)
            with concurrent.futures.ThreadPoolExecutor(args.workers) as executor:
                tasks_trees = [executor.submit(process_tree, i_file, file_path, full_cfg, input_cfg, output_dir, f"{input_cfg['outdir']}/jobs") for i_file, file_path in enumerate(file_paths)]

        else:
            logger("No sparses or trees found in the configuration for pre-processing", "ERROR")

        # Dump file_paths used
        dump_file_paths = f"{output_dir}/preprocess/file_paths_{type_name}.txt"
        os.makedirs(os.path.dirname(dump_file_paths), exist_ok=True)
        with open(dump_file_paths, "w") as dump_file:
            for file_path in file_paths:
                dump_file.write(f"{file_path}\n")

        logger(f"Finished processing {type_name}\n\n", "INFO")

    # use hadd to merge the results in the directories for the different pT bins
    for pt_dir in os.listdir(f"{output_dir}/preprocess/"):
        if not os.path.isdir(f"{output_dir}/preprocess/{pt_dir}"):
            continue
        for input_type_dir in os.listdir(f"{output_dir}/preprocess/{pt_dir}"):
            prep_dir = f"{output_dir}/preprocess/{pt_dir}/{input_type_dir}"
            files_to_merge = [f"./jobs/{file}" for file in os.listdir(f"{prep_dir}/jobs") if file.endswith(".root")]
            files_to_merge_str = ' '.join(files_to_merge)
            file_name = os.path.basename(files_to_merge[0]).split('_')[0]
            merged_file = f"{prep_dir}/{file_name}_{pt_dir}.root"
            log_merge = f"{prep_dir}/log_merge.txt"
            filelist_path = f"{prep_dir}/files_to_hadd.txt"
            with open(filelist_path, "w") as f: # To avoid crash if too many files
                for fpath in files_to_merge:
                    f.write(fpath + "\n")
            hadd = f"hadd -f {merged_file} @{filelist_path} > {log_merge}"
            logger(f"Running command: {hadd}\n\n", "INFO")
            try:
                os.system(f"cd {prep_dir} && {hadd}")
            except Exception as e:
                logger(f"Error while creating hadd command: {e}", "ERROR")
                os.system(f"rm {merged_file}") # Remove the partially merged file
                os.system(f"mv {log_merge} {prep_dir}/log_merge_error.txt")
