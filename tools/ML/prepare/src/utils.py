import psutil
import time
import yaml
import uproot
import pandas as pd
import os
import gc
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from multiprocessing import Pool
from ROOT import TFile, RDataFrame
import logging
def memory_limit(using_memory=8):
    '''Wait until memory usage is below a certain threshold to avoid memory issues during processing.'''
    # Ensure using_memory is within a safe range
    if using_memory < 1 or using_memory > 30:
        raise ValueError("using_memory should be between 1 and 30")
    waited_seconds = 0
    sleep_to_load = 0
    sleep_to_threshold = 0
    threshold_3 = max(0, 100 - 3 * using_memory)
    threshold_2 = max(0, 100 - 2 * using_memory)
    threshold_1 = max(0, 100 - using_memory)
    while True:
        mem = psutil.virtual_memory()
        if mem.percent < threshold_3:
            if mem.percent < 55:
                sleep_to_load = 0.1
                break
            elif mem.percent < 70:
                sleep_to_load = 0.2
                break
            elif mem.percent < 85:
                sleep_to_load = 0.4
                break
            else:
                time.sleep(1)
                waited_seconds += 1
                print(f"Memory usage above {mem.percent}%, waiting for memory to be freed... {waited_seconds} seconds", end='\r')
        elif mem.percent < threshold_2:
            time.sleep(1)
            waited_seconds += 1
            sleep_to_threshold = 0.6
            print(f"Memory usage {mem.percent}% above {threshold_2}%, waiting for memory to be freed... {waited_seconds} seconds", end='\r')
        else:
            time.sleep(1)
            waited_seconds += 1
            sleep_to_threshold = 1
            print(f"Memory usage {mem.percent}% above {threshold_1}%, waiting for memory to be freed... {waited_seconds} seconds", end='\r')
    print(f"Memory usage is now {mem.percent}%, proceeding with processing. Waited {waited_seconds} seconds.\n")
    return sleep_to_threshold if sleep_to_threshold > 0 else sleep_to_load

def _mutli_args_check(*args):
    """Check if all arguments have the same length."""
    for iArg, arg in enumerate(args):
        if len(arg) != len(args[0]):
            raise ValueError(
                f"Argument {iArg} has a different length {len(arg)} than the first argument {len(args[0])}."
            )

def multi_thread(func, *args, max_workers=4):
    """Run a function in multiple threads."""
    _mutli_args_check(*args)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures, results = [], []
        futures = [executor.submit(func, *arg) for arg in zip(*args)]
        for future in as_completed(futures):
            try:
                result = future.result()
                if isinstance(result, list):
                    results.extend(result)
                else:
                    results.append(result)
            except Exception as e:
                print(f"Error in thread: {e}")
    return results

def multi_process(func, *args):
    """Run a function in multiple processes using ProcessPoolExecutor."""
    _mutli_args_check(*args)

    args_list = list(zip(*args))
    results = []
    max_workers = 28
    using_memory = 3

    try:
        for iArg in range(0, len(args_list), max_workers):
            star_time = time.time()
            batch_args = args_list[iArg:iArg + max_workers]
            futures = []
            loading_time = 0.1

            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                ini_memory = None
                for iArg_tuple, arg_tuple in enumerate(batch_args):
                    sleep_to_load = memory_limit(using_memory)

                    if iArg_tuple == 1:
                        ini_memory = psutil.virtual_memory().percent

                    future = executor.submit(func, *arg_tuple)
                    
                    if iArg_tuple == 0:
                        result = future.result()
                        results.append(result)
                        running_time = time.time() - star_time
                        loading_time = running_time / max_workers
                        logging.info(f"Running time: {running_time:.2f} seconds for {func.__name__}, loading time: {loading_time:.2f} seconds")
                    elif iArg_tuple == 1:
                        i_mem_sample = 0
                        while i_mem_sample < 30:
                            time.sleep(loading_time)
                            i_mem_sample += 1
                            current_memory = psutil.virtual_memory().percent
                            if ((current_memory - ini_memory) > using_memory):
                                using_memory = current_memory - ini_memory
                        processor = (100 - 3 * using_memory - ini_memory) // using_memory
                        loading_time = running_time / processor
                        logging.info(f"Adjusted loading time: {loading_time:.2f} seconds, using_memory: {using_memory}, processor: {processor}")
                        results.append(future.result())
                    else:
                        futures.append(future)

                    time.sleep(loading_time * (0.85 + sleep_to_load))
                    
                    done_futures = [f for f in futures if f.done()]
                    for f in done_futures:
                        try:
                            results.append(f.result())
                        finally:
                            futures.remove(f)
                            del f
                            gc.collect()

                for future in as_completed(futures):
                    try:
                        results.append(future.result())
                    finally:
                        del future
                        gc.collect()

            one_batch_time = time.time() - star_time
            logging.info(
                f"Processed {len(results)}/{len(args_list)}, batch time: {one_batch_time:.2f} seconds, "
                f"left time: {loading_time * (len(args_list) - len(results)):.2f} seconds"
            )
        return results

    except Exception as e:
        print(f"Error in process pool: {e}")
        raise

# #TODO: delete the finished process tighter threshold
# def multi_process(func, *args):
#     """Run a function in multiple processes."""
#     _mutli_args_check(*args)

#     with Pool(processes=28) as pool:
#         try:
#             args_list = list(zip(*args))
#             results = []
#             using_memory = 8
            
#             for iArg in range(0, len(args_list), 28):
#                 star_time = time.time()
#                 async_results = []
#                 for iArg_tuple, arg_tuple in enumerate(args_list[iArg:iArg + 28]):
#                     sleep_to_load = memory_limit(using_memory)
#                     if iArg_tuple == 1:
#                         ini_memory = psutil.virtual_memory().percent
#                     async_result = pool.apply_async(func, arg_tuple)
#                     if iArg_tuple == 0:
#                         results.append(async_result.get())
#                         running_time = time.time() - star_time
#                         loading_time = running_time / 28
#                         print(f"Loading time: {loading_time:.2f} seconds")
#                     elif iArg_tuple == 1:
#                         i_mem_sample = 0
#                         while i_mem_sample < 30:
#                             time.sleep(loading_time)
#                             i_mem_sample += 1
#                             current_memory = psutil.virtual_memory().percent
#                             if ((current_memory - ini_memory) > using_memory):
#                                 using_memory = current_memory - ini_memory
#                         processor = (100 - 2 * using_memory - ini_memory) // using_memory
#                         loading_time = running_time / processor
#                         print(f"Adjusted loading time: {loading_time:.2f} seconds, using_memory: {using_memory}, processor: {processor}")
#                         results.append(async_result.get())
#                     else:
#                         async_results.append(async_result)
#                     time.sleep(loading_time * (0.75 + sleep_to_load))
#                     for r in async_results:
#                         if r.ready():
#                             results.append(r.get())
#                             r.remove()
#                             r.close()
#                             async_results.remove(r)
#                 results.extend([r.get() for r in async_results])
#                 print(f"Processed {len(results)}/{len(args_list)}, left time: {loading_time * (len(args_list) - len(results)):.2f} seconds")
#             return results
#         except Exception as e:
#             print(f"Error in process pool: {e}")
#             raise
#         finally:
#             pool.close()
#             pool.join()

def get_df_name(path_to_file):
    input_file = TFile(path_to_file)
    DFName = [key.GetName() for key in input_file.GetListOfKeys()]
    print(f"DFName: {DFName}")
    if len(DFName) > 2:
        print("Warning: More than one DF in the file")
        if 'parentFiles' in DFName:
            DFName.remove('parentFiles')
        else:
            print("Error: More than one DF in the file and no 'parentFiles' found")
            exit()
    if len(DFName) == 0:
        print("Error: No DF found in the file")
        exit()
    input_file.Close()
    del input_file
    return DFName[0]

def merge_dataframes(*args):
    """
    Merge dataframes from multiple ROOT files.
    
    Parameters:
        *args: Tuple containing input file path and output path:
            - input_file (str) : Path to the input ROOT file.
            - outputPath (str) : Path to save the merged output file (leave empty for same directory).
    Returns:
        output_file (str): Path to the merged output file.
    """
    input_file, outputPath = args
    output_file = os.path.join(os.path.dirname(input_file), f'temp_{os.path.basename(input_file)}'.replace('.root', '_DFmerged.root')) \
                    if outputPath == "" else os.path.join(outputPath, f'temp_{os.path.basename(input_file)}'.replace('.root', '_DFmerged.root'))
    input_txt = os.path.join(os.path.dirname(output_file), "input.txt")
    os.makedirs(os.path.dirname(input_txt), exist_ok=True)
    
    with open(input_txt, 'w') as f:
        f.write(input_file + '\n')
    
    cmd = f"o2-aod-merger --input {input_txt} --output {output_file} --max-size 10000000000"
    os.system(cmd)
    os.remove(input_txt)
    return output_file

def merge_tables(*args):

    dfMerged_file, isMC = args
    print(f"Merging tables from {dfMerged_file} with isMC={isMC}")
    df_name = get_df_name(dfMerged_file)
    input_file = uproot.open(dfMerged_file)
    table_merged = input_file.get(df_name)

    print('Loading tables...')

    try:
        if isMC:
            df_pt_eta_phi_mass_y = table_merged["O2hfd0base"].arrays(library="pd")
            df_mcmatchrec_origin = table_merged["O2hfd0mc"].arrays(library="pd")
            df_pars = table_merged["O2hfd0par"].arrays(library="pd")
            # df_cts = table_merged["O2hfd0pare"].arrays(library="pd")
            df_selflag = table_merged["O2hfd0sel"].arrays(library="pd")

            df_merged = pd.concat([
                df_pt_eta_phi_mass_y,
                df_mcmatchrec_origin,
                df_pars,
                # df_cts,
                df_selflag
            ], axis=1)
        else:
            df_pt_eta_phi_mass_y = table_merged["O2hfd0base"].arrays(library="pd")
            df_pars = table_merged["O2hfd0par"].arrays(library="pd")
            df_selflag = table_merged["O2hfd0sel"].arrays(library="pd")

            df_merged = pd.concat([
                df_pt_eta_phi_mass_y,
                df_pars,
                df_selflag
            ], axis=1)

        output_file = dfMerged_file.replace('_DFmerged.root', '_tableMerged.root')
        df_name = "TreeForML"
        print(f"Writing merged table to {output_file} with df_name={df_name}")

        with uproot.recreate(output_file) as root_file:
            root_file.mktree(df_name, df_merged.dtypes.to_dict())
            root_file[df_name].extend(df_merged.to_dict(orient='list'))

    finally:
        del root_file
        input_file.close()
        del input_file, table_merged

        for var in [
            "df_pt_eta_phi_mass_y", "df_mcmatchrec_origin", "df_pars", 
            "df_cts", "df_selflag", "df_merged"
        ]:
            if var in locals():
                del locals()[var]

        gc.collect()

    return output_file

def prepare_samples(*args):
    """
    Prepare samples for ML from ROOT files.
    Parameters:
        *args: Tuple containing configuration dictionary, input file path, output path, and isMC flag:
            - config (dict): Configuration dictionary with parameters for sample preparation.
                - ptMins (list): List of minimum pT values for each sample.
                - ptMaxs (list): List of maximum pT values for each sample.
                - low_edges (list): List of lower edges for invariant mass cuts, empty if not used.
                - high_edges (list): List of upper edges for invariant mass cuts, empty if not used.
                - ll_edges (list): List of left lower edges for invariant mass cuts.
                - hh_edges (list): List of right upper edges for invariant mass cuts.
                - data_filters (list): List of filters for data samples, empty if not used.
                - data_snapshot (str): variable names for data snapshot.
                - mc_prompt_filters (list): List of filters for prompt MC samples, empty if not used.
                - mc_prompt_snapshot (str): variable names for prompt MC snapshot.
                - mc_fd_filters (list): List of filters for FD MC samples, empty if not used.
                - mc_fd_snapshot (str): variable names for FD MC snapshot.
            - inputFile (str): Path to the input ROOT file.
            - outputPath (str): Path to save the prepared samples (leave empty for same directory).
            - isMC (bool): Flag indicating if the samples are from MC (True) or data (False).
    Returns:
        outputPath (str): Path to the directory where the prepared samples are saved.
    """
    inputFile, config, outputPath, isMC = args
    
    ptMins = config['ptMins']
    ptMaxs = config['ptMaxs']

    low_edges = config.get('low_edges', [])
    high_edges = config.get('high_edges', [])
    ll_edges = config['ll_edges']
    hh_edges = config['hh_edges']

    data_filters = config.get('data_filters', [])
    data_snapshot = config.get('data_snapshot')
    mc_prompt_filters = config.get('mc_prompt_filters', [])
    mc_prompt_snapshot = config.get('mc_prompt_snapshot')
    mc_fd_filters = config.get('mc_fd_filters', [])
    mc_fd_snapshot = config.get('mc_fd_snapshot')
    
    if outputPath == "":
        outputPath = os.path.dirname(inputFile)
    os.makedirs(outputPath, exist_ok=True)
    suffix = f'{min(ptMins)}_{max(ptMaxs)}'

    # join the inv mass filters
    if len(low_edges) == 0 and len(high_edges) == 0:
        pt_mass_cut = " || ".join([f"({ptMins[i]} < fPt && fPt < {ptMaxs[i]} && ({ll_edges[i]} < fM && fM < {hh_edges[i]}))" \
                                    for i in range(len(ptMins))])
    else:
        pt_mass_cut = " || ".join([f"({ptMins[i]} < fPt && fPt < {ptMaxs[i]} && (({ll_edges[i]} < fM && fM < {low_edges[i]}) || ({high_edges[i]} < fM && fM < {hh_edges[i]})))" \
                                    for i in range(len(ptMins))])
    print(f'pt_mass_cut: {pt_mass_cut}')

    # join the filters
    if data_filters:
        data_filter = " && ".join(data_filters)
        print('data_filter: ', data_filter)
    if mc_prompt_filters:
        mc_prompt_filter = " && ".join(mc_prompt_filters)
        print('mc_prompt_filter: ', mc_prompt_filter)
    if mc_fd_filters:
        mc_fd_filter = " && ".join(mc_fd_filters)
        print('mc_fd_filter: ', mc_fd_filter)

    treeFile = TFile(inputFile, "read")
    dataTree, mcTree, dfDataForMLApply, dfMcPromptForApply, dfMcFDForApply = None, None, None, None, None
    print(mc_fd_snapshot, mc_prompt_snapshot)
    if isMC:
        mcTree = treeFile.Get(config['treeName'])
        if mc_prompt_filters:
            logging.info(f'Processing prompt MC {inputFile} with tree {config["treeName"]}')
            dfMcPromptForApply = RDataFrame(mcTree)
            dfMcPromptForApply.Filter(mc_prompt_filter + " && " + pt_mass_cut) \
                .Snapshot("TreeForML", f'{outputPath}/McTreeForMLPromptTrain_{suffix}.root', mc_prompt_snapshot)
        if mc_fd_filters:
            logging.info(f'Processing FD MC {inputFile} with tree {config["treeName"]}')
            dfMcFDForApply = RDataFrame(mcTree)
            dfMcFDForApply.Filter(mc_fd_filter + " && " + pt_mass_cut) \
                .Snapshot("TreeForML", f'{outputPath}/McTreeForMLFDTrain_{suffix}.root', mc_fd_snapshot)
    else:
        dataTree = treeFile.Get(config['treeName'])
        logging.info(f'Processing data {inputFile} with tree {config["treeName"]}') 
        dfDataForMLApply = RDataFrame(dataTree)
        dfDataForMLApply.Filter(data_filter + " && " + pt_mass_cut) \
            .Snapshot("TreeForML", f'{outputPath}/DataTreeForMLTrain_{suffix}.root', data_snapshot)

    treeFile.Close()
    del treeFile, dataTree, mcTree, dfDataForMLApply, dfMcPromptForApply, dfMcFDForApply
    gc.collect()

    return outputPath

def prepare_multi_samples(config, inputFiles, outputPath="", isMC=False):
    if not isinstance(inputFiles, list):
        inputFiles = [inputFiles]
    if isinstance(config, dict):
        config = [config] * len(inputFiles)
        

# def PrepareSamples_multi(config):
#     with open(config, 'r') as file:
#         config = yaml.safe_load(file)

#     inputDataName = config.get("inputDataFile", '')
#     inputMcName = config.get("inputMcFile", '')
    
#     if inputDataName != "":
#         configs = []
#         inputDataFiles = []
#         for root, _, files in os.walk(inputDataName):
#             found_files = [os.path.join(root, name) for name in files if name.endswith('_tableMerged.root')]
#             inputDataFiles.extend(found_files)
#         if len(inputDataFiles) == 0:
#             print("No data files found.")
#         else:
#             print(f"Found {len(inputDataFiles)} data files.")
#             for inputDataFile in inputDataFiles:
#                 print(f"Processing data file: {inputDataFile}")
#                 config['inputDataFile'] = inputDataFile
#                 config['inputMcFile'] = ""
#                 config['outDir'] = os.path.dirname(inputDataFile)
#                 configs.append(config.copy())
#             for i in range(0, len(configs), 5):
#                 chunk = configs[i:i+5]
#                 multiprocess(PrepareSamples, [chunk], max_workers=5)
#                 print(f"Processed chunk {i//5 + 1}/{(len(inputDataFiles) + 5 - 1) // 5}")

#     if inputMcName != "":
#         configs = []
#         inputMcFiles = []
#         for root, _, files in os.walk(inputMcName):
#             found_files = [os.path.join(root, name) for name in files if name.endswith('_tableMerged.root')]
#             inputMcFiles.extend(found_files)
#         if len(inputMcFiles) == 0:
#             print("No MC files found.")
#         else:
#             print(f"Found {len(inputMcFiles)} MC files.")
#             for inputMcFile in inputMcFiles:
#                 print(f"Processing MC file: {inputMcFile}")
#                 config['inputMcFile'] = inputMcFile
#                 config['inputDataFile'] = ""
#                 config['outDir'] = os.path.dirname(inputMcFile)
#                 configs.append(config.copy())
#             for i in range(0, len(configs), 5):
#                 chunk = configs[i:i+5]
#                 multiprocess(PrepareSamples, [chunk], max_workers=5)