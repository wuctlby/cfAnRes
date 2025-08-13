import argparse
import yaml
from copy import deepcopy
import os
from datetime import datetime
import sys
sys.path.append('/home/wuct/ALICE/local/RTools/cfAnRes/src')
sys.path.append('/home/wuct/ALICE/local/RTools/cfAnRes/tools/ML/prepare/src')
from file_loader import FileLoader
from utils import multi_thread, merge_dataframes, multi_process, prepare_samples, merge_tables
import logging

def setup_logging(log_dir="/data/meta/BDT/DataPass4_MCPass4/Data/logs", log_level=logging.INFO):
    os.makedirs(log_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"run_{timestamp}.log")

    logging.basicConfig(
        level=log_level,
        format='[%(asctime)s] [%(levelname)s] %(funcName)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[
            logging.StreamHandler(),             # Console
            logging.FileHandler(log_file, mode='a')  # File
        ]
    )
    logging.info(f"Logging initialized. Writing to {log_file}")

def sample_config(config):
    ptSampleConfigs = []
    for ptmin, ptmax, low_edges, high_edges, ll_edges, hh_edges in zip(
        config.get("ptMins"), 
        config.get("ptMaxs"), 
        config.get("low_edges", [[]] * len(config.get("ptMins"))),
        config.get("high_edges", [[]] * len(config.get("ptMins"))),
        config.get("ll_edges"),
        config.get("hh_edges")
    ):
        ptSampleConfigs.append(deepcopy(config))
        ptSampleConfigs[-1]['ptMins'] = [ptmin]
        ptSampleConfigs[-1]['ptMaxs'] = [ptmax]
        ptSampleConfigs[-1]['low_edges'] = [low_edges] if low_edges else []
        ptSampleConfigs[-1]['high_edges'] = [high_edges] if high_edges else []
        ptSampleConfigs[-1]['ll_edges'] = [ll_edges]
        ptSampleConfigs[-1]['hh_edges'] = [hh_edges]
    return ptSampleConfigs

def process(config_sample, data=False, mc=False):
    
    with open(config_sample, 'r') as file:
        cfg = yaml.safe_load(file)
    if data:
        config = cfg.get("Data")
    elif mc:
        config = cfg.get("MC")

    stepFiles = []
    for iOpt, operator in enumerate(config.get("operations")['name']):
        if operator == "df_merge":
            if stepFiles == []:
                stepFiles = sorted(FileLoader(config.get("inputFiles"), keyFileName='AO2D.root').paths)

            dfMergedFiles = multi_thread(merge_dataframes, stepFiles, ['']*len(stepFiles), 
                                         max_workers=config.get("operations").get("maxWorkers")[iOpt])
            stepFiles = sorted(dfMergedFiles)

        elif operator == "table_merge":
            if stepFiles == []:
                stepFiles = sorted(FileLoader(config.get("inputFiles"), keyFileName='DFmerged.root').paths)
            
            tableMergedFiles = multi_process(merge_tables, stepFiles, [mc]*len(stepFiles))
            
            stepFiles = sorted(tableMergedFiles)

        elif operator == "prepare_samples_pt":
            print(config.get("prepare_samples"))
            ptSampleConfigs = sample_config(config.get("prepare_samples"))
            if stepFiles == []:
                stepFiles = sorted(FileLoader(config.get("inputFiles"), keyFileName='_tableMerged.root').paths)
            else:
                if len(stepFiles) != len(FileLoader(config.get("inputFiles"), keyFileName='_tableMerged.root').paths):
                    logging.warning("Input files for prepare_samples_pt do not match the expected number of files. "
                                    "Please check your input configuration.")
                    continue
                else:
                    for file in stepFiles:
                        if not os.path.exists(file):
                            logging.warning(f"File {file} does not exist. Skipping this file.")
                            stepFiles.remove(file)
            sampledFiles = []
            print(stepFiles)
            print(ptSampleConfigs)
            for ptConfig in ptSampleConfigs:
                sampledFiles = multi_process(prepare_samples, stepFiles, [ptConfig]*len(stepFiles), ['']*len(stepFiles), [mc]*len(stepFiles))

            stepFiles = sorted(sampledFiles)
        elif operator == "file_merge_pt":
            os.makedirs(f"{config.get('outputPath')}", exist_ok=True)
            for ptMin, ptMax in zip(config.get("file_merge").get("ptMins"), config.get("file_merge").get("ptMaxs")):
                if data:
                    preparedFiles = sorted(FileLoader(config.get("inputFiles"), keyFileName=f'_{ptMin}_{ptMax}.root').paths)
                    outputFile = os.path.join(config.get("outputPath"), f'DataTreeForMLTrain_{ptMin}_{ptMax}.root')
                    command = f"hadd -f -k -j 16 -n 2 {outputFile} " + " ".join(preparedFiles)
                    os.system(command)
                if mc:
                    preparedFiles = sorted(FileLoader(config.get("inputFiles"), keyFileName=f'_{ptMin}_{ptMax}.root').paths)
                    promptInputFile = [f for f in preparedFiles if 'Prompt' in f]
                    promptOutputFile = os.path.join(config.get("outputPath"), f'McTreeForMLPromptTrain_{ptMin}_{ptMax}.root')
                    command = f"hadd -f -k -j 16 -n 2 {promptOutputFile} " + " ".join(promptInputFile)
                    os.system(command)
                    fdInputFile = [f for f in preparedFiles if 'FD' in f]
                    fdOutputFile = os.path.join(config.get("outputPath"), f'McTreeForMLFDTrain_{ptMin}_{ptMax}.root')
                    command = f"hadd -f -k -j 16 -n 2 {fdOutputFile} " + " ".join(fdInputFile)
                    os.system(command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare Samples for ML")
    parser.add_argument('config_sample', metavar='text', default='config_sample.yml',
                        help='the configuration file for preparing samples')
    parser.add_argument('--data', '-d', action='store_true', default=False,
                        help='prepare data samples')
    parser.add_argument('--mc', '-m', action='store_true', default=False,
                        help='prepare MC samples')
    args = parser.parse_args()

    setup_logging()

    if args.data and args.mc:
        print("Please specify either --data or --mc, not both.")
        sys.exit(1)
    if args.data:
        logging.info("Preparing data samples...")
        process(args.config_sample, data=True, mc=False)
    if args.mc:
        logging.info("Preparing MC samples...")
        process(args.config_sample, data=False, mc=True)