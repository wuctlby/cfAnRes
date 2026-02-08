import os
import re
from numpy import size
import yaml
import argparse
from alive_progress import alive_bar
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor
from itertools import groupby

def filter_staged_files(staged_files):
    """Filter staged files to keep only those from Stage_(max-1)."""
    stage_lines = []
    no_stage_lines = []
    stage_pattern = re.compile(r"Stage_(\d+)")

    for line in staged_files:
        match = stage_pattern.search(line)
        if match:
            stage_number = int(match.group(1))
            stage_lines.append((stage_number, line))
        else:
            no_stage_lines.append(line)

    if not stage_lines: # no staged files found, return all raw files
        return no_stage_lines

    max_stage = max(stage_lines, key=lambda x: x[0])[0]
    if max_stage == 1: # only Stage_1 exists, return all raw files
        return no_stage_lines
    target_stage = max_stage - 1

    result = [line for num, line in stage_lines if num == target_stage]
    return result

def get_file_list_from_alien(args):
    '''Get the list of files from Alien directory.'''
    alien_dirs, Stage = args
    alien_dirs = [alien_dirs] if not isinstance(alien_dirs, list) else alien_dirs

    file_list = []
    for alien_dir in alien_dirs:
        temp_file_list = []
        RunNumber = alien_dir.strip().split('/')[-1]
        if RunNumber == "AOD": # in case it is slim data
            RunNumber = alien_dir.strip().split('/')[-2]
        if not Stage:
            os.system(f"alien_find {alien_dir} AnalysisResults.root -r > output_{RunNumber}.txt")
        else:
            with open(f"output_{RunNumber}.txt", "w") as f:
                f.write("")

        check_unmerged = False
        with open(f"output_{RunNumber}.txt", "r") as f:
            temp_file_list = [(line.strip(), RunNumber) for line in f if line.strip()]
            if len(temp_file_list) == 0:
                check_unmerged = True

        if check_unmerged:
            os.system(f"alien_find {alien_dir} */AnalysisResults.root > output_{RunNumber}.txt")
            staged_files = []
            with open(f"output_{RunNumber}.txt", "r") as uf:
                for line in uf:
                    staged_files.append(line.strip())

            # remove duplicates in case
            staged_files = list(set(staged_files))

            # keeping only Stage_(max-1)
            temp_file_list.extend([(line, RunNumber) for line in filter_staged_files(staged_files)])

        file_list.extend(temp_file_list)
        os.remove(f"output_{RunNumber}.txt")
    return file_list

def download_file(args):
    '''Download a file from Alien to local path.'''
    iFile, alien_file, RunNumber, LocalPath, FileName = args

    local_dir = os.path.join(LocalPath, RunNumber, f"{iFile:04d}")
    os.makedirs(local_dir, exist_ok=True)
    if FileName == "AO2D.root":
        os.system(f"alien_cp -T 16 {alien_file.replace('AnalysisResults.root', 'AO2D.root')} file:{local_dir}/AO2D.root > {os.path.join(local_dir, 'download.log')}")
    else:
        os.system(f"alien_cp -T 16 {alien_file} file:{local_dir}/{FileName} > {os.path.join(local_dir, 'download.log')}")

    with open(os.path.join(local_dir, 'download.log'), 'r') as log_file:
        log_content = log_file.read()
        if "TARGET VALID" in log_content: # already exists and is valid
            print(f"File {alien_file} ---> {local_dir}")
            print(f"    TARGET VALID (size match, source/lfn older than destination/local_file)")
        elif "No such file or directory" in log_content: # file does not exist on Alien
            print(f"File {alien_file} does not exist on Alien.")
            return None
        else: # downloaded successfully
            print(f"Downloaded {alien_file} to {local_dir}/{FileName}")
    os.remove(os.path.join(local_dir, 'download.log'))
    return f"{local_dir}/{FileName}"

def batch_download_files(download_file, LocalPath, NWorks, SortedAlienFileList, DownloadedFiles, fn):
    print(f"Starting download for file: {fn}")
    with alive_bar(sum(len(lst) for lst in SortedAlienFileList), title=f"Downloading files: {fn}") as bar:
        with ThreadPoolExecutor(max_workers=NWorks) as executor:
            tasks = (
                        (iFile, alien_file, RunNumber, LocalPath, fn)
                        for FileSingleRunList in SortedAlienFileList
                        for iFile, (alien_file, RunNumber) in enumerate(FileSingleRunList, start=1)
                    )
            futures = [executor.submit(download_file, arg) for arg in tasks]

            for f in as_completed(futures):
                try:
                    res = f.result()
                    if res:
                        DownloadedFiles.append(res)
                except Exception as e:
                    print(f"[ERROR] Task failed: {e}")
                bar()

def estimate_total_size(SortedAlienFileList):
    total_size = 0
    file_list = []
    for FileSingleRunList in SortedAlienFileList:
        for alien_file, _ in FileSingleRunList:
            file_list.append(alien_file)
    with ThreadPoolExecutor(max_workers=22) as executor:
        futures = [executor.submit(lambda f: os.popen(f"alien_ls -l {f}").read(), alien_file) for alien_file in file_list]
        for f in as_completed(futures):
            result = f.result()
            size_match = re.search(r"^[-drwxst]+\s+\w+\s+\w+\s+(\d+)", result)
            if size_match:
                total_size += int(size_match.group(1))
    print(f"Estimated total download size: {total_size / (1024**3):.2f} GB")

def download(config):
    if isinstance(config, dict):
        AlienOutputDirs = config.get('AlienOutputDirs')
        LocalPath = config.get('LocalPath')
        ForceStage = config.get('ForceStage', False)
        FileName = config.get('FileName', 'AnalysisResults.root')
    else:
        with open(args.config_download, 'r') as cfg:
            config = yaml.safe_load(cfg)
        AlienOutputDirs = config.get('AlienOutputDirs')
        LocalPath = config.get('LocalPath')
        ForceStage = config.get('ForceStage', False)
        FileName = config.get('FileName', 'AnalysisResults.root')

    NWorks = 22

    # Fetch file list needed to be downloaded from Alien
    with alive_bar(len(AlienOutputDirs), title="Fetching file list from Alien") as bar:
        with ProcessPoolExecutor(max_workers=NWorks) as executor:
            task = (
                (alien_dir, ForceStage)
                for alien_dir in AlienOutputDirs
            )
            futures = [executor.submit(get_file_list_from_alien, arg) for arg in task]
            AlienFileList = []
            for f in as_completed(futures):
                AlienFileList.extend(f.result())
                bar()
    SortedAlienFileList = [
        list(g) for _, g in groupby(sorted(AlienFileList, key=lambda x: (x[1], x[0])), key=lambda x: x[1])
    ]
    print(f"Outputs from Alien were found and sorted with Stage = {ForceStage}. Total runs: {len(SortedAlienFileList)}")
    estimate_total_size(SortedAlienFileList)

    # Download files from Alien
    DownloadedFiles = []
    if FileName == "AA":
        FileNames = ["AnalysisResults.root", "AO2D.root"]
        for fn in FileNames:
            batch_download_files(download_file, LocalPath, NWorks, SortedAlienFileList, DownloadedFiles, fn)
    else:
        batch_download_files(download_file, LocalPath, NWorks, SortedAlienFileList, DownloadedFiles, FileName)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download files from Alien")
    parser.add_argument('config_download', nargs='?', default='./config_download.yml',
                        help='Configuration file for downloading files')
    args = parser.parse_args()

    download(args)
