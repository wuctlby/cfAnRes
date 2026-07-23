import os
import re
import yaml
import argparse
import subprocess
from alive_progress import alive_bar
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor
from itertools import groupby

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def run_command(cmd):
    """Execute a shell command and return its stdout as a string."""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout.strip()

def filter_staged_files(staged_files):
    """Filter staged files to keep only those from Stage_(max-1)."""
    stage_lines = []
    no_stage_lines = []
    stage_pattern = re.compile(r"Stage_(\d+)")

    for line in staged_files:
        match = stage_pattern.search(line)
        if match:
            stage_lines.append((int(match.group(1)), line))
        else:
            no_stage_lines.append(line)

    if not stage_lines: # Return all if no staged files found
        return no_stage_lines

    max_stage = max(stage_lines, key=lambda x: x[0])[0]
    if max_stage == 1: # Return all if only Stage_1 exists
        return no_stage_lines
    
    target_stage = max_stage - 1
    return [line for num, line in stage_lines if num == target_stage]

def get_file_list_from_alien(args):
    """Get the list of files from Alien directory."""
    alien_dirs, force_stage, run_num, file_name = args
    alien_dirs = [alien_dirs] if not isinstance(alien_dirs, list) else alien_dirs

    file_list = []
    for alien_dir in alien_dirs:
        # Determine HyJobID safely
        parts = alien_dir.strip().strip('/').split('/')
        hy_job_id = parts[-2] if parts[-1] in ["AOD", "001"] else parts[-1]
        
        # Determine subdirectory name
        dir_id = str(run_num) if run_num is not None else hy_job_id
        temp_file_list = []

        if not force_stage:
            output = run_command(f"alien_find {alien_dir} {file_name} -r")
            lines = [line.strip() for line in output.split('\n') if line.strip()]
            temp_file_list = [(line, dir_id) for line in lines]
            
            if f"{alien_dir}/{file_name}" in [line[0] for line in temp_file_list]:
                temp_file_list = [(f"{alien_dir}/{file_name}", dir_id)]
            elif len(temp_file_list) > 0:
                print(f"Warning: Found {len(temp_file_list)} unmerged files under {alien_dir}, skipping since Stage is not requested.")
                temp_file_list = []
        
        # Check unmerged or force stage condition
        if force_stage and not temp_file_list:
            output = run_command(f"alien_find {alien_dir} */{file_name}")
            staged_files = list(set(line.strip() for line in output.split('\n') if line.strip()))
            
            filtered_staged = filter_staged_files(staged_files)
            temp_file_list.extend([(line, dir_id) for line in filtered_staged])
            
            if f"{alien_dir}/{file_name}" in staged_files:
                temp_file_list = [(f"{alien_dir}/{file_name}", dir_id)]

        file_list.extend(temp_file_list)
        
    return file_list

def download_file(args):
    """Download a file from Alien to local path."""
    iFile, alien_file, sub_dir, local_path, file_name = args
    if not alien_file:
        return None
        
    local_dir = os.path.join(local_path, str(sub_dir), f"{iFile:04d}")
    os.makedirs(local_dir, exist_ok=True)
    
    source_file = alien_file.replace('AnalysisResults.root', 'AO2D.root') if file_name == "AO2D.root" else alien_file
    target_file = os.path.join(local_dir, file_name)
    
    # Execute alien_cp and capture both stdout and stderr
    cmd = f"alien_cp -T 16 {source_file} file:{target_file}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    log_content = result.stdout + result.stderr

    if "TARGET VALID" in log_content:
        print(f"File {source_file} ---> {local_dir}\n    TARGET VALID (size match, source older than destination)")
    elif "No such file or directory" in log_content:
        print(f"File {source_file} does not exist on Alien.")
        return None
    else:
        print(f"Downloaded {source_file} to {target_file}")
        
    return target_file

def batch_download_files(local_path, n_works, sorted_alien_file_list, downloaded_files, fn, total_size=None):
    """Batch download files using ThreadPoolExecutor."""
    print(f"Starting download for file: {fn}")
    total_files = sum(len(lst) for lst in sorted_alien_file_list)
    size_str = f"({total_size / (1024**3):.4f} GB)" if total_size else ""
    
    with alive_bar(total_files, title=f"Downloading files {size_str}: {fn}") as bar:
        with ThreadPoolExecutor(max_workers=n_works) as executor:
            tasks = (
                (iFile, alien_file, sub_dir, local_path, fn)
                for run_list in sorted_alien_file_list
                for iFile, (alien_file, sub_dir) in enumerate(run_list, start=1)
            )
            futures = [executor.submit(download_file, arg) for arg in tasks]

            for f in as_completed(futures):
                try:
                    res = f.result()
                    if res:
                        downloaded_files.append(res)
                except Exception as e:
                    print(f"[ERROR] Task failed: {e}")
                bar()

def estimate_total_size(sorted_alien_file_list):
    """Estimate total size by parsing alien_ls output."""
    total_size = 0
    file_list = [alien_file.replace(',', '') for run_list in sorted_alien_file_list for alien_file, _ in run_list]
    chunk_size = 20
    
    with ThreadPoolExecutor(max_workers=chunk_size) as executor:
        futures = []
        for i in range(0, len(file_list), chunk_size):
            chunk = " ".join(file_list[i:i+chunk_size])
            futures.append(executor.submit(run_command, f"alien_ls -l {chunk}"))
            
        for f in as_completed(futures):
            output = f.result()
            # Parse sizes dynamically handling multiple files per chunk
            for line in output.split('\n'):
                parts = line.split()
                # alien_ls format: [perms] [user] [group] [size] [month] [day] [time/year] [filename]
                # Therefore, size is at index 3
                if len(parts) >= 4 and parts[0].startswith(('-', 'd', 'l')): 
                    try:
                        total_size += int(parts[3]) 
                    except ValueError:
                        pass
                        
    print(f"Estimated total download size: {total_size / (1024**3):.4f} GB")
    return total_size

def download(config_source):
    """Main function to parse config and initiate download."""
    if isinstance(config_source, dict):
        config = config_source
    else:
        with open(config_source, 'r') as cfg:
            config = yaml.safe_load(cfg)
            
    alien_output_dirs = config.get('AlienOutputDirs')
    local_path = config.get('LocalPath')
    force_stage = config.get('ForceStage', False)
    file_name = config.get('FileName', 'AnalysisResults.root')
    run_num_match = config.get('RunNumMatch', False)
    run_list = config.get('RunList', [])
    sample = config.get('sample', 1)  # Default to 1 if not specified
    allowance = config.get('allowance', -1)  # Default to -1 if not specified

    n_cpu = os.cpu_count() or 16
    n_works = min(int(n_cpu / 3), 16)

    if run_num_match and len(run_list) != len(alien_output_dirs):
        print("Error: RunList length does not match AlienOutputDirs length.")
        return

    # Determine exact file name to fetch
    file_name_to_fetch = "AO2D.root" if file_name == "AO2D.root" else ("AnalysisResults.root" if file_name == "AA" else file_name)

    # Fetch file list
    with alive_bar(len(alien_output_dirs), title="Fetching file list from Alien") as bar:
        with ProcessPoolExecutor(max_workers=n_works) as executor:
            if run_num_match:
                tasks = ((alien_dir, force_stage, run_num, file_name_to_fetch) for alien_dir, run_num in zip(alien_output_dirs, run_list))
            else:
                tasks = ((alien_dir, force_stage, None, file_name_to_fetch) for alien_dir in alien_output_dirs)
            
            futures = [executor.submit(get_file_list_from_alien, arg) for arg in tasks]
            alien_file_list = []
            for f in as_completed(futures):
                result = f.result()
                if result:
                    alien_file_list.extend(result)
                bar()

    # Sort upfront by group key (x[1]) and
    alien_file_list.sort(key=lambda x: (x[1], x[0]))

    sorted_alien_file_list = []

    # Calculate unique groups efficiently using set
    n_group = len(set(x[1] for x in alien_file_list))

    # Calculate base allowance per group
    group_allowance = 0
    if allowance > 0 and n_group > 0:
        group_allowance = allowance // n_group

    for _, group in groupby(alien_file_list, key=lambda x: x[1]):
        group_list = list(group)
        total_files = len(group_list)
        
        if allowance > 0:
            # Use a local variable to avoid overwriting global allowance
            current_allowance = min(group_allowance, total_files)
            
            if current_allowance == 0:
                continue
                
            # Distribute files evenly across the group
            dynamic_step = max(1, total_files // current_allowance) if sample > 0 else 1
            sampled_group = group_list[::dynamic_step][:current_allowance]
        else:
            # Use fixed step if no allowance is set
            fixed_step = sample if sample > 0 else 1
            sampled_group = group_list[::fixed_step]
            
        # Append as a sub-list directly to form the final grouped structure
        if sampled_group:
            sorted_alien_file_list.append(sampled_group)

    print(f"Outputs from Alien found (Stage = {force_stage}). Total runs/jobs: {len(sorted_alien_file_list)}")
    total_size = estimate_total_size(sorted_alien_file_list)

    # Download files
    downloaded_files = []
    files_to_download = ["AnalysisResults.root", "AO2D.root"] if file_name == "AA" else [file_name]
    
    for fn in files_to_download:
        batch_download_files(local_path, n_works, sorted_alien_file_list, downloaded_files, fn, total_size)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download files from Alien")
    parser.add_argument('config_download', nargs='?', default='./config_download.yml',
                        help='Configuration file for downloading files')
    args = parser.parse_args()

    download(args.config_download)