import os
import sys
import yaml
import argparse
import subprocess
import pathlib as pl

def get_tasks(log_path):
    tasks = {}
    suggested_tasks = None

    with open(log_path, 'r') as f:
        lines = f.readlines()
        for lnum, line in enumerate(lines):
            # prepare the tasks from 'adding' lines
            if line.startswith("Adding"):
                task = line.split("Adding")[1].strip()
            if line.startswith("However, if for local debugging"):
                suggested_tasks = lines[lnum + 1].strip()
    tasks["command"] = "-b --configuration json://configuration.json"
    for i, task in enumerate(suggested_tasks.split('|')):
        if i < len(suggested_tasks.split('|')) - 1:
            task = task.replace("-b --configuration json://configuration.json", "$OPTION").strip()
            tasks["command"] += f"{task} |\\ \n"
        else:
            # --aod-file @input_data.txt --aod-parent-access-level 1 --aod-parent-base-path-replacement "alien://;/home/alitrain/train-workdir/testdata/LFN"'
            tail_optiones = task.split("-b --configuration json://configuration.json")[-1].strip()
            task = task.replace("-b --configuration json://configuration.json", "$OPTION").replace(tail_optiones, "").strip()
            tasks["command"] += f"{task} {" ".join(tail_optiones.split(' ')[:2])} # {" ".join(tail_optiones.split(' ')[2:])}"

    return tasks["command"]


def main(config):
    with open(config, 'r') as f:
        config_data = yaml.safe_load(f)

    log_url = config_data.get('log_url').rstrip('/')
    task_path = config_data.get('task_path', '.')

    # Download the log files from the specified URL to the task_path/hyperloop directory
    os.makedirs(f"{task_path}/hyperloop", exist_ok=True)
    cmd = (
        f"curl -k --cert ~/.globus/usercert.pem --key ~/.globus/userkey.pem "
        f"--output-dir {task_path}/hyperloop "
        f"-O {log_url}/stdout.log "
        f"-O {log_url}/input_data.txt "
        f"-O {log_url}/full_config.json"
    )
    subprocess.run(cmd, shell=True)

    tasks = get_tasks(f"{task_path}/hyperloop/stdout.log")
    # Write the tasks to a file
    os.makedirs(f"{task_path}/job", exist_ok=True)
    with open(f"{task_path}/job/job.sh", 'w') as f:
        f.write(tasks)
    subprocess.run(f"cp {task_path}/hyperloop/input_data.txt {task_path}/job/", shell=True)
    subprocess.run(f"cp {task_path}/hyperloop/full_config.json {task_path}/job/", shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--config", '-c', type=str, help="configuration file",
                        default="config.yaml")
    args = parser.parse_args()

    main(
        config=args.config
    )