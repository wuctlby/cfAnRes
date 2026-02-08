import os
import yaml
import argparse
from alive_progress import alive_bar
from download import download

def express_task(tasks):
    '''Express multiple download tasks based on the provided configurations.'''
    print(f"Following {len(tasks)} download tasks will be executed:")
    for i, task_config in enumerate(tasks, start=1):
        print(f"\t Task {i}:")
        print(f"\t\t to local path: {task_config.get('LocalPath')}")
        print(f"\t\t with ForceStage: {task_config.get('ForceStage', False)}")
        print(f"\t\t downloading file: {task_config.get('FileName', 'AnalysisResults.root')}")
    input("Press Enter to continue...")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download multiple files from Alien directories.")
    parser.add_argument('--config', type=str, default="./config_multi_download.yaml",
                        help='Path to the configuration YAML file.')
    args = parser.parse_args()

    with open(args.config, 'r') as cfg:
        config = yaml.safe_load(cfg)

    tasks = [config[key] for key in config if key.startswith('task')]
    express_task(tasks)
    for i, task_config in enumerate(tasks, start=1):
        print(f"Starting download task {i}...")
        download(task_config)
        print(f"Completed download task {i}.")