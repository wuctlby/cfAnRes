'''
This script is to prepare the pT differential configuration for training BDT models and do the training.
'''
import os
import time
import sys
sys.path.append('/home/wuct/ALICE/local/RTools/cfAnRes/src')
sys.path.append('/home/wuct/ALICE/local/RTools/cfAnRes/tools/ML/prepare/src')
from cfanres.file_loader import FileLoader
from copy import deepcopy
import argparse
import yaml
import subprocess

import re
def extract_pt_range(path):
	match = re.search(r'_(\d+)_(\d+)', os.path.basename(path))
	if match:
		return (int(match.group(1)), int(match.group(2)))
	return (float('inf'), float('inf'))

def _get_input_paths(inputPath, keyFileName, ptMins, ptMaxs):
	"""
	Get the sorted paths of input files from the given input path/file
	"""
	inputPaths = []
	if isinstance(inputPath, list):
		inputPath = inputPath[0]

	if inputPath.endswith('.root'):
		paths = FileLoader(os.path.dirname(inputPath), keyFileName=keyFileName).paths
	else:
		paths = FileLoader(inputPath, keyFileName=keyFileName).paths

	sorted_paths = sorted(paths, key=extract_pt_range)
	for sorted_path in sorted_paths:
		if extract_pt_range(sorted_path) in zip(ptMins, ptMaxs):
			inputPaths.append(sorted_path)

	return inputPaths

def pt_diff_training(config):
	"""
	training the BDT model for a given pT bined configuration
	"""
	cmd = (
		f'python3 /home/wuct/ALICE/local/RTools/RTools/ML/Trainning/MLClassfication_run3_derived.py {config} --train'
	)
	# print(f"Running command: {cmd}")
	# os.system(cmd)
	try:
		subprocess.run(cmd, shell=True, check=True)
	except subprocess.CalledProcessError as e:
		print(f"Error occurred while training BDT model: {e}")
	finally:
		print("BDT training process completed.")

	return cmd

if __name__ == "__main__":
	argparser = argparse.ArgumentParser(description="Prepare pT differential configuration for BDT training")
	argparser.add_argument("config_training", metavar="text",
						  default="config_training.yml", help="the training configuration file including all pT bins")
	args = argparser.parse_args()

	with open(args.config_training, 'r') as f:
		config = yaml.load(f, Loader=yaml.FullLoader)

	# I/O paths
	inputPromptPath = config['input']['prompt']
	inputFDPath = config['input']['FD']
	inputDataPath = config['input']['data']
	outputPath = config['output']['dir']

	ptMins = config.get("pt_ranges").get("min")
	ptMaxs = config.get("pt_ranges").get("max")
	nPtBins = len(ptMins)

	# handle I/O paths
	promptPath = _get_input_paths(inputPromptPath, 'McTreeForMLPrompt', ptMins, ptMaxs)
	fdPath = _get_input_paths(inputFDPath, 'McTreeForMLFD', ptMins, ptMaxs)
	dataPath = _get_input_paths(inputDataPath, 'DataTreeForML', ptMins, ptMaxs)

	# safety check
	if len(ptMins) != len(ptMaxs) or len(ptMins) != len(promptPath) or len(ptMins) != len(fdPath) or len(ptMins) != len(dataPath):
		raise ValueError(f"ptMins: {len(ptMins)}, ptMaxs: {len(ptMaxs)}, promptPath: {len(promptPath)}, fdPath: {len(fdPath)}, dataPath: {len(dataPath)}. Please check your input configuration.")

	cmds = []
	# prepare the configuration for each pT bin
	os.makedirs(f"{outputPath}/config", exist_ok=True)
	for ptMin, ptMax, prompt, fd, data in zip(ptMins, ptMaxs, promptPath, fdPath, dataPath):
		with open(f"{outputPath}/config/config_training_pT_{ptMin}_{ptMax}.yml", 'w') as f:
			config_pt = deepcopy(config)
			config_pt['input']['prompt'] = [prompt]
			config_pt['input']['FD'] = [fd]
			config_pt['input']['data'] = [data]
			config_pt['pt_ranges']['min'] = [ptMin]
			config_pt['pt_ranges']['max'] = [ptMax]
			config_pt['ml']['saved_models'] = [f"{outputPath}/pt{ptMin}_{ptMax}/ModelHandler_D0ToKPi_pT_{ptMin}_{ptMax}.pickle"]
			if ptMin >= 0 and ptMax <= 5:
				config_pt['ml']['hyper_pars_opt']['ntrials'] = 25
				config_pt['ml']['hyper_pars_opt']['njobs'] = 1
			elif ptMin > 5 and ptMax <= 10:
				config_pt['ml']['hyper_pars_opt']['ntrials'] = 25
				config_pt['ml']['hyper_pars_opt']['njobs'] = 1
				config_pt['data_prep']['filt_bkg_mass'] = 'fCpa > 0.92'
			else:
				config_pt['ml']['hyper_pars_opt']['ntrials'] = 25
				config_pt['ml']['hyper_pars_opt']['njobs'] = 2
				config_pt['data_prep']['filt_bkg_mass'] = 'fCpa > 0.92'

			yaml.dump(config_pt, f, default_flow_style=False)
		print(f"Configuration for pT bin {ptMin}-{ptMax} saved to {outputPath}/config/config_training_pT_{ptMin}_{ptMax}.yml")
		print(
			f"Training BDT model for pT bin {ptMin}-{ptMax}...\n"
			f"with prompt: {prompt}, \nfd: {fd}, \ndata: {data}"
		)
		# train the BDT model for each pT bin
		cmds.append(pt_diff_training(f"{outputPath}/config/config_training_pT_{ptMin}_{ptMax}.yml"))
	print('\n'.join(cmds))