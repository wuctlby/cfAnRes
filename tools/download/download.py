import yaml
import argparse

parser = argparse.ArgumentParser(description="Download files from Alien")
parser.add_argument('config_download', default='config_download.yml',
                    help='Configuration file for downloading files')
args = parser.parse_args()

with open(args.config_download, 'r') as cfg:
    config = yaml.safe_load(cfg, Loader=yaml.FullLoader)

alienOutputDirs = config.get('alienOutputDirs')
localPath = config.get('localPath')
