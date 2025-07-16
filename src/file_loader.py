'''
File Loader Module
This module provides functionality to load files from a specified path.
'''

import os
from ROOT import TFile

class FileLoader:
    """
    FileLoader class to handle loading ROOT files from a specified path.
    """
    def __init__(self, path, **kwargs):
        """
        Initializes the FileLoader with a path or list of paths to ROOT files.

        Parameters:
            path (str or list of str): Path or list of paths to ROOT files.

            **kwargs: Additional keyword arguments:
                - keyFileName (str): Optional keyword to filter files by name.
        """
        keyFileName = kwargs.get('keyFileName', None)
        self.paths = []
        self.keyFileName = keyFileName

        if isinstance(path, str):
            self.paths = [path]
        elif isinstance(path, list):
            self.paths = path
        else:
            raise ValueError("Unsupported type for path. Expected str or list of str.")
    

    def get_files(self):
        self.files = []
        for p in self.paths:
            if isinstance(p, str) and p.endswith('.root'):
                self.files.append(TFile(p, 'READ'))
            elif isinstance(p, str) and os.path.isdir(p):
                for root, _, files in os.walk(p):
                    for file in files:
                        if file.endswith('.root') and (self.keyFileName is None or self.keyFileName in file):
                            self.files.append(TFile(os.path.join(root, file), 'READ'))