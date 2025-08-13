'''
File Loader Module
This module provides functionality to load files from a specified path.
'''

import os
from typing import Union, List
from ROOT import TFile

class FileLoader:
    """
    FileLoader class to handle loading ROOT files from a specified path.
    """
    def __init__(self, path: Union[str, List[str]], **kwargs):
        """
        Initializes the FileLoader with a path or list of paths to ROOT files.

        Parameters:
            path (str or list of str): Path or list of paths to ROOT files.

            **kwargs: Additional keyword arguments:
                - keyFileName (str): Optional keyword to filter files by name.
        """
        keyFileName = kwargs.get('keyFileName', None)
        filetype = kwargs.get('filetype', '.root')

        self.keyFileName = keyFileName
        self.filetype = filetype
        self.raw_paths = path
        self.paths = self._get_paths()

    def _get_paths(self) -> List[str]:
        """
        Retrieves the paths of ROOT/{filetype} files from the specified path(s).
        """
        paths = []
        if not isinstance(self.raw_paths, list):
            self.raw_paths = [self.raw_paths]

        for p in self.raw_paths:
            if isinstance(p, str) and p.endswith(self.filetype):
                paths.append(p)
            elif isinstance(p, str) and os.path.isdir(p):
                for root, _, files in os.walk(p):
                    for f in files:
                        if f.endswith(self.filetype) and (self.keyFileName is None or self.keyFileName in f):
                            paths.append(os.path.join(root, f))
        return paths

    def get_files(self) -> List[TFile]:
        """
        Retrieves the opened ROOT files.
        """
        self.files = []
        for path in self.paths:
            if os.path.exists(path) and path.endswith('.root'):
                self.files.append(TFile(path, 'READ'))
            else:
                raise FileNotFoundError(f"File {path} does not exist or is not a valid ROOT file.")
