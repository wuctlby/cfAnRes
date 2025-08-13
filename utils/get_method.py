import os
from typing import Union, List

def get_paths(rawPaths: Union[str, List[str]], fileType: str = '.root', keyFileName: str = None) -> List[str]:
    """
    Retrieve paths of files from the specified raw paths with the given filetype and optional keyFileName.
    """
    paths = []
    if not isinstance(rawPaths, list):
        rawPaths = [rawPaths]

    for p in rawPaths:
        # Check if the path is a string and ends with the specified filetype and matches the keyFileName if provided
        if isinstance(p, str) and p.endswith(fileType) and (keyFileName is None or all(k in p for k in keyFileName.split(';'))):
            paths.append(p)

        # If the path is a directory, walk through it to find files
        elif isinstance(p, str) and os.path.exists(p):
            for root, _, files in os.walk(p):
                for f in files:
                    if f.endswith(fileType) and (keyFileName is None or all(k in f for k in keyFileName.split(';'))):
                        paths.append(os.path.join(root, f))
    return paths