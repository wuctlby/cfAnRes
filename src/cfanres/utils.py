import cfanres
from pathlib import Path

def doc_path():
    """
    Returns the path to the doc directory.
    """
    print(cfanres.__file__)
    cfAnResRoot = Path(cfanres.__file__).parent
    
    devDocPath = cfAnResRoot.parent.parent / 'doc'
    print(f"Checking for documentation in: {devDocPath}")
    if devDocPath.exists():
        return devDocPath
    
    installDocPath = cfAnResRoot / 'data'
    if installDocPath.exists():
        return installDocPath

    raise FileNotFoundError("doc directory not found in cfAnRes package or parent directory.")