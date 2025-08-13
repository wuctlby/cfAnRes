from ROOT import TFile
from itertools import product
from concurrent.futures import ThreadPoolExecutor
import gc
from cfanres.doc_loader import hytask_thn_loader

def multi_proj(inputFiles, thn, axisProj, axisFilter={}): #TODO: thn could be the task name and thnkey, or the thn name
    
    if not isinstance(inputFiles, list):
        inputFiles = [inputFiles]

    initHistos = get_proj_thn(inputFiles[0], thn, axisProj, axisFilter)
    
    if len(inputFiles) <= 1:
        return initHistos
    else:
        with ThreadPoolExecutor() as executor:
            results = list(executor.map(lambda f: get_proj_thn(f, thn, axisProj, axisFilter), inputFiles[1::]))
        for result in results:
            for i, histo in enumerate(result):
                initHistos[i].Add(histo)
        return initHistos

def get_proj_thn(inputFile, thn, axisProj, axisFilters={}):
    histos = []
    file = TFile(inputFile, 'READ')
    thnName = hytask_thn_loader.get_thn(thn[0], thn[1]).name
    print(f"DEBUG: Projecting THnSparse {thnName} from file {inputFile} onto axis {axisProj} with filters {axisFilters}.")
    thnsparse = file.Get(thnName)
    thnAxes = hytask_thn_loader.get_thn(thn[0], thn[1]).axes

    if axisFilters == {}:
        histos.append(project_thn(thnsparse, axisProj, thnAxes))
    else:
        combFilters = get_filters(axisFilters)
        for iHisto, (filterName, filters) in enumerate(combFilters.items()):
            histos.append(project_thn(thnsparse.Clone(f'temp_thn_{iHisto}'), axisProj, thnAxes, filterName, filters))
        thnsparse.Delete()
        del thnsparse

    file.Close()
    gc.collect()  # Clean up memory

    return histos

def get_task_histo(inputFiles, task_histo):
    if not isinstance(inputFiles, list):
        inputFiles = [inputFiles]

    for iFile, inputFile in enumerate(inputFiles):
        file = TFile(inputFile, 'READ')
        histo = file.Get(task_histo)
        histo.SetDirectory(0)  # Detach from file
        if iFile == 0:
            histoTemp = histo.Clone(f'{histo.GetName()}s')
            histoTemp.SetDirectory(0)  # Detach from file
            histoTemp.Reset()
            histoTemp.Add(histo)
        else:
            histoTemp.Add(histo)
        file.Close()
    histoTemp.Sumw2()
    return histoTemp

def project_thn(thnToProj, axisProj, thnAxes={}, filterName='', filters={}):
    """
    Project a THnSparse (will be deleted) onto specified axes with optional filtering.
    Args:
        thn (THnSparse): The THnSparse object to project.
        axisProj (int): The axis number to project.
        thnAxes (dict): Dictionary of axes names and their corresponding numbers.
        filterName (str): Optional name for the filter.
        filters (dict, {axis1: [min, max], axis2: [min, max], ...}): Optional dictionary of filters
    Returns:
        TH1: The projected histogram with the name 'proj_{axisProj}_{filterName}'.
    """
    if not isinstance(axisProj, int) and thnAxes == {}:
        raise ValueError("axisProj must be an integer or thnAxes must be provided.")
    elif isinstance(axisProj, str):
        axisProj = thnAxes.get(axisProj, -1)
        if axisProj == -1:
            raise ValueError(f"Axis {axisProj} not found in {thnAxes}.")

    if filters != {}:
        for axisToFilter, filterRange in filters.items():
            if isinstance(axisToFilter, str):
                axisToFilter = thnAxes.get(axisToFilter, -1)
                if axisToFilter == -1:
                    raise ValueError(f"Axis {axisToFilter} not found in {thnAxes}.")
            elif not isinstance(axisToFilter, int):
                raise ValueError(f"Axis {axisToFilter} must be an integer or a valid axis name.")
            thnToProj.GetAxis(axisToFilter).SetRangeUser(filterRange[0], filterRange[1])

    print(f"DEBUG: Projecting THnSparse {thnToProj.GetName()} onto axis {axisProj} with filters {filters}.")
    histo_temp = thnToProj.Projection(axisProj)
    histo_temp.SetDirectory(0)  # Detach from file
    histo = histo_temp.Clone()
    histo.SetDirectory(0)  # Detach from file
    histo.Reset()
    histo.Add(histo_temp)
    histo.SetName(f'proj_{axisProj}_{filterName}')

    histo_temp.Delete()
    del histo_temp
    thnToProj.Delete()
    del thnToProj
    gc.collect()  # Clean up memory

    return histo

def get_filters(axisFilter):
    """
    Generate all combinations of axis filters and the suffix text for each combination.
    Args:
        axisFilter (dict, axis: [min, max] for single range or [[min, max], [min, max], ...]): Dictionary where keys are axis numbers and values are filter ranges.
    Returns:
        dict:
            str (['axis1_min1_max1_axis2_min2_max2', 'axis1_min1_max1_axis2_min2_max2', ...]): Suffix text for the filter combination.
            list ([{'axis1': [min1, max1], 'axis2': [min2, max2], ...}, {...}, ...]): List of dictionaries with filter combinations.
    """
    if not axisFilter:
        return None
    
    # list for each axis and its filter ranges
    axes = list(axisFilter.keys())
    filterRanges = list(axisFilter.values())

    # Normalize filter ranges, make sure they are all a list of lists [[min, max], ...]
    normalizedRanges = []
    for fr in filterRanges:
        if isinstance(fr, list) and len(fr) == 2 and not any(isinstance(x, list) for x in fr):
            normalizedRanges.append([fr]) # allow the case like [min, max] for single range
        elif not fr:
            print(f"Skipping empty range for axis: {axes[filterRanges.index(fr)]}")
            continue # skip empty ranges
        else:
            normalizedRanges.append(fr)

    combinations, filtersNames = [], []
    for combo in product(*normalizedRanges):
        filterDict = {}
        filtersName = ''
        for iAxis, (axis, filter_range) in enumerate(zip(axes, combo)):
            filterDict[axis] = filter_range
            filtersName += ('' if iAxis == 0 else '_') + f"axis{axis}_{filter_range[0]}_{filter_range[1]}"
        combinations.append(filterDict)
        filtersNames.append(filtersName)
        print(f"DEBUG: {filtersName} -> {filterDict}")

    return dict(zip(filtersNames, combinations))