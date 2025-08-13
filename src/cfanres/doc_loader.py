"""
This module provides classes and utilities for loading and accessing information
defined in the hyTask YAML configuration file.

## hytask_thn_loader (Global):
    HyTaskThnLoader: An instance of the HyTaskThnLoader class containing all tasks and their
    thnsparse informations. Initialized at module load and used globally across the project.
"""
from dataclasses import dataclass, field
from pathlib import Path
import yaml
from cfanres.utils import doc_path
from typing import Union

@dataclass
class ThnInfo:
    """
    Data class to hold thnsparse information.
    """
    key: str # the key of the thnsparse
    name: str # the name of the thnsparse
    axes: dict[str, int] = field(default_factory=dict)
    
    def get_axis(self, axisName: str) -> int:
        """
        Get the axis number by its name.
        Args:
            axisName (str): The name of the axis.
        Returns:
            int: The axis number, or -1 if not found.
        """
        return self.axes.get(axisName, -1)

    def get_axis_name(self, axisNum: int) -> str:
        """
        Get the axis name by its number.
        Args:
            axisNum (int): The number of the axis.
        Returns:
            str: The name of the axis, or None if not found.
        """
        for name, num in self.axes.items():
            if num == axisNum:
                return name
        return None
    
    # def get_axis_var(self, axis: Union[int, str]) -> str:
    #     """
    #     Get the physical variable name by axis number or name.
    #     Args:
    #         axis (Union[int, str]): The axis number or name.
    #     Returns:
    #         str: The physical variable name, or None if not found.
    #     """
    #     if isinstance(axis, int):
    #         axisName = self.get_axis_name(axis)
    #     elif isinstance(axis, str):
    #         axisName = axis

@dataclass
class HyTaskThn:
    """
    Data class to hold all thnsparse information for a task.
    """
    taskName: str
    thns: list[ThnInfo] = field(default_factory=list)
    names: list[str] = field(default_factory=list)
    keys: list[str] = field(default_factory=list)

class HyTaskThnLoader:
    """
    Class to load and manage thnsparse information from the hyTask/task_thn YAML file.
    """
    
    def __init__(self):
        self.hytask_thn_doc = Path(doc_path()) / 'hyTask' / 'task_thn.yml'
        self._task_map: dict[str, HyTaskThn] = {}
        self._load_tasks()

    def _load_tasks(self):
        """
        Load all tasks from the YAML file and store them in the task map.
        """
        if not self.hytask_thn_doc.exists():
            raise FileNotFoundError(f"HyTask thnsparse document {self.hytask_thn_doc} does not exist.")
        
        with open(self.hytask_thn_doc) as f:
            taskDoc = yaml.safe_load(f)
        
        for taskName, thnItems in taskDoc.items():
            thns, names, keys = [], [], []
            for key, info in thnItems.items():
                thns.append(ThnInfo(key=key, name=info['name'], axes=info['axes']))
                names.append(info['name'])
                keys.append(key)
            self._task_map[taskName] = HyTaskThn(taskName, thns, names, keys)

    def _get_thn_by_name(self, thnName: str) -> ThnInfo:
        """
        Get a thnsparse information by its name.
        Args:
            thnName (str): The name of the thnsparse.
        Returns:
            ThnInfo: The thnsparse information.
        """
        for task in self._task_map.values():
            for thn in task.thns:
                if thn.name == thnName:
                    return thn
        raise KeyError(f"Thn '{thnName}' not found in any task.")

    def _get_thn_by_key(self, taskName: str, thnKey: str) -> ThnInfo:
        """
        Get a thnsparse information by task name and thn key.
        Args:
            taskName (str): The name of the task.
            thnKey (str): The key of the thnsparse.
        Returns:
            ThnInfo: The thnsparse information.
        """
        task = self.get_task(taskName)
        for thn in task.thns:
            if thn.key == thnKey:
                return thn
        raise KeyError(f"Thn with key '{thnKey}' not found in task '{taskName}'.")

    def get_task(self, taskName: str) -> HyTaskThn:
        """
        Retrieve a specific task by its name.
        Args:
            taskName (str): The name of the task.
        Returns:
            HyTaskThn: An instance containing the task name and its thnsparse information.
        """
        if taskName not in self._task_map:
            raise KeyError(f"Task '{taskName}' not found in {self.hytask_thn_doc}.")
        return self._task_map[taskName]

    def get_thn(self, *args):
        """
        Retrieve a thnsparse information by its name or key.
        Args:
            *args: Either the task name and thn key, or just the thn name.
        Returns:
            ThnInfo: An instance containing the thnsparse information.
        """
        if len(args) == 2:
            taskName, thnKey = args
            task = self.get_task(taskName)
            return task.thns[task.keys.index(thnKey)]
        elif len(args) == 1:
            thnName = args[0]
            return self._get_thn_by_name(thnName)
        else:
            raise ValueError("Invalid arguments. Provide either (taskName, thnKey) or (thnName).")
hytask_thn_loader = HyTaskThnLoader()