"""
targets.py
Functions and classes for handling the target data.
"""

import os
import shutil
import yaml
import numpy as np
import pandas as pd

from . import __path__, targets, ligands, edges, utils


def clean_yaml_file(path):
    if os.path.exists(path):
        count = 1
        bk_path = f"{path}.bk.{count}"
        while os.path.exists(bk_path):
            count += 1
            bk_path = f"{path}.bk.{count}"
        shutil.copy(path, bk_path)
        with open(path, "r") as file:
            yaml_dict = [d for d in yaml.full_load_all(file)]
            print(yaml_dict)
        with open(path, "w") as file:
            if len(yaml_dict) == 1:
                file.write(yaml.dump(yaml_dict[0]))
            else:
                file.write(yaml.dump_all(yaml_dict))
    else:
        raise ValueError(f"Path for file {path} not found.")


def clean_metadata():
    clean_yaml_file(os.path.join(targets.data_path, "targets.yml"))
    for target in targets.target_dict:
        clean_yaml_file(
            os.path.join(targets.get_target_data_path(target), "target.yml")
        )

        clean_yaml_file(
            os.path.join(targets.get_target_data_path(target), "ligands.yml")
        )

        clean_yaml_file(os.path.join(targets.get_target_data_path(target), "edges.yml"))
