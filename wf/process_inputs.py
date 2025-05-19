import csv
import pickle
import os
from pathlib import Path
from typing import Tuple
import typing


def input_task(
        project_name: str,
    # .csv file with target/isoform info
    target_file: str,
    # Optional: directory of folders with Genbank files for each target
    genbank_dir: typing.Optional[str]=None,
    # Optional: directory of ADAPT outputs for targets
    adapt_dir: typing.Optional[str]=None,
):
    """
    Does the following:
    1. Converts input file into dictionary
    2. Maps any Genbank files to a target
    3. Maps any ADAPT files to a target
    """
    # Pass to Path object and resolve
    target_list = Path(target_file).resolve()
    target_dict = create_target_dict(target_list)

    # Make a copy of any existing ADAPT / Genbank dirs
    if genbank_dir is not None:
        for gb_folder in os.listdir(genbank_dir):
            target_dict[gb_folder]["gb_dir"] = genbank_dir + "/" + gb_folder

    if adapt_dir is not None:
        target_dict = map_adapt_dirs(target_dict, adapt_dir)

    # Pickle the target_dict
    target_pickled = "./results/"+project_name+"/targets.pkl"
    with open(target_pickled, "wb") as f:
        pickle.dump(target_dict, f)

    print("Target pickle file:", target_pickled)

# Helper function that returns a dictionary of target info
def create_target_dict(target_path):
    """
    Creates a target dictionary object that looks like the following:
    {
        "FGFR1": {
            gb_dir: "/gbdir/",
            target_groups: [
                {
                    # If we have isoforms: [], then we target all isoforms
                    isoforms: [],
                    identifier: "FGFR1_1a",
                    adapt_file: "",
                },
                {
                    isoforms: ["NC129929"],
                    identifier: "FGFR1_2c",
                    adapt_file: "",
                },
            ]
        }
    }
    """
    target_dict = {}

    # Process CSV file
    with open(target_path, encoding="utf-8-sig") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        # Append to gene_names list
        for row in csv_reader:
            group_obj = {
                "isoforms": row[2].split(",") if "," in row[2] else [],
                "identifier": row[1],
                "adapt_file": None,
            }
            if row[0] in target_dict.keys():
                target_dict[row[0]]["target_groups"].append(group_obj)
            else:
                target_dict[row[0]] = {
                    "gb_dir": None,
                    "target_groups": [group_obj]
                }

    return target_dict

def map_genbank_dirs(
    target_dict,
    genbank_dir,
):
    """
    Loops through provided genbank directory and maps any genbank folders to targets in target_dict.
    """
    for gb_folder in os.listdir(genbank_dir):
        target_dict[gb_folder]["gb_dir"] = genbank_dir + "/" + gb_folder

    return target_dict

def map_adapt_dirs(
    target_dict,
    adapt_dir,
):
    """
    Loops through provided ADAPT directory and maps any ADAPT files to targets in target_dict.
    """
    for adapt_file in os.listdir(adapt_dir):

        try:
            isoform_ids = adapt_file.split("_")[1].split(",")
            target_name = adapt_file.split("_")[0]
        except:
            target_name = adapt_file.split(".")[0]
            isoform_ids = []


        # Loop through target groups and find the one that matches
        for target_group in target_dict[target_name]["target_groups"]:
            if set(target_group["isoforms"]) == set(isoform_ids):
                target_group["adapt_file"] = adapt_dir + "/" + adapt_file
                break

    return target_dict