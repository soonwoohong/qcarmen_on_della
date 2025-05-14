import csv
import os
from pathlib import Path
import pickle
import typing


from typing import List

from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing
from datetime import datetime

# Library imports
from .lib.adapt_lib import complete_targets, sliding_window, align_seqs, select_spacer

def allocate_cpu(
    fastas: LatchDir,
    **kwargs
) -> int: # number of cores to allocate
    fasta_dir = Path(fastas).resolve()
    # Get a static list of all fasta files (aka gene names) in the directory
    all_fastas = os.listdir(fasta_dir)
    return min(len(all_fastas), 32)

# @custom_task(cpu=allocate_cpu, memory=128)
@custom_task(cpu=24, memory=128)
def adapt_task(
    target_obj: LatchFile, 
    fastas: LatchDir,
    output_dir: LatchDir,
    adapt_dir: typing.Optional[LatchDir] = None,
    specificity: bool = False,
    dt_string: typing.Optional[str] = None,
) -> LatchDir:
    """
    Starts a subprocess for all genes of interest...

    ADAPT designs are returned as a LatchDir, can be parsed through for design data.
    """
    # Start by unpickling target object
    target_path = Path(target_obj).resolve()
    with open(target_path, "rb") as f: target = pickle.load(f)

    # Pass to Path object and resolve
    fasta_dir = Path(fastas).resolve()
    # Get a static list of all fasta files (aka gene names) in the directory
    all_fastas = os.listdir(fasta_dir)

    print("All FASTAs:", all_fastas)

    # ADAPT outputs
    if adapt_dir is None:
        adapt_designs = "/root/adapt_designs/"
        os.mkdir(adapt_designs)
    else:
        adapt_designs = Path(adapt_dir).resolve()

    # Specificity directory
    spec_dir = "/root/specificity_files/"
    os.mkdir(spec_dir)

    # Go ahead and add all FASTAs to a dictionary
    fasta_dict = {}
    for fasta_file in all_fastas:
        # Parse into a list with file name as a key
        fasta_dict[fasta_file] = list(SeqIO.parse(str(fasta_dir) + "/" + fasta_file, "fasta"))

    print(fasta_dict)

    # ADAPT commands for each gene of interest
    complete_cmds = []
    sliding_cmds = []
    # Generate commands for each gene
    for target_key in target.keys():
        for target_group in target[target_key]["target_groups"]:
            adapt_cmd = generate_adapt_cmd(target_key, target_group, fasta_dict, adapt_designs)
            if len(target_group["isoforms"]) > 0: sliding_cmds.append(adapt_cmd)
            else: complete_cmds.append(adapt_cmd)

    # Multiprocessing
    num_cores = os.cpu_count()
    print("CPU Count:", num_cores)
    with multiprocessing.Pool(num_cores) as pool:
        # print(pool.starmap(sliding_window, sliding_cmds))
        print(pool.starmap(complete_targets, complete_cmds))

    print("Design process complete.")

    outdir_path = output_dir.remote_path
    print("Output directory:", outdir_path)
    return LatchDir(adapt_designs, f"{outdir_path}/{dt_string}/adapt/")

def generate_adapt_cmd(
    target_key: str,
    target_group: dict,
    fasta_dict: dict,
    adapt_designs: str,
):
    """
    Creates ADAPT command based on target group info.

    We have two cases: design for specific isoforms (sliding window) or design for all isoforms (complete targets).
    1. Specific Isoforms: select only the isoforms of interest and run sliding window, no specificity requirement
    2. All Isoforms: run complete targets on every isoform present and design against all other full genes
    """
    # If there are no FASTAs...
    if len(fasta_dict[target_key + ".fasta"]) == 0:
        # Add None to adapt_cmds: this will cause complete_targets to error out before running
        return (None, None, target_key, "", None, adapt_designs)
    
    # If isoforms are specified...
    if len(target_group["isoforms"]) > 0:
        # First grab all FASTAs
        all_isoforms = fasta_dict[target_key + ".fasta"]
        # Now remove the isoform from the list of all isoforms
        other_isoforms = [fasta for fasta in all_isoforms if fasta.description not in target_group["isoforms"]]
        # Grab the isoform of interest and put it in a list
        isoform_fasta = [fasta for fasta in all_isoforms if fasta.description in target_group["isoforms"]]
    else:
        other_isoforms = None
        isoform_fasta = fasta_dict[target_key + ".fasta"]

    # Start by obtaining an alignment for the target transcripts
    aligned_targets = align_seqs(isoform_fasta)
    other_genes = [indiv_fasta for fasta_name in fasta_dict.keys() if fasta_name != target_key + ".fasta" for indiv_fasta in fasta_dict[fasta_name]]

    return (aligned_targets, other_genes, target_key, target_group["identifier"], other_isoforms, adapt_designs)