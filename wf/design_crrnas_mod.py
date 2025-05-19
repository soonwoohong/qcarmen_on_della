import os
import pickle
from Bio import SeqIO
#import multiprocessing
from joblib import Parallel, delayed


# Library imports
from .lib.adapt_lib_mod import complete_targets, sliding_window, align_seqs, select_spacer



def adapt_task(
        project_name: str
):
    """
    Starts a subprocess for all genes of interest...

    ADAPT designs are returned as a LatchDir, can be parsed through for design data.
    """
    # Start by unpickling target object
    target_path = "./results/" + project_name + "/targets.pkl"
    with open(target_path, "rb") as f:
        target = pickle.load(f)


    # Pass to Path object and resolve
    fasta_dir = "./results/" + project_name + "/fasta/"
    # Get a static list of all fasta files (aka gene names) in the directory
    all_fastas = os.listdir(fasta_dir)

    print("All FASTAs:", all_fastas)

    # ADAPT outputs
    adapt_designs = "./results/" + project_name + "/adapt_designs/"
    os.makedirs(adapt_designs, exist_ok=True)


    # Specificity directory
    spec_dir = "./results/" + project_name + "/specificity_files/"
    os.makedirs(spec_dir, exist_ok=True)


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

    for cmd in complete_cmds:
        complete_targets(*cmd)


    print("Design process complete.")



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