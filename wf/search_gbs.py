from pathlib import Path
from typing import Tuple
import typing
import pickle
import os


from Bio import SeqIO
from datetime import datetime

from .lib.search_lib import get_gbs, save_gbs


def gene_search_task(
    project_name: str,
    ncbi_key: typing.Optional[str] = None,
    genbank_dir: typing.Optional[str] = None,
    organism: typing.Optional[str] = "Homo sapiens",
) :
    """
    Find Genbank files for each gene and return a directory of FASTAs and a 
    directory with individual Genbank files.
    """
    # Unpickle target object
    target_path = "./results/"+project_name+"/targets.pkl"
    with open(target_path, "rb") as f:
        target = pickle.load(f)

    # FASTA directory
    fasta_dir = "./results/"+project_name+"/fasta/"
    os.makedirs(fasta_dir, exist_ok=True)

    # Genbank directory
    if genbank_dir is None:
        gb_dir = "./results/"+project_name+"/gbs/"
        os.makedirs(gb_dir, exist_ok=True)
    else:
        gb_dir = Path(genbank_dir).resolve()

    # Loop through each target to check if Genbank files are present
    for target_key in target.keys():
        gb_path = target[target_key]["gb_dir"]
        if gb_path is None:
            # Get genbank files
            gbs = get_gbs(target_key, organism, ncbi_key)
            # Write these to fasta files
            SeqIO.write(gbs, fasta_dir + target_key + ".fasta", "fasta")
            # Save the genbank files
            save_gbs(target_key, gbs, gb_dir)

            # Update the target object with the directories
            target[target_key]["gb_dir"] = gb_dir + target_key + "/"
        else:
            # Otherwise, create FASTA equivalents from provided gb files
            gbs = []
            for gb_file in os.listdir(gb_dir):
                gbs.append(SeqIO.read(gb_dir + "/" + gb_file, "genbank"))
            
            # Write these to fasta files
            SeqIO.write(gbs, fasta_dir + target_key + ".fasta", "fasta")

    print("Directories:", os.listdir(gb_dir))
    print("FASTA Directories:", os.listdir(fasta_dir))

    # Rewrite the target object
    target_pickled = "./results/"+project_name+"/targets.pkl"
    with open(target_pickled, "wb") as f:
        pickle.dump(target, f)


