from pathlib import Path
from typing import Tuple
import typing
import pickle
import os

from latch.types.file import LatchFile
from latch.types.directory import LatchDir
from latch.resources.tasks import small_task
from latch.ldata.path import LPath

from Bio import SeqIO
from datetime import datetime

from .lib.search_lib import get_gbs, save_gbs

@small_task
def gene_search_task(
    target_obj: LatchFile,
    output_dir: LatchDir,
    ncbi_key: typing.Optional[str] = None,
    genbank_dir: typing.Optional[LatchDir] = None,
    organism: typing.Optional[str] = "Homo sapiens",
    dt_string: typing.Optional[str] = None,
) -> Tuple[LatchDir, LatchDir, LatchFile]:
    """
    Find Genbank files for each gene and return a directory of FASTAs and a 
    directory with individual Genbank files.
    """
    # Unpickle target object
    target_path = Path(target_obj).resolve()
    with open(target_path, "rb") as f:
        target = pickle.load(f)

    # FASTA directory
    fasta_dir = "/root/fastas/"
    os.mkdir(fasta_dir)

    # Genbank directory
    if genbank_dir is None:
        gb_dir = "/root/gbs/"
        os.mkdir(gb_dir)
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
    target_pickled = "/root/targets.pkl"
    with open(target_pickled, "wb") as f:
        pickle.dump(target, f)

    outdir_path = output_dir.remote_path
    return LatchDir(gb_dir, f"{outdir_path}/{dt_string}/gbs/"), \
        LatchDir(fasta_dir, f"{outdir_path}/{dt_string}/fastas/"), \
        LatchFile(target_pickled, f"{outdir_path}/{dt_string}/tmp/targets.pkl")