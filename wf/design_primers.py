import csv
import os

import pickle
from joblib import Parallel, delayed


from .lib.search_lib import read_gbs
from .lib.primer_lib import design_candidates
from .lib.parallelization_lib import run_parallel_processes


def primer_task(
    project_name: str,
    # Timeout in seconds (set to 1 hour by default)
    timeout: int = 3600
):
    """
    Designs primers given a set of crRNA candidates and all sequences.

    Design processes are all parallelized using multiprocessing.
    """
    # Start by unpickling target object
    target_path = "./results/" + project_name + "/targets.pkl"
    with open(target_path, "rb") as f:
        target = pickle.load(f)

    # Flatten the target dictionary
    all_targets = []
    for target_key in target.keys():
        for target_group in target[target_key]["target_groups"]:
            all_targets.append(target_group)

    # Get genbank files
    gb_path = "./results/" + project_name + "/gbs"

    # Get guides
    adapt_path = "./results/" + project_name + "/adapt_designs"

    # Loop through each of the targets and their respective isoforms
    target_data = []
    guides = get_crrnas(adapt_path)
    # guides: {'Rpl13a': [('TTGGAGCACTGCCTTGCTCCTTCCCAGC', '7.562244534492493')], ...}
    for target in all_targets:

        all_seqs = read_gbs(str(gb_path) + "/" + target["identifier"])
        # target: {'isoforms': [], 'identifier': 'Rtp4', 'adapt_file': './results/alex5/adapt_designs/Rtp4.tsv'}
        # all_seq =[SeqRecord(seq=Seq('AGATTTGAGGAAGACAGTGAGCTGCTGAGTAACCCGGGGCCTCAGCTTCCACAC...GAA'), id='NM_023386.6', name='NM_023386', description='Mus musculus receptor transporter protein 4 (Rtp4), mRNA', dbxrefs=[])]

        #print("Target Isoforms:", target["isoforms"])
        if len(target["isoforms"]) == 0:
            target_indices = list(range(len(all_seqs)))
        else:
            target_indices = [ind for ind, seq in enumerate(all_seqs) if seq.id in target["isoforms"]]


        # Design primers + crRNA
        target_data.append((target["identifier"], [all_seqs, target_indices, guides[target["identifier"]]]))


    design_res = Parallel(n_jobs=-1, verbose=5, backend="loky")(delayed(design_function)(cmd)
                                                         for cmd in target_data)
    design_res_in_dict = {idx: res for (idx, _), res in zip(target_data, design_res)}
    #design_res = run_parallel_processes(design_function, target_data, timeout)

    print("Design process complete.", design_res_in_dict)

    # Pickle the design result
    designs_pickled = "./results/" + project_name + "/primer_designs.pkl"
    with open(designs_pickled, "wb") as f:
        pickle.dump(design_res_in_dict.copy(), f)


def design_function(params):
    try:

        return design_candidates(*params[1])
    except:
        pass


def get_crrnas(adapt_dir):
    """
    Returns crRNAs from ADAPT file.
    """
    # Get and process guides
    guides = {}
    for adapt_file in os.listdir(adapt_dir):
        isoform_guides = []
        with open(str(adapt_dir) + "/" + adapt_file) as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter="\t")
            for line in tsv_reader:
                # print(line)
                # Case: sliding window search
                if line[2] == "1": isoform_guides.append((line[-2], line[3]))
                # Case: complete target search
                if line[12] == "1": isoform_guides.append((line[-2], line[0]))

        guides[str(adapt_file).split(".")[0]] = isoform_guides

    return guides