from pathlib import Path
import pickle
import typing



from .lib.mod_lib import design_gen1, add_t7_promoter
from .lib.search_lib import read_gbs

def mod_task(
    project_name: str
):
    """
    Adds T7 promoter, 3' blockers, and required mismatches for Gen1 primers.

    Starts with blockers. Then the T7 promoters.
    """


    # Start by unpickling primer designs object
    primer_path = "./results/" + project_name + "/primer_designs.pkl"
    with open(primer_path, "rb") as f: primers = pickle.load(f)

    # Unpickle target object
    target_path = "./results/" + project_name + "/targets.pkl"
    with open(target_path, "rb") as f:
        targets = pickle.load(f)


    print("Primers:", primers)

    # Get genbank files
    gb_path = "./results/" + project_name + "/gbs"

    # Loop through each of the targets and their respective isoforms
    for gene_key in targets.keys():
        all_seqs = read_gbs(str(gb_path) + "/" + gene_key)
        for target_group in targets[gene_key]["target_groups"]:
            if primers[target_group["identifier"]] is None:
                target_group["fw_primers"] = []
                target_group["rev_primers"] = []

                # Add crRNA to target_group as well
                target_group["crRNA"] = ""
            else:
                # Get unmodified forward and reverse primers
                fw_primers = primers[target_group["identifier"]][1]
                rev_primers = primers[target_group["identifier"]][2]

                # If primers are present, we write to fw_primers and rev_primers in target_group
                target_group["fw_primers"] = get_modded(fw_primers, all_seqs, 1)
                target_group["rev_primers"] = get_modded(rev_primers, all_seqs, -1)

                # Add crRNA to target_group as well
                target_group["crRNA"] = primers[target_group["identifier"]][0]

    print("Targets:", targets)

    # Pickle the design result
    targets_pickled = "./results/" + project_name + "/targets.pkl"
    with open(targets_pickled, "wb") as f:
        pickle.dump(targets.copy(), f)



# Mods with blocker/mismatch and T7 promoter if forward primer
def get_modded(primers, all_seqs, direction=1):
    modded = []
    for primer in primers:
        modded.append(design_gen1(primer, all_seqs, direction))

    if direction == 1:
        modded = [add_t7_promoter(primer_seq) for primer_seq in modded]

    return modded