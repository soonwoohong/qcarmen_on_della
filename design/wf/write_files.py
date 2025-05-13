from pathlib import Path
import pickle
import typing
import csv

from latch.types.file import LatchFile
from latch.types.directory import LatchDir
from latch.resources.tasks import small_task

@small_task
def write_task(
    target_obj: LatchFile,
    output_dir: LatchDir,
    dt_string: typing.Optional[str] = None,
) -> typing.Tuple[LatchFile, LatchFile]:
    """
    Writes to CSV files.
    """
    # Start by unpickling primer designs object
    target_path = Path(target_obj).resolve()
    with open(target_path, "rb") as f: targets = pickle.load(f)

    # Write primers and crRNAs to file
    # Now, output as a file
    with open("/root/final.csv", mode='w') as output_file:
        output_writer = csv.writer(output_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        
        # Header row
        output_writer.writerow([
            "Gene",
            "Target ID", 
            "Cas13 Spacer",
            # "Guide Score",
        ])

        for gene_key in targets.keys():
            for target_group in targets[gene_key]["target_groups"]:
                output_writer.writerow([
                    gene_key,
                    target_group["identifier"],
                    target_group["crRNA"],
                    # target_group["guide_score"],
                ])

    output_file.close()

    # Now, we want just the primer output
    with open("/root/primers_only.csv", mode="w") as primer_file:
        output_writer = csv.writer(primer_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        # Header row
        output_writer.writerow([
            "Gene",
            "Target ID",
            "Direction", 
            "Sequence", 
        ])

        for gene_key in targets.keys():
            for target_group in targets[gene_key]["target_groups"]:
                num_fw = len(target_group["fw_primers"])
                num_rev = len(target_group["rev_primers"])
                for ind, primer in enumerate(target_group["fw_primers"]):
                    output_writer.writerow([
                        gene_key,
                        target_group["identifier"] if num_fw == 1 else target_group["identifier"] + "_" + str(ind + 1),
                        "Forward",
                        primer,
                    ])
                for ind, primer in enumerate(target_group["rev_primers"]):
                    output_writer.writerow([
                        gene_key,
                        target_group["identifier"] if num_rev == 1 else target_group["identifier"] + "_" + str(ind + 1),
                        "Reverse",
                        primer,
                    ])

    primer_file.close()
    
    # return LatchFile("/root/final.csv"), LatchFile("/root/primers_only.csv")
    outdir_path = output_dir.remote_path
    return LatchFile("/root/final.csv", f"{outdir_path}/{dt_string}/results/crRNAs.csv"),\
        LatchFile("/root/primers_only.csv", f"{outdir_path}/{dt_string}/results/primers_only.csv")