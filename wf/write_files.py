
import pickle
import typing
import csv
import os


def write_task(
    project_name: str):
    """
    Writes to CSV files.
    """
    # Start by unpickling primer designs object

    target_path = "./results/" + project_name + "/targets.pkl"
    with open(target_path, "rb") as f: targets = pickle.load(f)

    # Write primers and crRNAs to file
    # Now, output as a file
    output_dir = "./results/" + project_name + "/final_results/"
    os.makedirs(output_dir, exist_ok=True)

    with open(output_dir+"final.csv", mode='w') as output_file:
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
    with open(output_dir+"primers_only.csv", mode="w") as primer_file:
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
    
