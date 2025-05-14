from datetime import datetime

import argparse
import os


from wf.process_inputs import input_task
from wf.search_gbs import gene_search_task
from wf.design_crrnas import adapt_task
#from wf.design_primers import primer_task
#from wf.add_mods import mod_task
#from wf.write_files import write_task


def main():

    # Get start time
    now = datetime.now()
    dt_string = now.strftime("%y_%m_%d_%H_%M")

    # command-line interface
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--project_name", help="project_name", type=str, required=True)
    parser.add_argument("-t", "--target", help="input file", type = str, required=True)
    parser.add_argument("-n", "--ncbi", help="API key for NCBI", type = str, required=True)
    parser.add_argument("-org", "--organism", help="organism species", type = str, required=False)
    parser.add_argument("-nc", "--num_cpu", help="number of CPU cores", type=int, required=False, default=os.cpu_count())


    args = parser.parse_args()
    project_name = args.project_name
    target_file = args.target
    ncbi_key = args.ncbi
    organism_species = args.organism
    num_cpu = args.num_cpu

    os.makedirs("results/"+project_name, exist_ok=True)

    input_task(project_name = project_name,
               target_file = target_file)
    gene_search_task(project_name = project_name,
                     ncbi_key = ncbi_key,
                     organism=organism_species)
    adapt_task(project_name = project_name)
    adapt_task(project_name="test1")

    #primer_obj = primer_task(target_obj=target_obj, gb_dir=gb_dir, adapt_dir=adapt_dir)
    #mod_obj = mod_task(primer_obj=primer_obj, target_obj=target_obj, gb_dir=gb_dir)
    #write_task(target_obj=mod_obj, dt_string=dt_string)

if __name__ == "__main__":
    main()