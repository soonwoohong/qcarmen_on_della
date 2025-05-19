from datetime import datetime

import argparse
import os


from wf.design_crrnas import adapt_task
from wf.design_primers import primer_task
from wf.add_mods import mod_task
from wf.write_files import write_task
from wf.process_inputs import input_task


def main():

    # Get start time
    now = datetime.now()
    dt_string = now.strftime("%y_%m_%d_%H_%M")

    # command-line interface
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--project_name", help="project_name", type=str, required=True)
    parser.add_argument("-t", "--target", help="input file", type = str, required=True)
    parser.add_argument("-nc", "--num_cpu", help="number of CPU cores", type=int, required=False, default=os.cpu_count())

    args = parser.parse_args()
    project_name = args.project_name
    target_file = args.target


    num_cpu = args.num_cpu


    os.makedirs("results/"+project_name, exist_ok=True)

    input_task(project_name=project_name,
               target_file=target_file,
               genbank_dir="./results/"+project_name+"/gbs",
               adapt_dir="./results/"+project_name+"/adapt_designs")

    primer_task(project_name = project_name, timeout= 3600) # timeout in seconds

    mod_task(project_name = project_name)

    write_task(project_name = project_name)

if __name__ == "__main__":
    main()