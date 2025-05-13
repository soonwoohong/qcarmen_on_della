from datetime import datetime

from wf.process_inputs import input_task
from wf.search_gbs import gene_search_task
from wf.design_crrnas import adapt_task
from wf.design_primers import primer_task
from wf.add_mods import mod_task
from wf.write_files import write_task
from latch.types import LatchFile, LatchDir

# build_index(ref_genome=LatchFile("latch:///wgs/ref_genome/ecoli_rel606.fasta"))

# Get start time
now = datetime.now()
dt_string = now.strftime("%y_%m_%d_%H_%M")

# target_obj = input_task(
#     target_file = LatchFile("latch:///test_data/short_list.csv"),
# )

# gb_dir, fasta_dir = gene_search_task(target_obj=target_obj, dt_string=dt_string)

# adapt_dir = adapt_task(target_obj=target_obj, fastas=fasta_dir, specificity=False)

target_obj = LatchFile("/root/targets.pkl")
gb_dir = LatchDir("/root/gbs/")
adapt_dir = LatchDir("/root/adapt_designs/")

primer_obj = primer_task(target_obj=target_obj, gb_dir=gb_dir, adapt_dir=adapt_dir)

mod_obj = mod_task(primer_obj=primer_obj, target_obj=target_obj, gb_dir=gb_dir)

write_task(target_obj=mod_obj, dt_string=dt_string)