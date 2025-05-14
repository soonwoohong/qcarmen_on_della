from datetime import datetime
from pathlib import Path
from typing import Tuple

from typing import Optional

from h5py.h5t import string_info

#from .process_inputs import input_task
#from .search_gbs import gene_search_task
#from .design_crrnas import adapt_task
#from .design_primers import primer_task
#from .add_mods import mod_task
#from .write_files import write_task

def design(
    # Genes/isoforms in .csv
    target_file: str,
    # output_directory
    output_dir: str,
    # API key for NCBI
    ncbi_key: str,
    # Organism name
    organism: str = "Homo sapiens",
    # Optional: Genbank files for targets
    genbank_dir: str = None
):
    """Description...

    qCARMEN Design Pipeline
    ----

    This workflow takes a list of genes (and optional isoforms) as input and outputs primers 
    and crRNAs for each target.
    """
    # Get start time
    now = datetime.now()
    dt_string = now.strftime("%y_%m_%d_%H_%M")

    target_obj = {"target_file": target_file,
                  "output_dir":output_dir,
                  "ncbi_key": ncbi_key,
                  "organism": organism,
                  "genbank_dir": genbank_dir,
                  "dt_string": dt_string,
    }

