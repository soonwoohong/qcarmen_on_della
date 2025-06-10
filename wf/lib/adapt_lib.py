from adapt import guide_search
from adapt import alignment
from adapt import guide_search
from adapt.utils import predict_activity
# from adapt.utils import weight
from adapt.utils.version import get_project_path, get_latest_model_version
from adapt.specificity import alignment_query

from adapt import primer_search
from adapt import target_search

# General python imports
import os
import subprocess
from io import StringIO
from Bio import SeqIO
import csv
import sys
from datetime import datetime

"""
Designs Cas13 crRNAs using the sliding-window method.
"""
def sliding_window(aligned_targets, non_target_genes, gene_name, isoform_name, non_target_isoforms, output_dir):
    # If there were no FASTAs found, cop out
    if aligned_targets is None:
        return None
    # Get the path of the ADAPT directory
    dir_path = get_project_path()
    cla_path_all = os.path.join(dir_path, 'models', 'classify', 'cas13a')
    reg_path_all = os.path.join(dir_path, 'models', 'regress', 'cas13a')

    cla_version = get_latest_model_version(cla_path_all)
    reg_version = get_latest_model_version(reg_path_all)
    cla_path = os.path.join(cla_path_all, cla_version)
    reg_path = os.path.join(reg_path_all, reg_version)
    # Model defaults
    cla_thres, reg_thres = None, None

    # Set up predictor
    predictor = predict_activity.Predictor(cla_path, reg_path,
        classification_threshold=cla_thres,
        regression_threshold=reg_thres)

    # I believe we can just remove the isoform of interest here...
    # Get the sequences themselves from our alignment array
    gene_isoforms = [str(sequence.seq).upper() for sequence in aligned_targets]

    # Weights and whatnot
    # norm_weights = weight.normalize(defaultdict(lambda: 1), seq_names)
    # seq_norm_weights = [norm_weights[seq_name] for seq_name in seq_names]

    # Turn the aligned sequences into an alignment object that ADAPT can use
    aln = alignment.Alignment.from_list_of_seqs(gene_isoforms)

    # Just set other_isoforms to None right now
    other_isoforms = None

    # Specificity
    other_genes = alignment.SequenceList([str(seq.seq).upper() for seq in non_target_genes])
    if non_target_isoforms:
        other_isoforms = alignment.SequenceList([str(seq.seq).upper() for seq in non_target_isoforms])
        # alns array for ADAPT
        alns = [aln, other_isoforms, other_genes]
    else: alns = [aln, other_genes]

    # AlignmentQuerier parameters
    guide_length = 28
    # Number of mismatches allowed in the guide to be considered binding
    diff_id_mismatches = 4
    allow_gu_pairs = True
    aq = alignment_query.AlignmentQuerierWithKmerSharding(alns, guide_length, diff_id_mismatches, allow_gu_pairs)
    aq.setup()

    # Guide checking parameters... doing specificity + 1 to avoid divide by zero errors
    diff_id_frac = min(0.01, 1 / (len(non_target_genes) + 1))
    # diff_id_frac = 0.2
    alns_in_same_taxon = [0]
    guide_is_specific = aq.guide_is_specific_to_alns_fn(alns_in_same_taxon, diff_id_frac)

    def guide_is_suitable(guide):
        # Return True if guide does not hit too many sequences in alignments other than aln
        return guide_is_specific(guide)

    # Guide searcher params
    soft_guide_constraint = 1
    hard_guide_constraint = 5
    penalty_strength = 0.25
    missing_data_params = (0.5, 0.05, 1.5)

    # Set up the guide searcher
    gs = guide_search.GuideSearcherMaximizeActivity(
        aln,
        guide_length,
        soft_guide_constraint,
        hard_guide_constraint,
        penalty_strength,
        missing_data_params,
        guide_is_suitable_fn=guide_is_suitable,
        predictor=predictor,
    )

    # Get current time
    now = datetime.now()
    dt_string = now.strftime("%y_%m_%d_%H_%M")

    isoform_id = "_" + isoform_name if other_isoforms != None else ""

    print("Beginning search for " + gene_name + isoform_id + " at " + dt_string + ".")
    sys.stdout.flush()

    output_file = output_dir + gene_name + isoform_id + ".tsv"

    gs.find_guides_with_sliding_window(200,
        output_file,
        window_step=1,
        sort=False,
        print_analysis=True)


"""
Designs Cas13 crRNAs using the complete-targets method.

I've set the objective function parameters both to 0, so netiher the
amplicon length nor the number of primers matters.
"""
def complete_targets(aligned_targets, non_target_genes, gene_name, isoform_name, non_target_isoforms, output_dir):
    # If there were no FASTAs found, cop out
    if aligned_targets is None:
        return None
    # Get the path of the ADAPT directory
    dir_path = get_project_path()
    cla_path_all = os.path.join(dir_path, 'models', 'classify', 'cas13a')
    reg_path_all = os.path.join(dir_path, 'models', 'regress', 'cas13a')

    cla_version = get_latest_model_version(cla_path_all)
    reg_version = get_latest_model_version(reg_path_all)
    cla_path = os.path.join(cla_path_all, cla_version)
    reg_path = os.path.join(reg_path_all, reg_version)
    # Model defaults
    cla_thres, reg_thres = None, None

    # Set up predictor
    predictor = predict_activity.Predictor(cla_path, reg_path,
        classification_threshold=cla_thres,
        regression_threshold=reg_thres)

    # I believe we can just remove the isoform of interest here...
    # Get the sequences themselves from our alignment array
    gene_isoforms = [str(sequence.seq).upper() for sequence in aligned_targets]

    # Weights and whatnot
    # norm_weights = weight.normalize(defaultdict(lambda: 1), seq_names)
    # seq_norm_weights = [norm_weights[seq_name] for seq_name in seq_names]

    # Turn the aligned sequences into an alignment object that ADAPT can use
    aln = alignment.Alignment.from_list_of_seqs(gene_isoforms)

    # Set other_isoforms to None
    other_isoforms = None

    # Specificity
    other_genes = alignment.SequenceList([str(seq.seq).upper() for seq in non_target_genes])
    if non_target_isoforms:
        other_isoforms = alignment.SequenceList([str(seq.seq).upper() for seq in non_target_isoforms])
        # alns array for ADAPT
        alns = [aln, other_isoforms, other_genes]
    else: alns = [aln, other_genes]

    # AlignmentQuerier parameters
    guide_length = 28
    diff_id_mismatches = 4
    allow_gu_pairs = True
    aq = alignment_query.AlignmentQuerierWithKmerSharding(alns, guide_length, diff_id_mismatches, allow_gu_pairs)
    aq.setup()

    # Guide checking parameters... doing specificity + 1 to avoid divide by zero errors
    # diff_id_frac = min(0.01, 1 / (len(non_target_genes) + 1))
    diff_id_frac = 0.2
    alns_in_same_taxon = [0]
    guide_is_specific = aq.guide_is_specific_to_alns_fn(alns_in_same_taxon, diff_id_frac)

    def guide_is_suitable(guide):
        # Return True if guide does not hit too many sequences in alignments other than aln
        return guide_is_specific(guide)

    # PrimerSearcher parameters
    primer_length = 30
    primer_mismatches = 3
    primer_cover_frac = 0.1
    missing_thres = (0.5, 0.05, 1.5)
    # Set up the primer searcher, which we don't care about anyways
    ps = primer_search.PrimerSearcherMinimizePrimers(
        aln, 
        primer_length,
        primer_mismatches,
        primer_cover_frac,
        missing_thres,
        seq_groups=None,
        primer_gc_content_bounds= (0.35, 0.65))# <- got error __init__() got multiple values for argument 'primer_gc_content_bounds'


    # Guide searcher params
    soft_guide_constraint = 2
    hard_guide_constraint = 10
    penalty_strength = 0.2
    missing_data_params = (0.5, 0.05, 1.5)

    # Set up the guide searcher
    gs = guide_search.GuideSearcherMaximizeActivity(
        aln,
        guide_length,
        soft_guide_constraint,
        hard_guide_constraint,
        penalty_strength,
        missing_data_params,
        #guide_is_suitable_fn=guide_is_suitable, # it was removed from the current adapt version
        predictor=predictor,
    )

    # Target searcher parameters
    obj_type = "max"
    max_primers_at_site = None
    max_target_length = None
    obj_fn_weights = (0, 0)

    # Target searcher (combines primer search with guide search)
    ts = target_search.TargetSearcher(ps, ps, gs,
        #obj_type=obj_type, # it was removed from the current adapt version
        max_primers_at_site=max_primers_at_site,
        max_target_length=max_target_length,
        obj_weights=obj_fn_weights,
    )

    # Get current time
    now = datetime.now()
    dt_string = now.strftime("%y_%m_%d_%H_%M")

    isoform_id = "_" + isoform_name if other_isoforms != None else ""

    print("Beginning search for " + gene_name + isoform_id + " at " + dt_string + ".")
    sys.stdout.flush()

    output_file = output_dir + gene_name + isoform_id + ".tsv"

    # Actually find guides / targets and save to an output directory
    ts.find_and_write_targets(
        output_file,
    )

    print("Search for " + gene_name + isoform_id + " crRNA spacer sequence is complete.")

    # Returns the output directory when done
    return output_file

# Takes unaligned SeqIO object and returns aligned SeqIO object
def align_seqs(seqs):
    # Placeholder FASTA string
    seq_str = ''
    # Loop through sequences and add formatted FASTA to string
    for seq in seqs:
        seq_str += '>' + seq.description + '\n'
        seq_str += str(seq.seq) + '\n'
    
    # Start child subprocess for sequence alignment
    child = subprocess.Popen(['mafft', '--quiet', '-'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    # Write to STDIN
    child.stdin.write(seq_str.encode())
    # Get STDOUT and put into child_out
    child_out = child.communicate()[0].decode('utf8')
    # Convert aligned file to SeqIO object
    seq_ali = list(SeqIO.parse(StringIO(child_out), 'fasta'))
    # Close STDIN
    child.stdin.close()

    return seq_ali

"""
Gets passed a path to an ADAPT file and parses the TSV.

Should return a list of all valid spacers.

In our case, a valid spacer design is one that doesn't require multiple guides
"""
def select_spacer(adapt_path):
    valid_guides = []

    # Open TSV file
    with open(adapt_path) as file:
        # Start reading it in
        tsv_file = csv.reader(file, delimiter="\t")
        # Skip the header row
        next(tsv_file)
        # Loop through lines
        for line in tsv_file:
            # If number of guides is 1, return that entire line/guide
            if line[12] == "1": valid_guides.append(line)

    # If we loop through all of these and there's nothing, return []
    return valid_guides