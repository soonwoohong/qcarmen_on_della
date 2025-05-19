import primer3
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import numpy as np
import itertools
import math
from typing import List

def design_candidates(
    # GenBank records of all isoforms
    all_seqs,
    # Indices of targets within all_seqs
    target_indices,
    # List of tuples with (guide_sequences, predicted_activity), these should all be unique
    guides,
    desired_amplicon_length=300,
):
    """
    Finds primer pool + crRNA candidates. End-to-end search function. Runs candidates through optimizers.

    Looking for average amplicon length, primer design, crRNA activity.

    Final consideration that I need to implement: search_primer_sets will accurately return
    crRNA + exon 1 + exon 2 tuples that will SPECIFICALLY amplify our targets. Here's the issue, exons
    can have similarities between them. So just because the exon set is "specific" for our targets,
    the ultimate primer set + crRNA that we choose within those junctions may not necessarily be so.

    So, we need to do a final check of the primers we end up choosing for specificity.
    """
    # Get target GenBank records
    target_gbs = [seq for ind, seq in enumerate(all_seqs) if ind in target_indices]


    # Create exon map
    exon_map = create_exon_map(all_seqs)

    # Get primer matrix
    primer_mat, primer_locs = get_primer_loc_matrix(target_gbs, all_seqs, exon_map)

    # Get shared points across targets
    shared_points = get_shared_points(target_gbs, all_seqs, exon_map)

    print("Shared points:", shared_points, primer_mat, primer_locs, exon_map, target_gbs)

    # Search for primers, this returns indices of primer_locs
    # Convert target_indices to a list of 0s and 1s
    primer_search_res = search_primer_sets(
        primer_mat, 
        primer_locs, 
        [1 if ind in target_indices else 0 for ind in range(len(all_seqs))], 
        shared_points
    )


    # What exons or exon-junctions are the guides in?
    # Create a dictionary with exon/exon junction as keys and (guide, predicted_activity) as values
    guide_dict = {}
    for guide in guides:
        # Get exon/exon junction
        exon = get_exon_from_seq(guide[0], all_seqs[target_indices[0]], exon_map)
        if exon in guide_dict:
            guide_dict[exon].append(guide)
        else:
            guide_dict[exon] = [guide]

    print("Guides:", guides)
    print("Guide Dict:", guide_dict)

    # Prune primer sets that don't have a corresponding crRNA in that exon, some function here
    valid_crRNAs = []
    # Loop through results
    for crRNA_loc in primer_search_res:
        if crRNA_loc[0] in guide_dict: valid_crRNAs.append(crRNA_loc)

    print("Valid crRNAs:", valid_crRNAs)

    # Flatten the valid_crRNAs list so that we have groupings of (crRNA, forward_primers, reverse_primers)
    combos_flat = []
    # Simultaneously get exons/junctions that we should be covering with primers (note that these are also indices of primer_locs)
    covered_points = []
    for crRNA_loc in valid_crRNAs:
        # Might have multiple sets of primer combos e.g. 1 left + 3 right, 2 left + 2 right, 3 left + 1 right
        lr_sets = crRNA_loc[1]

        for lr_set in lr_sets:
            # Loop through each primer location in the combined set of forward and reverse primers
            for item in lr_set[0] + lr_set[1]:
                if item in covered_points: continue
                covered_points.append(item)
            # Add each left-right set to combos_flat with its corresponding crRNA location as a tuple
            combos_flat += [(crRNA_loc[0], lr_set)]

    # Try to design primers for each exon or exon junction and return average amplicon length, GC content, Tm
    primer_designs = {}
    # Design for all the points in covered_points
    for point_ind in covered_points:
        # Get the point it's referring to
        point = primer_locs[point_ind]

        # If the point is a tuple, design primer across the junction, provide the exons to the left and right
        if isinstance(point, tuple):
            first_exon = exon_map[point[0]]
            second_exon = exon_map[point[1]]
            # Combine the first and second exons for search
            combined_seq = Seq(str(first_exon) + str(second_exon))
            primer_designs[point] = find_primers(combined_seq, target_gbs, required_coverage=len(first_exon))
        else:
            primer_designs[point] = find_primers(exon_map[point], target_gbs)

    # We also want to store in memory primer sets that we've already tested
    # Keys are tuples of exon / junction pairs, values are objects
    # Objects hold optimal primer pair for that point pair and some score that assesses
    # things like GC content, Tm, etc.
    tested_point_pairs = {}

    print("combos flat:", combos_flat)

    # Store results of all combos
    combo_res = []
    # A combo contains a crRNA location, and one possible pool of n forward primers, and one possible pool of m reverse primers
    #combos flat: [(1, ([0], [1])), (1, ([1], [1])), (1, ([2], [1]))]

    for combo in combos_flat:
        # First loop through all of the primers and see if any of them don't have any candidate primers
        missing_primers = False
        # Because find_primers may not always return primers e.g. if Tm is invalid for everything in the window
        for primer in set(combo[1][0] + combo[1][1]):
            if len(primer_designs[primer_locs[primer]]) == 0:
                missing_primers = True
                break

        # Skip testing this combo if one of the primer locations doesn't have any primer candidates
        if missing_primers: continue

        # Otherwise, create all combinations of forward and reverse primers from the combo
        point_combos = list(itertools.product(combo[1][0], combo[1][1]))

        # Loop through point combos (combination of two primer locations from the set of valid primer locations)
        for point_combo in point_combos:            
            # Check first to see if we've tested it already
            if point_combo in tested_point_pairs: continue

            first_primers = [design[0] for design in primer_designs[primer_locs[point_combo[0]]]]
            second_primers = [str(Seq(design[0]).reverse_complement()) for design in primer_designs[primer_locs[point_combo[1]]]]

            if not (len(first_primers) > 0 and len(second_primers) > 0): continue

            # Get the primer pair with the most optimal primer length
            optimal_primer_pair = optimize_primer_length(
                first_primers,
                second_primers,
                target_gbs
            )

            # This is a tuple of the primer pair and the amplicon length
            tested_point_pairs[point_combo] = optimal_primer_pair

        # Dictionary storing scores for individual primers for this particular combo
        fw_primer_dict = {}
        rev_primer_dict = {}

        # Loop through combos again and calculate composite score
        for point_combo in point_combos:
            """
            To explain what's going on here:

            We're looking at possible combinations of a forward pool and reverse pool. The forward and reverse
            pools are primer locations. In a previous loop, we looked at different combinations of primer locations
            and identified an optimal primer pair based on amplicon length for each combination of forward LOCATION
            and reverse LOCATION.

            Now what we're doing is we're looking at those tested point pairs and we're giving each specific PRIMER 
            SEQUENCE a "point/award" if it appears in a combination of primer locations. 
            """
            if point_combo not in tested_point_pairs: continue

            # Calculate a score for the primer pair
            pair_score = int(np.sum(primer_mat[:, point_combo[0]] * primer_mat[:, point_combo[1]])) * (-abs(desired_amplicon_length - tested_point_pairs[point_combo][1]))

            # Initialize if needed
            if point_combo[0] not in fw_primer_dict: fw_primer_dict[point_combo[0]] = {}
            if point_combo[1] not in rev_primer_dict: rev_primer_dict[point_combo[1]] = {}

            primer_1 = tested_point_pairs[point_combo][0][0]
            primer_2 = tested_point_pairs[point_combo][0][1]

            # Add scores for primers to dictionary
            if primer_1 not in fw_primer_dict[point_combo[0]]: fw_primer_dict[point_combo[0]][primer_1] = 0
            if primer_2 not in rev_primer_dict[point_combo[1]]: rev_primer_dict[point_combo[1]][primer_2] = 0
                                                                                                
            # Add pair score 
            fw_primer_dict[point_combo[0]][primer_1] += pair_score
            rev_primer_dict[point_combo[1]][primer_2] += pair_score


        """
        It's then in these subsequent loops that we loop through each primer in the forward and reverse pools
        and choose the PRIMER SEQUENCE that had the highest score. A high score means that that primer sequence
        appeared as an "optimal" sequence in the most combinations of forward + reverse primers.
        """
        # For each forward position and reverse position, get the primer with the highest score
        forward_final = []
        forward_score = 0
        for forward_pos in combo[1][0]:
            # Sort the primers at the location by score and add to forward list
            sorted_best = sorted(fw_primer_dict[forward_pos].items(), key=lambda x: x[1], reverse=True)[0]
            forward_final.append(
                sorted_best[0]
            )
            forward_score += sorted_best[1]

        # Same for reverse primers
        reverse_final = []
        reverse_score = 0
        for reverse_pos in combo[1][1]:
            sorted_best = sorted(rev_primer_dict[reverse_pos].items(), key=lambda x: x[1], reverse=True)[0]
            reverse_final.append(
                sorted_best[0]
            )
            reverse_score += sorted_best[1]

        """
        We need to then store this composite score in a dictionary of "flat" combos (crRNA, forward_pool, reverse_pool)
        and then things like composite score, average amplicon length, etc.

        Honestly, I think just a single composite score should be fine, then we can rank.
        """
        combo_res.append((combo, forward_final, reverse_final, forward_score + reverse_score))

    # Sort primer sets by score
    combo_res.sort(key=lambda x: x[3], reverse=True)

    non_target_seqs = [seq for ind, seq in enumerate(all_seqs) if ind not in target_indices]
    print('non_target_seqs', non_target_seqs)
    print("Combo Res:", combo_res)

    # Choose the primer set with the highest score, then find the crRNA with the highest score and use that
    for res in combo_res:
        # Primer pair
        fw_primers = res[1]
        rev_primers = res[2]

        # Check to see if we have a crRNA that exists within the result, search dictionary of crRNAs
        if res[0][0] in guide_dict:
            # If we do, sort crRNAs by score
            current_crRNAs = guide_dict[res[0][0]]
            current_crRNAs.sort(key=lambda x: x[1], reverse=True)

            # Loop through crRNAs and check specificity
            for crRNA in current_crRNAs:
                crRNA_sequence = crRNA[0]
                specificity_check = design_specificity(crRNA_sequence, fw_primers, rev_primers, non_target_seqs)

                print("Check params:", crRNA_sequence, fw_primers, rev_primers, non_target_seqs, target_gbs)

                # Check that the crRNA is in between the forward and reverse primers
                crrna_between_primers = check_between_primers(crRNA_sequence, fw_primers, rev_primers, target_gbs)

                if specificity_check and crrna_between_primers:
                    return (crRNA_sequence, fw_primers, rev_primers)

    # If we go through all of the results and none of them have a valid crRNA, just return None
    return None

def check_between_primers(
    crRNA: str, 
    fw_pool: List[str], 
    rev_pool: List[str], 
    target_gbs: list # List of SeqRecord objects
) -> bool:
    """
    Checks that the crRNA is in between the forward and reverse primers.
    """
    fw_rev_combos = list(itertools.product(fw_pool, rev_pool))
    between_all = True
    for gb in target_gbs:
        cr_ind = str(gb.seq).upper().find(crRNA.upper())
        if cr_ind == -1:
            between_all = False
            break

        for fw, rev in fw_rev_combos:
            fw_ind = str(gb.seq).upper().find(fw.upper())
            rev_ind = str(gb.seq).upper().find(str(Seq(rev).reverse_complement()).upper())

            if fw_ind == -1 or rev_ind == -1:
                between_all = False
                break

            if cr_ind < fw_ind or cr_ind > rev_ind:
                between_all = False
                break

    return between_all

def find_primers(
    # Search region provided 5' to 3'
    search_seq,
    # Target sequence provided 5' to 3'
    target_seqs,
    # Window parameters
    window_size=200,
    step_size=100,
    melting_temp=61,
    # Maximum amount Tm can differ from desired melting temp
    max_tm_diff=1.5,
    # GC content requirements
    min_gc=40,
    max_gc=60,
    # None or index in search_seq
    required_coverage=None,
    # How much the primer needs to cover the flanking sides by
    min_flank=5,
    # Minimum primer length
    min_primer_length=18,
    # Maximum possible primer length
    max_primer_length=30,
):
    """
    Given a stretch of sequence, finds all primers that fit the parameters.

    If a particular point needs to be covered, then function ensures primer covers that point. 
    """
    primer_res = []

    # Set the window start to 0
    current_window_start = 0
    # Loop through every possible window by step size and ensure that we have coverage
    while current_window_start < len(search_seq):
        # Set search start, 0 by default unless required_coverage
        search_start = current_window_start
        search_end = current_window_start + window_size - min_primer_length

        # If there's a junction we need to cover, adjust accordingly
        if required_coverage and (search_end < required_coverage or search_start > required_coverage):
            current_window_start += step_size
            continue
        elif required_coverage:
            search_start = required_coverage - max_primer_length
            search_end = required_coverage

        valid_primers = []
        # Loop through search range
        for start in range(search_start, search_end):
            # Check if we cover the junction if needed
            if required_coverage:
                if not (
                    start <= required_coverage - min_flank and 
                    start + max_primer_length >= required_coverage + min_flank
                ): continue

            # Check to make sure we aren't out of the bounds of the sequence
            if start + max_primer_length > len(search_seq): continue

            # Loop through different primer lengths until we get some valid primers
            for primer_len in range(min_primer_length, max_primer_length):
                # Get sequence, Tm, and GC content
                primer_seq = search_seq[start:start+primer_len]
                primer_tm = primer3.calcTm(str(primer_seq))
                gc_content = GC(primer_seq)

                # If the primer melting temperature is within bounds, add to valid_primers
                if abs(primer_tm - melting_temp) < max_tm_diff and min_gc <= gc_content <= max_gc:
                    valid_primers.append((str(primer_seq), primer_tm, gc_content))

        # After searching, run primers found in this window through the primer scorer
        if len(valid_primers) > 0:
            primer_res.append(score_primer_set([primer[0] for primer in valid_primers], target_seqs))

        # Update current window start
        current_window_start += step_size
    
    # Return primer results, should only have ~len(search_seq) / step_size total in here
    return primer_res

# Determines average length of FASTAs
def amplicon_length_isoforms(primers, fastas):
    total_len = 0
    min_len = max([len(fasta.seq) for fasta in fastas])
    max_len = 0
    num_hits = 0

    # Loop through all isoforms
    for fasta in fastas:
        # Calculate end/start
        fw_start = str(fasta.seq).upper().find(primers[0].upper())
        rev_end = str(fasta.seq).upper().find(str(Seq(primers[1]).reverse_complement()).upper()) + len(primers[1])

        if fw_start != -1 and rev_end != -1:
            num_hits += 1
            amplicon_len = rev_end - fw_start

            # Add to total length
            total_len += amplicon_len

            if (amplicon_len > max_len): max_len = amplicon_len

            if (amplicon_len < min_len): min_len = amplicon_len 

    # Returns average length of amplicon, min length, max length
    return [round(float(total_len) / num_hits), min_len, max_len]

def score_primer_set(
    primer_set,
    target_seqs,
    positive_strand=True,
):
    """
    This function looks across windows of some size and steps some amount to find the optimal primer
    in that window.

    This'll be better for runtime because we don't want to keep searching for unique forward + primer
    combinations for every single possible point-point pair. It would be far better
    """
    test_primer = primer_set[0] if positive_strand else str(Seq(primer_set[0]).reverse_complement())
    gb_first = [seq for seq in target_seqs if seq.seq.find(test_primer) != -1][0]

    def normalize_arr(data):
        if np.max(data) - np.min(data) == 0: return np.array(data)

        return np.array((data - np.min(data)) / (np.max(data) - np.min(data)))

    # Gets the RNA base based on a given sequence and the primer sequence
    def get_rna_base(full_sequence, primer):
        primer_ind = full_sequence.find(primer)
        rna_seq = full_sequence
        if primer_ind == -1:
            rna_seq = [seq.seq for seq in target_seqs if seq.seq.find(primer) != -1][0]
            primer_ind = rna_seq.find(primer)

        # Returns the base immediately after the end of the primer
        return rna_seq[primer_ind + len(primer)]

    # Literally just checks the ends of the primers to see if an rU base is required
    # Scores are either 0, 0.5, or 1. (0.5 points subtracted if a primer requires an rU)  
    def rna_base_optimizer(arr):
        # Find the first genbank file that contains 
        # gb_first = target_seqs[0]

        # Sequence from Genbank file, we set based on primer direction
        search_seq = gb_first.seq if positive_strand else gb_first.seq.reverse_complement()

        # Check forward
        rna_bases = [get_rna_base(search_seq, primer) for primer in arr]
        base_scores = [0 if base == "T" else 0.5 for base in rna_bases]

        return normalize_arr(base_scores)
    
    # Scores mono/double repeats
    def repeat_scorer(primer):
        # Things that can be repeated
        repeated_strings = ["A", "T", "C", "G", 
                            "AT", "TA", "CG", "GC", "AC", "CA", "GT", "TG"]
        # Loop through types of repeats
        for repeat_type in repeated_strings:
            repeat_present = primer.find(repeat_type * 4)
            if repeat_present != -1:
                return 0

        return 1

    # Checks for extended repeats in primers, both mono-/di-nucleotide
    def repeat_optimizer(arr):
        # All repeats
        repeats = np.array([repeat_scorer(primer) for primer in arr])

        return normalize_arr(repeats)
    
    # Score GC content
    def gc_scorer(primer):
        # Measure GC content
        g_content = str(primer).upper().count("G")
        c_content = str(primer).upper().count("C")
        gc_content = (g_content + c_content) / len(primer)

        # Score: if <0.4 or >0.6, set to 0
        score = 0 if (gc_content > 0.6 or gc_content < 0.4) else 1 - 10 * np.abs(gc_content - 0.5)

        return score

    # Optimizes for GC content: aim for ~50%
    def gc_optimizer(arr):
        # Calculate forward/reverse gc contents
        gcs = np.array([gc_scorer(primer) for primer in arr])

        return normalize_arr(gcs)

    # Calculate scores for each primer pair
    base_scores = rna_base_optimizer(primer_set)
    repeat_scores = repeat_optimizer(primer_set)
    gc_scores = gc_optimizer(primer_set)

    """
    Check this... incrementing?
    """
    weighted_scores = + 0.2 * base_scores + 0.2 * repeat_scores + 0.6 * gc_scores

    # Return primers and weighted score
    return (primer_set[np.argmax(weighted_scores)], np.max(weighted_scores))

def design_specificity(crRNA, fw_pool, rev_pool, non_target_seqs):
    """
    Checks specificity of primer set + crRNA using non-target sequences.

    I imagine the more efficient way to do this is essentially what I do in design_candidates.

    Search for each primer first, then create a matrix. So instead of constantly searching for sequences
    in each combination of forward, reverse, and crRNA, we just search once and then use vector operations
    to calculate specificity.
    """
    # Search for forward primer across all non-target sequences
    # fw_specificity = [seq.seq.find(fw_primer) for seq in non_target_seqs]
    # rev_specificity = [seq.seq.find(rev_primer) for seq in non_target_seqs]
    cr_specificity = [seq.seq.find(crRNA) for seq in non_target_seqs]

    # Create combinations of forward and reverse primers from pools
    fw_rev_combos = list(itertools.product(fw_pool, rev_pool))

    for combo in fw_rev_combos:
        fw_specificity = [seq.seq.find(combo[0]) for seq in non_target_seqs]
        rev_specificity = [seq.seq.find(str(Seq(combo[1]).reverse_complement())) for seq in non_target_seqs]

        # Loop through indices of non-target sequences
        for ind in range(len(non_target_seqs)):
            # If any of the sequences have the combination of primer + crRNA across all, return False
            if fw_specificity[ind] != -1 and rev_specificity[ind] != -1 and cr_specificity[ind] != -1: return False

    # Otherwise, return true
    return True

def optimize_primer_length(
    fw_set, 
    rev_set,
    target_seqs,
    desired_length=300,
):
    """
    This function takes a set of forward primers and a set of reverse primers and returns
    a pair of primers that yields the most optimal length amplicon.
    """
    # Make all possible primer pair combinations
    primer_combos = list(itertools.product(fw_set, rev_set))

    # Loop through all primer combos and get the amplicon length
    amplicon_lengths = [amplicon_length_isoforms(primer_pair, target_seqs)[0] for primer_pair in primer_combos]

    """
    Update this code later where we don't just return the primer pair with the amplicon length
    closest to our desired, let's also incorporate scoring with things like GC content, rhPCR bases, etc.
    """
    # Get the index of the primer pair with the closest amplicon length to the desired length
    closest_index = amplicon_lengths.index(min(amplicon_lengths, key=lambda x: abs(x - desired_length)))

    # Return the ideal primer combo and its length
    return (primer_combos[closest_index], amplicon_lengths[closest_index])

def get_exon_from_seq(search_seq, ref_seq, exon_map):
    """
    Provided any sequence, a reference sequence, and exon map, returns the exon that the sequence is in.

    The reference sequence should contain the sequence we're searching for.
    """
    # Get location of sequence
    seq_loc =  ref_seq.seq.find(search_seq)
    seq_len = len(search_seq)

    # Get exons of ref_seq
    ref_exons = [exon for exon in ref_seq.features if exon.type == "exon"]

    # Get exon signature of ref_seq
    ref_exon_sig = get_exon_signature(ref_seq, exon_map)

    # Loop through exons
    for ind, exon in enumerate(ref_exons):
        # Get exon start and end
        e_start = int(exon.location.start)
        e_end = int(exon.location.end)

        start_in_exon = seq_loc >= e_start and seq_loc <= e_end
        end_in_exon = seq_loc + seq_len >= e_start and seq_loc + seq_len <= e_end

        if start_in_exon and end_in_exon: return ref_exon_sig[ind]
        if start_in_exon and not end_in_exon: return (ref_exon_sig[ind], ref_exon_sig[ind + 1])

    return None

def search_primer_sets(
    # Matrix of exon and exon junction locations across isoforms
    loc_matrix,
    # Candidates for primer locations
    loc_candidates,
    # Indices of target sequences
    target_indices,
    # Points shared across target sequences
    shared_points,
    # Maximum number of primers in subset, set to 5 by default
    max_primers=5,
):
    """
    We want to start by creating a matrix of possible primer locations (exon or exon junction, 
    these are columns) across isoforms (rows). A value of 1 in the matrix represents an available
    exon or junction for that isoform.

    We'll create this matrix for all isoforms using the exon map.
    """
    # Create target matrix: keep rows that match a 1 in target_indices
    target_matrix = loc_matrix[[ind for ind, target_ind in enumerate(target_indices) if target_ind == 1], :]
    flipped_target_matrix = np.logical_not(target_matrix).astype(int)

    # Get non-target matrix
    non_target_matrix = loc_matrix[[ind for ind, target_ind in enumerate(target_indices) if target_ind == 0], :]
    print('non_target',non_target_matrix)
    # Returns nCr for n and r
    def ncr(n, r):
        return int(np.math.factorial(n) / (np.math.factorial(r) * np.math.factorial(n - r)))
    
    # Sees if a given primer set is valid
    def test_coverage(loc_inds):
        # Invert 0s and 1s of target matrix and multiply across rows
        # If sum of resulting vector is 0, then all targets are covered
        targets_covered = np.sum(np.prod(flipped_target_matrix[:, loc_inds], axis=1)) == 0

        # Next we want to check specificity
        # Don't invert this time, just multiply across rows and add up
        # If sum of resulting vector is equal to the length of the vector, then primer set is specific
        # primers_specific = np.sum(np.prod(non_target_matrix[:, loc_inds], axis=1)) == non_target_matrix.shape[0]

        # Return both values, if the second is False, then we can mark this as an invalid parent for future levels
        # return (targets_covered, primers_specific)
        return targets_covered
    
    # Inputs are all primer_loc indices, fw_pool and rev_pool are pools of primers (both lists)
    def test_detection_specificity(crRNA, fw_pool, rev_pool):
        # Get columns and reformat accordingly for tensor operations
        crRNA_cols = non_target_matrix[:, crRNA][:, np.newaxis, np.newaxis]
        fw_cols = non_target_matrix[:, fw_pool][:, :, np.newaxis]
        rev_cols = non_target_matrix[:, rev_pool][:, np.newaxis, :]

        # Just check non-target matrix to see if we get 1 * 1 * 1 = 1 for any rows
        specificity = np.sum(crRNA_cols * fw_cols * rev_cols) == 0

        return specificity

    # Given an exon or exon junction, returns all indices of exons or exon junctions that flank to the left or right
    def get_flanking_candidates(restriction_point, left=True):
        # Is this an exon or exon junction?
        is_exon = isinstance(restriction_point, int)
        comp_point = restriction_point if is_exon else (restriction_point[0] + 1 if left else restriction_point[1] - 1)

        flanking_inds = []
        # for ind, loc in enumerate(loc_candidates):
        #     loc_is_exon = isinstance(loc, int)
        #     # Exon and left is True
        #     if (loc_is_exon) and left:
        #         if loc < comp_point: flanking_inds.append(ind)
        #     elif (loc_is_exon) and not left:
        #         if loc > comp_point: flanking_inds.append(ind)
        #     elif (not loc_is_exon) and left:
        #         if loc[1] < comp_point: flanking_inds.append(ind)
        #     elif (not loc_is_exon) and not left:
        #         if loc[0] > comp_point: flanking_inds.append(ind)

        for ind, loc in enumerate(loc_candidates):
            loc_is_exon = isinstance(loc, int)
            # Exon and left is True
            if (loc_is_exon) and left:
                if loc <= comp_point: flanking_inds.append(ind)
            elif (loc_is_exon) and not left:
                if loc >= comp_point: flanking_inds.append(ind)
            elif (not loc_is_exon) and left:
                if loc[1] <= comp_point: flanking_inds.append(ind)
            elif (not loc_is_exon) and not left:
                if loc[0] >= comp_point: flanking_inds.append(ind)
        
        return flanking_inds
    
    # Converts an index into a combination of list indices
    def index_to_combination(index, n, r):
        if r > n or index >= math.comb(n, r):
            raise ValueError("Invalid input values")

        combination = []
        current = 0

        while r > 0:
            if index < math.comb(n - current - 1, r - 1):
                combination.append(current)
                r -= 1
            else:
                index -= math.comb(n - current - 1, r - 1)
            current += 1

        return combination

    valid_results = []
    # Loop through valid primer combination lengths, start with 2
    for num_primers in range(2, max_primers + 1):
        level_results = []

        # Loop through shared points
        for shared_point in shared_points:
            # Lists of lists of tuples (left_inds, right_inds)
            combos_for_shared_point = []

            # Get all primers to the left and right of the crRNA exon/junction
            left_candidates = get_flanking_candidates(shared_point, left=True)
            right_candidates = get_flanking_candidates(shared_point, left=False)

            # print(shared_point, [loc_candidates[ind] for ind in left_candidates], [loc_candidates[ind] for ind in right_candidates])

            if len(left_candidates) == 0 or len(right_candidates) == 0: continue

            # Loop through different left_num + right_num = num_primers combinations
            for num_left in range(1, num_primers):
                # Number of right primers to search
                num_right = num_primers - num_left

                # If we're searching a num of left primers greater than candidates, break
                if num_left > len(left_candidates): break
                # Since right nums decrease, just continue if greater
                if num_right > len(right_candidates): continue

                # Start with all left combinations
                for ind in range(ncr(len(left_candidates), num_left)):
                    # Create a list of indices for the left primers
                    left_target_inds = index_to_combination(ind, len(left_candidates), num_left)

                    left_test = [left_candidates[ind] for ind in left_target_inds]

                    # Test validity of this current primer set
                    forward_coverage = test_coverage(left_test)

                    # If this pool of forward primers doesn't work, skip the reverse primer search
                    if not forward_coverage: continue

                    # If we have coverage, search for reverse primers
                    for ind in range(ncr(len(right_candidates), num_right)):
                        # Create a list of indices for the right primers
                        right_target_inds = index_to_combination(ind, len(right_candidates), num_right)

                        right_test = [right_candidates[ind] for ind in right_target_inds]

                        # Test reverse coverage
                        reverse_coverage = test_coverage(right_test)

                        # If no coverage, skip
                        if not reverse_coverage: continue

                        # Otherwise, go ahead and test specificity across pools
                        specificity = test_detection_specificity(loc_candidates.index(shared_point), left_test, right_test)
                        # if specificity: 
                        #     print("Testing specificity for:", shared_point, loc_candidates.index(shared_point), left_test, right_test, specificity)

                        # If we have specificity, add to valid results
                        if specificity: combos_for_shared_point.append((left_test, right_test))

            # if len(all_left) > 0 and len(all_right) > 0:
            if len(combos_for_shared_point) > 0:
                level_results.append((shared_point, combos_for_shared_point))

        if len(level_results) > 0:
            valid_results.append(level_results)
            return level_results

    # Return valid results here
    return None

def get_primer_loc_matrix(target_seqs, all_seqs, exon_map):
    """
    Creates primer location matrix with rows as isoforms, columns as possible primer locations
    """
    # Get exon signatures for all sequences
    all_signatures = [get_exon_signature(seq, exon_map) for seq in all_seqs]

    # Add junctions to signatures
    sig_w_junctions = [insert_junctions(signature) for signature in all_signatures]

    # Get target signatures indices by indexing all_signatures
    target_indices = get_target_indices(target_seqs, all_seqs)

    # Loop through target signatures and create a set of possible primer locations
    primer_locs = []
    for ind, target_bool in enumerate(target_indices):
        # Skip non-targets
        if target_bool == 0: continue

        # Add primer locations for target
        primer_locs += sig_w_junctions[ind]
        primer_locs = list(set(primer_locs))
        
    # print(primer_locs)

    # Create location matrix, order of primer locs doesn't really matter, so long as it maps to primer_locs
    loc_matrix = []
    for sig in sig_w_junctions:
        loc_matrix.append([1 if loc in sig else 0 for loc in primer_locs])

    return (np.array(loc_matrix), primer_locs)


# Returns a list of 0s and 1s representing target or non-target
def get_target_indices(targets, all):
    # Get target signatures indices by indexing all_signatures
    indices = []
    for ind, seq in enumerate([all_seq.seq for all_seq in all]):
        if seq in [target_seq.seq for target_seq in targets]: indices.append(1)
        else: indices.append(0)

    return indices

def get_shared_points(target_seqs, all_seqs, exon_map):
    """
    Inputs are target gb records in a list, list of all gb records for that gene, as well as an exon map.

    Function uses get_exon_signature to get signatures for all isoforms in all_seqs.

    Function then identifies exons that are shared across target_seqs. These exons or exon junctions are
    candidates for crRNAs.

    We can then loop through these exons and look for primer sets that are either within the exons or flank
    the exons.

    Primer + crRNA combinations should not detect any non-target sequences. However, 

    Returns a list of tuples of the form (point_1, point_2, point_3). Points can be either an exon or an
    exon-exon junction. point_1, point_2, and point_3 can be the same exon or they can be different, but they should
    be in order.
    """
    # Get exon signatures for all sequences
    all_signatures = [get_exon_signature(seq, exon_map) for seq in all_seqs]

    # Get target signatures indices by indexing all_signatures
    target_indices = []
    for ind, seq in enumerate([all_seq.seq for all_seq in all_seqs]):
        if seq in [target_seq.seq for target_seq in target_seqs]: target_indices.append(ind)

    # Get signatures
    target_signatures = [all_signatures[ind] for ind in target_indices]

    # Start by assuming that all exons are shared
    exon_indices = list(range(len(exon_map)))
    shared_exons = exon_indices.copy()

    # Loop through all exons in exon_map and check if they are shared across all target_seqs
    for exon in exon_indices:
        # Loop through all target signatures and check if exon is present, get target exon signature from all_signatures
        for target_signature in target_signatures:
            # If exon is not present in target signature, remove it from shared exons
            if exon not in target_signature and exon in shared_exons:
                shared_exons.remove(exon)
                break

    # Initialize junctions
    ee_junctions = []

    # Loop through all target_signatures and check for exon-exon junctions
    for target_signature in target_signatures:
        # Loop through all exon-exon junctions in target_signature
        for ind in range(len(target_signature) - 1):
            # Get exon-exon junction
            junction = (target_signature[ind], target_signature[ind + 1])

            # If junction is not in ee_junctions, add it
            if junction not in ee_junctions: ee_junctions.append(junction)

    # Loop through ee_junctions and remove the ones that aren't shared across all target_signatures
    shared_junctions = ee_junctions.copy()
    for junction in ee_junctions:
        # Loop through all target_signatures and check if junction is present
        for target_signature in target_signatures:
            # Get junctions for that target_signature
            target_junctions = [(target_signature[ind], target_signature[ind + 1]) for ind in range(len(target_signature) - 1)]

            # If junction is not present, remove it
            if junction not in target_junctions and junction in shared_junctions: shared_junctions.remove(junction)

    # Combine shared_exons and shared_junctions into a single list in order
    shared_points = shared_exons.copy()
    for junction in shared_junctions:
        # Figure out where to insert junction into shared_points
        for ind, point in enumerate(shared_points):
            # If junction is to the right of point, insert it
            if junction[0] == point:
                shared_points.insert(ind + 1, junction)
                break

    return shared_points

def insert_junctions(exon_signature):
    """
    Takes an exon signature as input and inserts exon junctions in order.
    """
    # Initialize junctions
    ee_junctions = []

    # Loop through all exon-exon junctions in target_signature
    for ind in range(len(exon_signature) - 1):
        # Get exon-exon junction
        junction = (exon_signature[ind], exon_signature[ind + 1])

        # If junction is not in ee_junctions, add it
        if junction not in ee_junctions: ee_junctions.append(junction)

    # Combine shared_exons and shared_junctions into a single list in order
    shared_points = exon_signature.copy()
    for junction in ee_junctions:
        # Figure out where to insert junction into shared_points
        for ind, point in enumerate(shared_points):
            # If junction is to the right of point, insert it
            if junction[0] == point:
                shared_points.insert(ind + 1, junction)
                break

    return shared_points

def get_exon_signature(gb_record, exon_map):
    """
    This function should pass in a Genbank record as well as an exon "map" for that gene.

    It should return a list of exons and exon junctions and this list should consist of IDs that correspond
    to the exons and their order in the exon map.
    """
    # Get exon sequences from gb_record
    gb_exons = [gb_record.seq[exon[0]:exon[1]] for exon in get_exons(gb_record)]

    return [exon_map.index(exon) for exon in gb_exons]

def create_exon_map(isoforms):
    """
    Creating a map of exons. The idea is that we want to start with the largest transcript in the set of
    isoforms and start creating a list of exons. We then want to iterate through the other isoforms and
    look at the exons in each of those. If an exon is already present, we skip over it, if we come across
    a new exon, we add it to the list based on its position relative to flanking exons.

    It gets a little tricky with the ordering because we might have a case where we know that an exon is to
    the left of another exon. But in another we isoform, we might find another exon that's to the left of
    a known flanking exon, but we don't know if it's to the left or the right of the last one that we just inserted
    into the list. For example, if we use integers from 0 onwards to represent unique exons, we might have the case:

    [0, 1, 4, 5, 6] (first iteration, found these exons)

    -> 0-1-3-4-5-6 (looking at next isoform that has these exons)

    [0, 1, 3, 4, 5, 6] (adding 3 to the list, we know its next to 4 and to the right of 1)

    -> 0-1-2-4-5-6 (looking at next isoform that has these exons)

    [0, 1, 3, 2, 4, 5, 6] (adding 2 to the list, we know its to the left of 4 and generally to the right of 1, but not sure if it's to the left or right of 3)

    -> 0-1-2-3-4-5 (looking at this isoform, and now we know that 2 is to the left of 3)
    """
    # Make a copy of isoforms
    isoforms = isoforms.copy()

    # Iterate through isoforms and find the transcript with the most exons
    longest_id = -1
    longest_length = -1
    for ind, isoform in enumerate(isoforms):
        isoform_exons = get_exons(isoform)
        if len(isoform_exons) > longest_length:
            longest_id = ind
            longest_length = len(isoform_exons)

    # Start by adding exon sequences from the longest isoform
    longest_exon_indices = get_exons(isoforms[longest_id])
    exons = [isoforms[longest_id].seq[exon[0]:exon[1]] for exon in longest_exon_indices]

    # Remove the longest isoform from list
    isoforms.pop(longest_id)

    # Then loop through the rest of the isoforms and start adding unique exons in order
    for is_ind, isoform in enumerate(isoforms):
        # print("Currently looking at isoform:", isoform.description)
        isoform_exons = [isoform.seq[exon[0]:exon[1]] for exon in get_exons(isoform)]

        # Update the order based on isoform_exons, because the order there is definitively true
        # checked_map = id_sort(exons, isoform_exons)

        # Iterate through exons in isoform
        for ex_ind, exon in enumerate(isoform_exons):
            # print("Looking at:", ex_ind, end=", ")
            # If exon is already in list, skip it
            if exon in exons: continue

            # Otherwise, find where it should be inserted
            insertion_index = find_insertion_index(exons, exon, isoform_exons)
            # print("Calculated insertion index:", insertion_index)

            # Insert exon into list
            exons.insert(insertion_index, exon)

        # print(exons)
        # print(isoform_exons)
        exons = id_sort(exons, isoform_exons)
        # print(exons)

    return exons

"""
Parameters: current_map (current exon map), exon_seq (sequence of exon being inserted), 
all_exons (all exons for the isoform being analyzed)
"""
def find_insertion_index(current_list, sequence, all_sequences):
    # Get the index of the sequence in all_sequences
    sequence_index = all_sequences.index(sequence)

    # Find the position where the sequence should be inserted in current_list
    insertion_index = len(current_list)  # Initialize to the end of the list
    for i, seq in enumerate(current_list):
        if seq in all_sequences:
            if all_sequences.index(seq) > sequence_index:
                insertion_index = i
                break
        else:
            continue

    return insertion_index

"""
Given array_1, which contains an accumulation of semi-ordered exons, and array_2, which contains
a set of definitely ordered exons, we want to sort array_1 such that the semi-ordered exons are
in the same order as the definitely ordered exons.
"""
def id_sort(arr_1, arr_2):
    array_1 = list(range(len(arr_1)))
    array_2 = [arr_1.index(val) for val in arr_2]

    positions = {id_: index for index, id_ in enumerate(array_2)}
    # print(positions)
    sorted_array_1 = []
    not_in_array_2 = []

    for id_ in array_1:
        if id_ not in positions:
            not_in_array_2.append(id_)
            continue

        insert_index = len(sorted_array_1)
        while insert_index > 0 and positions[sorted_array_1[insert_index - 1]] > positions[id_]:
            insert_index -= 1
        sorted_array_1.insert(insert_index, id_)

    for id_ in not_in_array_2:
        insert_index = len(sorted_array_1)
        while insert_index > 0 and sorted_array_1[insert_index - 1] > id_:
            insert_index -= 1
        sorted_array_1.insert(insert_index, id_)

    return [arr_1[id_] for id_ in sorted_array_1]

"""
Returns a list of all exons in a given Genbank file by index.
"""
def get_exons(sequence_gb):
    return [(int(feature.location.start), int(feature.location.end)) for feature in sequence_gb.features if feature.type == "exon"]
