"""
Rename this to mod_lib.
"""

def design_gen1(sequence, genbank_files, direction):
    # Just grab first genbank file and use as template
    gb_first = genbank_files[0]
    # Sequence from Genbank file, we set based on primer direction
    strand_sequence = gb_first.seq if direction > 0 else gb_first.seq.reverse_complement()

    add_on = ""
    fw_primer_ind = strand_sequence.index(sequence)
    # Get the RNA residue
    rna_residue = strand_sequence[fw_primer_ind + len(sequence)]
    add_on += "r" + ("U" if rna_residue == "T" else rna_residue)
    # Then add on DNA residue
    add_on += strand_sequence[fw_primer_ind + len(sequence) + 1]
    add_on += strand_sequence[fw_primer_ind + len(sequence) + 2]
    add_on += strand_sequence[fw_primer_ind + len(sequence) + 3]
    add_on += strand_sequence[fw_primer_ind + len(sequence) + 4]
    # Add a mismatch
    add_on += calculate_mismatch(strand_sequence[fw_primer_ind + len(sequence) + 5])
    # Add one blocker
    add_on += "/3SpC3/"
        
    return sequence + add_on

def design_gen2(sequence, genbank_files, direction):
    # Just grab first genbank file and use as template
    gb_first = genbank_files[0]
    # Sequence from Genbank file, we set based on primer direction
    strand_sequence = gb_first.seq if direction > 0 else gb_first.seq.reverse_complement()

    add_on = ""
    fw_primer_ind = strand_sequence.index(sequence)
    # Get the RNA residue
    rna_residue = strand_sequence[fw_primer_ind + len(sequence)]
    add_on += "r" + ("U" if rna_residue == "T" else rna_residue)
    # Then add on DNA residue
    add_on += strand_sequence[fw_primer_ind + len(sequence) + 1]
    # Add the two blockers
    add_on += "/iSpC3/" + "/iSpC3/"
    # Add another DNA residue
    add_on += strand_sequence[fw_primer_ind + len(sequence) + 4]
    # Finally, add a mismatch
    add_on += calculate_mismatch(strand_sequence[fw_primer_ind + len(sequence) + 5])
        
    return sequence + add_on

# Calculates a DNA mismatch for a given base
def calculate_mismatch(base):
    if base.upper() == "A":
        return "G"
    if base.upper() == "T":
        return "C"
    if base.upper() == "C":
        return "A"
    if base.upper() == "G":
        return "T"
    else:
        return "A"

def add_t7_promoter(sequence):
    # Most optimal T7 promoter sequence
    promoter_sequence = "TAATACGACTCACTATAGGG" + "AAATA"

    return promoter_sequence + sequence