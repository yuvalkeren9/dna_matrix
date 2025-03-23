import random
from utils import DNA_BASES, shift_string, get_next_base


def consume_logic_random(dna_list, base):
    """
    A consume logic function, which picks a strand at random to consume a base from.

    :param dna_list: A list of dna strands.
    :param base: The base to consume.
    :return: The strand to consume
    """
    assert dna_list

    # Filter the strands that begin with the given base
    filtered_strands = [temp_strand for temp_strand in dna_list if temp_strand.startswith(base)]

    if len(filtered_strands) == 0:
        return None

    return random.choice(filtered_strands)


# This version just looks if a selected strand has the next base in its lookahead
def consume_logic_lookahead_one_v1(dna_list, base):
    assert dna_list

    # Filter the strands that begin with the given base
    filtered_strands = [temp_strand for temp_strand in dna_list if temp_strand.startswith(base)]

    if len(filtered_strands) == 0:
        return None

    next_base = get_next_base(base)

    # Look for a strand with the next base as a second character
    for strand in filtered_strands:
        if len(strand) > 1 and strand[1] == next_base:
            return strand

    # Did not find any strand with the next base as the next character, returning a random strand
    return random.choice(filtered_strands)


# Can be written with better complexity
def consume_logic_lookahead_one_v2(dna_list, base):
    assert dna_list

    # Filter the strands that begin with the given base
    filtered_strands = [temp_strand for temp_strand in dna_list if temp_strand.startswith(base)]

    if len(filtered_strands) == 0:
        return None

    next_base = get_next_base(base)
    shifted_dna_bases = shift_string(DNA_BASES, next_base)

    is_base_needed = []
    for curr_base in shifted_dna_bases[:-1]: # Last char will have special treatment
        is_base_needed.append(not any(s[0] == curr_base for s in dna_list))
    temp_counter = 0
    for strand in dna_list:
        if strand[0] == shifted_dna_bases[-1]:
            temp_counter += 1
            if temp_counter == 2:
                break
    is_base_needed.append(temp_counter > 1)

    for i, curr_base in enumerate(shifted_dna_bases):
        for strand in filtered_strands:
            if len(strand) > 1 and strand[1] == curr_base and (is_base_needed[i] is True):
                return strand

    # Did not find any strand with the next base as the next character, returning a random strand
    return random.choice(filtered_strands)

def consume_logic_lookahead_one_v3(dna_list, base):
    assert dna_list

    # Filter the strands that begin with the given base
    filtered_strands = [temp_strand for temp_strand in dna_list if temp_strand.startswith(base)]

    if len(filtered_strands) == 0:
        return None

    bases_counter_d = {}
    for curr_base in DNA_BASES:
        bases_counter_d[curr_base] = 0

    for strand in dna_list:
        bases_counter_d[strand[0]] += 1

    # Removing 1 from current base, since it will be deleted by one anyway
    bases_counter_d[base] -= 1

    # Edge case: last base in last strand
    max_key = max(bases_counter_d, key=bases_counter_d.get)
    if max_key == 0:
        print("edge case!")
        return random.choice(filtered_strands)

    # Define a custom order
    custom_order = shift_string(DNA_BASES, get_next_base(base))
    # Create a mapping from base to its priority in custom_order
    order_rank = {base: rank for rank, base in enumerate(custom_order)}

    # Pick a strand that will unlock a base that is least common in current synthesizeable (funny word) indices
    sorted_by_value = sorted(bases_counter_d.items(), key=lambda x: (x[1], order_rank[x[0]]))
    for curr_base, _ in sorted_by_value:
        for strand in filtered_strands:
            if len(strand) > 1 and strand[1] == curr_base:
                return strand

    # Did not find any strand with the next base as the next character, returning a random strand
    return random.choice(filtered_strands)


# Can be written with better complexity
# Same as V2, but instead of picking randomly the next base, pick the longest strand (with the base)
def consume_logic_lookahead_one_v4(dna_list, base):
    assert dna_list

    # Filter the strands that begin with the given base
    filtered_strands = [temp_strand for temp_strand in dna_list if temp_strand.startswith(base)]

    if len(filtered_strands) == 0:
        return None

    next_base = get_next_base(base)
    shifted_dna_bases = shift_string(DNA_BASES, next_base)

    is_base_needed = []
    for curr_base in shifted_dna_bases[:-1]: # Last char will have special treatment
        is_base_needed.append(not any(s[0] == curr_base for s in dna_list))
    temp_counter = 0
    for strand in dna_list:
        if strand[0] == shifted_dna_bases[-1]:
            temp_counter += 1
            if temp_counter == 2:
                break
    is_base_needed.append(temp_counter > 1)

    curr_longest_strand = ""
    curr_max_len = 0
    for i, curr_base in enumerate(shifted_dna_bases):
        for strand in filtered_strands:
            if len(strand) > 1 and strand[1] == curr_base and (is_base_needed[i] is True):
                if len(strand) > curr_max_len:
                    curr_longest_strand = strand
                    curr_max_len = len(strand)

    if curr_longest_strand != "":
        return curr_longest_strand

    # Did not find any strand with the next base as the next character, returning a random strand
    return random.choice(filtered_strands)