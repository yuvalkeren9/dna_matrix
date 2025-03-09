# DNA project

# Imports
import random
import sys
import copy
import numpy as np
from direct_solvers import get_optimal_row_solution

# Globals:
DNA_BASES = 'ACGT'
DEBUG = False

#trying sphnix
def my_function(param1, param2):
    """
    This is an example function.

    :param param1: The first parameter.
    :param param2: The second parameter.
    :return: A string describing the operation.
    """
    return f"Result of {param1} and {param2}"

# Functions to set up experiments
def create_random_strand(length):
    """
    Generates a random DNA strand of the given length.

    Parameters:
    :param length:  The length of the DNA strand.
    :return: A string describing the DNA strand.
    """
    return ''.join(random.choice("ACGT") for _ in range(length))


def create_dna_list(num_of_strands, length_of_strands):
    """
    Generates a list of DNA strands, to stimulate a row in the synthesis machine.

    :param num_of_strands: The number of strands.
    :param length_of_strands: The length of the strands.
    :return: A list of strings representing the DNA strands consisting of 'A', 'C', 'G', 'T'.
    """
    return [create_random_strand(length_of_strands) for _ in range(num_of_strands)]


# Functions for consuming
def consume_base(index, dna_list):
    """
    Remove a base from a dna strand, dna_list[index]

    :param index: The strand to pick from dna_list.
    :param dna_list: The list of the dna strands.
    """
    if 0 <= index <= len(dna_list):
        dna_list[index] = dna_list[index][1:]
        if dna_list[index] == "":
            dna_list.pop(index)
    else:
        print("index is out of range!")


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
        #        print("Debug. No base available!") TODO: remove
        return None

    return random.choice(filtered_strands)


# This version just looks if a selected strand has the next base in its lookahead
def consume_logic_lookahead_one_v1(dna_list, base):
    assert dna_list

    # Filter the strands that begin with the given base
    filtered_strands = [temp_strand for temp_strand in dna_list if temp_strand.startswith(base)]

    if len(filtered_strands) == 0:
        #        print("Debug. No base available!") TODO: remove
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
        #        print("Debug. No base available!") TODO: remove
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
        #        print("Debug. No base available!") TODO: remove
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
        #        print("Debug. No base available!") TODO: remove
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


def consume_with_logic(dna_list, base, func):

    assert dna_list

    selected_strand = func(dna_list, base)

    if selected_strand is None:
        return

    # Find the index of the selected strand in the original list
    selected_index = dna_list.index(selected_strand) if selected_strand else None

    consume_base(selected_index, dna_list)

    return


# Functions for DNA arithmetics

def get_next_base(current_base):
    index = DNA_BASES.find(current_base)

    if index == -1:
        raise ValueError(f"Base '{current_base}' is not a valid DNA base.")

    # Return the next base in the string, wrapping around cyclically
    next_index = (index + 1) % len(DNA_BASES)
    return DNA_BASES[next_index]


def shift_string(s, char):
    """Shifts string `s` to begin with the first occurrence of `char`."""
    index = s.find(char)

    if index == -1:
        raise ValueError(f"Character '{char}' not found in string.")

    return s[index:] + s[:index]  # Rotate the string


# DNA synthesiser machine TODO: currently only a single row is supported
def synthezise(dna_list, func, should_visualize = False):
    base = DNA_BASES[0]
    num_of_cycles = 0
    while dna_list:
        if should_visualize:
            clear_console()
            for strand in dna_list:
                print(strand)
            time.sleep(0.1)
        consume_with_logic(dna_list, base, func)

        num_of_cycles += 1
        base = get_next_base(base)

    return num_of_cycles


def run_experiment(num_of_strands, len_of_strands, logic_func):
    cycle_results = []
    percent_results = []
    for _ in range(300):
        cycles, percent = run_single_experiment(num_of_strands, len_of_strands, logic_func)
        cycle_results.append(cycles)
        percent_results.append(percent)

    mean = np.mean(cycle_results)
    percent_mean = np.mean(percent_results)
    print(f"The mean is: {mean}")
    print(f"The percent is is: {percent_mean}")


def run_single_experiment(num_of_strands, len_of_strands, logic_func, own_list=None):
    if own_list is None:
        dna_list = create_dna_list(num_of_strands, len_of_strands)
    else:
        dna_list = own_list

    dna_list_copy = copy.deepcopy(dna_list)

    num_of_cycles = synthezise(dna_list, logic_func, should_visualize=False)
    # print(f"The num of cycles for selected logic: {num_of_cycles}")

    num_of_cycles_control_group = synthezise(dna_list_copy, consume_logic_random, should_visualize=False)
    # print(f"The num of cycles for random logic: {num_of_cycles_control_group}")

    cycles_saved = num_of_cycles_control_group - num_of_cycles
    percentage_save = num_of_cycles_control_group / num_of_cycles

    return cycles_saved, percentage_save

import os
import sys
import time

def clear_console():
    # Works for both Windows and Linux/Mac
    os.system('cls' if os.name == 'nt' else 'clear')

if __name__ == '__main__':
    run_experiment(3, 100, consume_logic_lookahead_one_v1)

    #run_single_experiment(5, 20, consume_logic_lookahead_one_v3)

